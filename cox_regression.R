library(tidyverse)
library(sas7bdat)

df <- read.csv("NHANES2.csv", stringsAsFactors = FALSE)

# ----1. Clean data and create variables --------------------

#create a function to better cleaning the data
to_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))


df_clean <- df %>%
  # death
  mutate(DEATH_clean = trimws(DEATH)) %>%
  filter(DEATH_clean %in% c("0", "1")) %>%
  mutate(DEATH_num = as.integer(DEATH_clean)) %>%
  
  # calculate age at exit for our time scale
  mutate(
    BORN_YR_n = to_num(BORN_YR),
    LAST_YR_n = to_num(LAST_YR),
    DIE_YR_n  = to_num(DIE_YR),
    age_exit  = if_else(
      DEATH_num == 1,
      (DIE_YR_n  + 1900) - (BORN_YR_n + 1900),
      (LAST_YR_n + 1900) - (BORN_YR_n + 1900)
    )
  ) %>%
  
  # exclude age_exit <= AGEYRS
  filter(age_exit > AGEYRS) %>%
  
  # BMI
  mutate(BMI = WT / (HEIGHT / 100)^2) %>%
  
  # NIACIN quartiles
  mutate(NIACIN_q = cut(
    NIACIN,
    breaks = c(-Inf, quantile(NIACIN, probs = c(0.25, 0.50, 0.75)), Inf),
    labels = c("Q1", "Q2", "Q3", "Q4"),
    right  = TRUE
  )) %>%
  
  # Recode covariates
  mutate(
    Sex          = factor(SEX, levels = c(1, 2),
                          labels = c("Male", "Female")),
    Race         = factor(RACE, levels = c(1, 2, 3),
                          labels = c("White", "Black", "Other")),
    Exercise     = factor(RECEX, levels = c(1, 2, 3),
                          labels = c("High", "Medium", "Low")),
    Diabetes     = factor(to_num(DIAB), levels = c(0, 1),
                          labels = c("No", "Yes")),
    Hypertension = factor(to_num(HTN_REP), levels = c(0, 1),
                          labels = c("No", "Yes")),
    Cholesterol  = to_num(SERCHOL),
    HDL          = to_num(HDL),
    Triglycerides= to_num(TRIGLYC),
    Death        = factor(DEATH_num, levels = c(0, 1),
                          labels = c("Alive", "Dead"))
  ) %>%
  rename(
    Age       = AGEYRS,
    Education = GRADES,
    Smoking   = AVGSMK,
    Alcohol   = BOOZE
  )

cat("Final analytic sample:", nrow(df_clean), "\n")
cat("Deaths:", sum(df_clean$DEATH_num), "\n")
cat("Censored:", sum(df_clean$DEATH_num == 0), "\n")

library(dplyr)
library(skimr)
library(tidyr)
library(survival)
library(survminer)
library(haven)
library(broom)
library(rms)
library(stringr)
library(tableone)

# ----2. Table 1 --------------------------------------------

# ----3. Create survival time variables --------------------

# follow-up time in years (for KM/log-rank only)
df_clean <- df_clean %>%
  mutate(fu_time = age_exit - Age)

# main Surv object for Cox models (age as time scale)
surv_age <- with(df_clean,
                 Surv(time = Age, time2 = age_exit, event = DEATH_num))

# Surv object for KM / log-rank (time since baseline)
surv_fu  <- with(df_clean,
                 Surv(time = fu_time, event = DEATH_num))

# ----4. Unadjusted Kaplan-Meier curves and log-rank test --------------------

# KM curves by niacin quartile (time since baseline)
fit_km <- survfit(surv_fu ~ NIACIN_q, data = df_clean)
ggsurvplot(fit_km, 
           data = df_clean, 
           risk.table = TRUE,
           risk.table.height = 0.25,
           risk.table.fontsize = 3,
           legend.title = "Niacin Quartile",
           legend.labs  = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
           xlim = c(0, 18),
           break.time.by = 5,
           ggtheme = theme_minimal(base_family = "sans", base_size = 11),
           tables.theme = theme_cleantable(base_family = "sans", base_size = 9))

# Log-rank test across niacin quartiles (unadjusted)
logrank <- survdiff(surv_fu ~ NIACIN_q, data = df_clean)
logrank

# ----5. Primary Cox models (age-adjusted and multivariable) --------------------

# 4a. Age-adjusted Cox model (age as time scale, no covariates)
surv_obj <- with(df_clean,
                 Surv(time = Age, time2 = age_exit, event = DEATH_num)
)

cox_age <- coxph(
  surv_obj ~ NIACIN_q,
  data = df_clean
)
summary(cox_age)

# 4b. Multivariable Cox model (primary model)
cox_multiv <- coxph(
  surv_obj ~ NIACIN_q +
    Sex + Race + Education +
    Smoking + Alcohol + Exercise +
    BMI + Diabetes + Hypertension,
  data = df_clean,
  ties = "efron"
)
summary(cox_multiv)

# ----6. Check nonlinearity of niacin (primary Cox model) --------------------

# Compare linear vs quartile niacin in multivariable Cox model
cox_linear <- coxph(
  surv_obj ~ NIACIN +
    Sex + Race + Education +
    Smoking + Alcohol + Exercise +
    BMI + Diabetes + Hypertension,
  data = df_clean
)

# Likelihood ratio test for nonlinearity
anova(cox_linear, cox_multiv, test = "LRT")

# ----7. Check proportional hazards assumption --------------------

# Schoenfeld residual tests for NIACIN_q and global PH
ph_test <- cox.zph(cox_multiv)
ph_test

ggcoxzph(ph_test[1])

# ----8. Create Table 2 --------------------

library(broom)
library(dplyr)
library(stringr)

# ---- Get formatted HR(95% CI) for niacin quartiles from Cox models ----

format_niacin <- function(fit) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(str_detect(term, "^NIACIN_q")) %>%
    mutate(
      quart = str_replace(term, "NIACIN_q", ""),
      hr    = round(estimate, 2),
      lcl   = round(conf.low, 2),
      ucl   = round(conf.high, 2),
      HR_CI = sprintf("%.2f (%.2f, %.2f)", hr, lcl, ucl)
    ) %>%
    arrange(quart) %>%
    pull(HR_CI)
}

hr_age    <- format_niacin(cox_age)      # Q2–Q4
hr_multiv <- format_niacin(cox_multiv)   # Q2–Q4

# ---- Build Table 2 skeleton with real Cox numbers, Poisson placeholders ----

table2_data <- tibble(
  Characteristic = c(
    # Cox
    "Cox regression",
    "   Q1 (reference)",
    "   Q2",
    "   Q3",
    "   Q4",
    # Poisson
    "Poisson regression",
    "   Q1 (reference)",
    "   Q2",
    "   Q3",
    "   Q4"
  ),
  Age_adjusted = c(
    # Cox
    "",
    "1.00 (ref)",
    hr_age[1],
    hr_age[2],
    hr_age[3],
    # Poisson
    "",
    "1.00 (ref)",
    "x.xx (x.xx, x.xx)",
    "x.xx (x.xx, x.xx)",
    "x.xx (x.xx, x.xx)"
  ),
  Multivariable_adjusted = c(
    # Cox
    "",
    "1.00 (ref)",
    hr_multiv[1],
    hr_multiv[2],
    hr_multiv[3],
    # Poisson
    "",
    "1.00 (ref)",
    "x.xx (x.xx, x.xx)",
    "x.xx (x.xx, x.xx)",
    "x.xx (x.xx, x.xx)"
  )
)

section_rows <- c(1, 7)

ft2 <- flextable(table2_data) %>%
  set_header_labels(
    Characteristic         = "Characteristic",
    Age_adjusted           = "Age-adjusted\nHR (95% CI)",
    Multivariable_adjusted = "Multivariable-adjusted\nHR (95% CI)"
  ) %>%
  theme_booktabs() %>%
  bold(i = section_rows, j = 1, part = "body") %>%
  bold(part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(j = 1, align = "left",   part = "all") %>%
  align(j = 2:3, align = "center", part = "all") %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2:3, width = 1.5) %>%
  add_footer_lines(values = c(
    "HR denotes hazard ratio; IRR, incidence rate ratio.",
    "Age-adjusted model: adjusted for age (as timescale for Cox; as covariate for Poisson).",
    "Multivariable-adjusted model: additionally adjusted for sex, race, education, smoking, alcohol, exercise, BMI, diabetes, and hypertension.",
    "Q1: NIACIN ≤ 11.17 mg/day; Q2: 11.17–16.40 mg/day; Q3: 16.40–23.42 mg/day; Q4: > 23.42 mg/day."
  )) %>%
  add_header_lines("Table 2. Association of dietary niacin intake with all-cause mortality in Cox and Poisson regression models.") %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

save_as_docx(ft2, path = "/Users/cutiepie/epi204/table2.docx")
