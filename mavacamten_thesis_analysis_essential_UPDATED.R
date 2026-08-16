# =============================================================================
# Genetic Determinants of Treatment Response to Mavacamten in Obstructive HCM
# Essential verified statistical analysis code
#
# MSc Genomic Medicine dissertation
#
# This script documents the FINAL statistical models reported in the thesis.
# Patient-level SHaRe data are not included because of data-access and
# confidentiality restrictions.
#
# The script assumes that the final analytic datasets have already been created
# by the protected data-preparation workflow:
#
#   rest_03m_new, vals_03m_new, lvef_03m_new, nyha_03m_new
#   rest_9to15_new, vals_9to15_new, lvef_9to15_new, nyha_9to15_new
#   rest_gene_03m_new, vals_gene_03m_new,
#   lvef_gene_03m_new, nyha_gene_03m_new
#   lvot_long_new, lvef_long_new
#   persistence3
#
# Analysis definitions:
#   Baseline = most recent eligible assessment at or before mavacamten start.
#   0–3 month follow-up = latest eligible assessment >0 to <=3 months.
#   9–15 month follow-up = latest eligible assessment >9 to <=15 months.
#   Longitudinal sensitivity analyses = all post-treatment observations
#   from 0 to 15 months after mavacamten initiation.
#
# Genotype-negative and female are the reference categories.
# Age covariate = age at mavacamten initiation (mava_start_age).
# LVEF changes are expressed in percentage points.
# =============================================================================

library(dplyr)
library(emmeans)
library(lme4)
library(lmerTest)
library(broom)

# -----------------------------------------------------------------------------
# 0. Helper functions
# -----------------------------------------------------------------------------

set_primary_refs <- function(dat) {
  dat %>%
    mutate(
      GeneticGroup3 = relevel(
        factor(GeneticGroup3),
        ref = "Genotype-negative"
      ),
      Sex = relevel(
        factor(Sex),
        ref = "Female"
      )
    )
}

set_gene_refs <- function(dat) {
  dat %>%
    mutate(
      GeneGroup = relevel(
        factor(GeneGroup),
        ref = "Genotype-negative"
      ),
      Sex = relevel(
        factor(Sex),
        ref = "Female"
      )
    )
}

overall_p_lm <- function(model, term) {
  tab <- drop1(model, test = "F")
  tab[term, "Pr(>F)"]
}

overall_p_glm <- function(model, term) {
  tab <- drop1(model, test = "Chisq")
  tab[term, "Pr(>Chi)"]
}

show_linear_model <- function(model, group_term, group_variable) {
  cat("\nCoefficient table\n")
  print(broom::tidy(model, conf.int = TRUE))
  cat("\nOverall group effect P value\n")
  print(overall_p_lm(model, group_term))
  cat("\nAdjusted group means\n")
  print(
    summary(
      emmeans::emmeans(
        model,
        specs = as.formula(paste("~", group_variable))
      ),
      infer = c(TRUE, TRUE)
    )
  )
}

show_logistic_model <- function(model, group_term, group_variable) {
  cat("\nAdjusted odds ratios\n")
  print(
    broom::tidy(
      model,
      conf.int = TRUE,
      exponentiate = TRUE
    )
  )
  cat("\nOverall group effect P value\n")
  print(overall_p_glm(model, group_term))
  cat("\nAdjusted probabilities\n")
  print(
    summary(
      emmeans::emmeans(
        model,
        specs = as.formula(paste("~", group_variable)),
        type = "response"
      ),
      infer = c(TRUE, TRUE)
    )
  )
}


# =============================================================================
# 1. PRIMARY ANALYSES: 0–3 MONTHS
# =============================================================================

rest_03m_new <- set_primary_refs(rest_03m_new)
vals_03m_new <- set_primary_refs(vals_03m_new)
lvef_03m_new <- set_primary_refs(lvef_03m_new)
nyha_03m_new <- set_primary_refs(nyha_03m_new)

# Final analytic sample sizes
stopifnot(
  nrow(rest_03m_new) == 234,
  nrow(vals_03m_new) == 231,
  nrow(lvef_03m_new) == 235,
  nrow(nyha_03m_new) == 142
)


# -----------------------------------------------------------------------------
# 1.1 Resting LVOT gradient change
# Expected overall genotype effect: P = 0.658
# -----------------------------------------------------------------------------

rest_final <- lm(
  LVOT_change_rest ~
    GeneticGroup3 +
    LVOT_rest_baseline +
    mava_start_age +
    Sex,
  data = rest_03m_new
)

show_linear_model(
  rest_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 1.2 Valsalva LVOT gradient change
# Expected overall genotype effect: P = 0.515
# -----------------------------------------------------------------------------

vals_final <- lm(
  LVOT_change_valsalva ~
    GeneticGroup3 +
    LVOT_valsalva_baseline +
    mava_start_age +
    Sex,
  data = vals_03m_new
)

show_linear_model(
  vals_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 1.3 Continuous LVEF change
# Expected overall genotype effect: P = 0.097
# Main inference is the OVERALL genotype effect.
# -----------------------------------------------------------------------------

lvef_final <- lm(
  LVEF_change ~
    GeneticGroup3 +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_03m_new
)

show_linear_model(
  lvef_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 1.4 LVEF decline >=10 percentage points
# Expected overall genotype effect: P = 0.826
# -----------------------------------------------------------------------------

lvef_03m_new <- lvef_03m_new %>%
  mutate(
    LVEF_decline10 = if_else(
      LVEF_change <= -10,
      1L,
      0L
    )
  )

lvef10_final <- glm(
  LVEF_decline10 ~
    GeneticGroup3 +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_03m_new,
  family = binomial
)

show_logistic_model(
  lvef10_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)

# Observed event counts
lvef_03m_new %>%
  group_by(GeneticGroup3) %>%
  summarise(
    n = n(),
    decline_n = sum(LVEF_decline10),
    decline_percent = 100 * mean(LVEF_decline10),
    .groups = "drop"
  ) %>%
  print()


# -----------------------------------------------------------------------------
# 1.5 NYHA improvement >=1 functional class
# Expected overall genotype effect: P = 0.578
# -----------------------------------------------------------------------------

nyha_03m_new <- nyha_03m_new %>%
  mutate(
    NYHA_improved = if_else(
      NYHA_03m <= NYHA_Baseline - 1,
      1L,
      0L
    )
  )

nyha_final <- glm(
  NYHA_improved ~
    GeneticGroup3 +
    NYHA_Baseline +
    mava_start_age +
    Sex,
  data = nyha_03m_new,
  family = binomial
)

show_logistic_model(
  nyha_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)

# Observed event counts
nyha_03m_new %>%
  group_by(GeneticGroup3) %>%
  summarise(
    n = n(),
    improved_n = sum(NYHA_improved),
    improved_percent = 100 * mean(NYHA_improved),
    .groups = "drop"
  ) %>%
  print()


# =============================================================================
# 2. SECONDARY ANALYSES: 9–15 MONTHS
# =============================================================================

rest_9to15_new <- set_primary_refs(rest_9to15_new)
vals_9to15_new <- set_primary_refs(vals_9to15_new)
lvef_9to15_new <- set_primary_refs(lvef_9to15_new)
nyha_9to15_new <- set_primary_refs(nyha_9to15_new)

stopifnot(
  nrow(rest_9to15_new) == 228,
  nrow(vals_9to15_new) == 220,
  nrow(lvef_9to15_new) == 228,
  nrow(nyha_9to15_new) == 154
)


# -----------------------------------------------------------------------------
# 2.1 Resting LVOT gradient change
# Expected overall genotype effect: P = 0.728
# -----------------------------------------------------------------------------

rest_9to15_final <- lm(
  rest_change ~
    GeneticGroup3 +
    rest_baseline +
    mava_start_age +
    Sex,
  data = rest_9to15_new
)

show_linear_model(
  rest_9to15_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 2.2 Valsalva LVOT gradient change
# Expected overall genotype effect: P = 0.751
# -----------------------------------------------------------------------------

vals_9to15_final <- lm(
  vals_change ~
    GeneticGroup3 +
    vals_baseline +
    mava_start_age +
    Sex,
  data = vals_9to15_new
)

show_linear_model(
  vals_9to15_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 2.3 Continuous LVEF change
# Expected overall genotype effect: P = 0.538
# -----------------------------------------------------------------------------

lvef_9to15_final <- lm(
  LVEF_change ~
    GeneticGroup3 +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_9to15_new
)

show_linear_model(
  lvef_9to15_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)


# -----------------------------------------------------------------------------
# 2.4 LVEF decline >=10 percentage points
# Expected overall genotype effect: P = 0.938
# -----------------------------------------------------------------------------

lvef_9to15_new <- lvef_9to15_new %>%
  mutate(
    LVEF_decline10 = if_else(
      LVEF_change <= -10,
      1L,
      0L
    )
  )

lvef10_9to15_final <- glm(
  LVEF_decline10 ~
    GeneticGroup3 +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_9to15_new,
  family = binomial
)

show_logistic_model(
  lvef10_9to15_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)

# Observed event counts
lvef_9to15_new %>%
  group_by(GeneticGroup3) %>%
  summarise(
    n = n(),
    decline_n = sum(LVEF_decline10),
    decline_percent = 100 * mean(LVEF_decline10),
    .groups = "drop"
  ) %>%
  print()


# -----------------------------------------------------------------------------
# 2.5 NYHA improvement >=1 functional class
# Expected overall genotype effect: P = 0.696
# -----------------------------------------------------------------------------

# In the final prepared dataset NYHA_improved is defined as improvement
# of at least one functional class from baseline to the latest >9 to <=15
# month assessment. Recreate it here only if the follow-up variable exists.
if ("NYHA_9to15" %in% names(nyha_9to15_new)) {
  nyha_9to15_new <- nyha_9to15_new %>%
    mutate(
      NYHA_improved = if_else(
        NYHA_9to15 <= NYHA_Baseline - 1,
        1L,
        0L
      )
    )
}

stopifnot("NYHA_improved" %in% names(nyha_9to15_new))

nyha_9to15_final <- glm(
  NYHA_improved ~
    GeneticGroup3 +
    NYHA_Baseline +
    mava_start_age +
    Sex,
  data = nyha_9to15_new,
  family = binomial
)

show_logistic_model(
  nyha_9to15_final,
  group_term = "GeneticGroup3",
  group_variable = "GeneticGroup3"
)

nyha_9to15_new %>%
  group_by(GeneticGroup3) %>%
  summarise(
    n = n(),
    improved_n = sum(NYHA_improved),
    improved_percent = 100 * mean(NYHA_improved),
    .groups = "drop"
  ) %>%
  print()


# =============================================================================
# 3. EXPLORATORY GENE-LEVEL ANALYSES: 0–3 MONTHS
#    Genotype-negative vs MYBPC3 vs MYH7
# =============================================================================

rest_gene_03m_new <- set_gene_refs(rest_gene_03m_new)
vals_gene_03m_new <- set_gene_refs(vals_gene_03m_new)
lvef_gene_03m_new <- set_gene_refs(lvef_gene_03m_new)
nyha_gene_03m_new <- set_gene_refs(nyha_gene_03m_new)

stopifnot(
  nrow(rest_gene_03m_new) == 178,
  nrow(vals_gene_03m_new) == 177,
  nrow(lvef_gene_03m_new) == 179,
  nrow(nyha_gene_03m_new) == 107
)


# -----------------------------------------------------------------------------
# 3.1 Resting LVOT
# Expected overall gene-group effect: P = 0.440
# -----------------------------------------------------------------------------

rest_gene_final <- lm(
  LVOT_change_rest ~
    GeneGroup +
    LVOT_rest_baseline +
    mava_start_age +
    Sex,
  data = rest_gene_03m_new
)

show_linear_model(
  rest_gene_final,
  group_term = "GeneGroup",
  group_variable = "GeneGroup"
)


# -----------------------------------------------------------------------------
# 3.2 Valsalva LVOT
# Expected overall gene-group effect: P = 0.651
# -----------------------------------------------------------------------------

vals_gene_final <- lm(
  LVOT_change_valsalva ~
    GeneGroup +
    LVOT_valsalva_baseline +
    mava_start_age +
    Sex,
  data = vals_gene_03m_new
)

show_linear_model(
  vals_gene_final,
  group_term = "GeneGroup",
  group_variable = "GeneGroup"
)


# -----------------------------------------------------------------------------
# 3.3 Continuous LVEF
# Expected overall gene-group effect: P = 0.070
# Tukey-adjusted pairwise comparisons are exploratory and non-significant.
# -----------------------------------------------------------------------------

lvef_gene_final <- lm(
  LVEF_change ~
    GeneGroup +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_gene_03m_new
)

show_linear_model(
  lvef_gene_final,
  group_term = "GeneGroup",
  group_variable = "GeneGroup"
)

lvef_gene_emm <- emmeans(
  lvef_gene_final,
  ~ GeneGroup
)

summary(
  pairs(
    lvef_gene_emm,
    adjust = "tukey"
  ),
  infer = c(TRUE, TRUE)
)


# -----------------------------------------------------------------------------
# 3.4 LVEF decline >=10 percentage points
# Expected overall gene-group effect: P = 0.441
# -----------------------------------------------------------------------------

lvef_gene_03m_new <- lvef_gene_03m_new %>%
  mutate(
    LVEF_decline10 = if_else(
      LVEF_change <= -10,
      1L,
      0L
    )
  )

lvef10_gene_final <- glm(
  LVEF_decline10 ~
    GeneGroup +
    LVEF_baseline +
    mava_start_age +
    Sex,
  data = lvef_gene_03m_new,
  family = binomial
)

show_logistic_model(
  lvef10_gene_final,
  group_term = "GeneGroup",
  group_variable = "GeneGroup"
)


# -----------------------------------------------------------------------------
# 3.5 NYHA improvement >=1 class
# Expected overall gene-group effect: P = 0.256
# -----------------------------------------------------------------------------

stopifnot("NYHA_improved" %in% names(nyha_gene_03m_new))

nyha_gene_final <- glm(
  NYHA_improved ~
    GeneGroup +
    NYHA_Baseline +
    mava_start_age +
    Sex,
  data = nyha_gene_03m_new,
  family = binomial
)

show_logistic_model(
  nyha_gene_final,
  group_term = "GeneGroup",
  group_variable = "GeneGroup"
)


# =============================================================================
# 4. LONGITUDINAL MIXED-EFFECTS SENSITIVITY ANALYSES
#    Post-treatment observations only: 0–15 months
# =============================================================================

lvot_long_new <- set_primary_refs(lvot_long_new)
lvef_long_new <- set_primary_refs(lvef_long_new)

stopifnot(
  nrow(lvot_long_new) == 1521,
  dplyr::n_distinct(lvot_long_new$HCM_DCC_PID) == 364,
  nrow(lvef_long_new) == 1523,
  dplyr::n_distinct(lvef_long_new$HCM_DCC_PID) == 366,
  min(lvot_long_new$months_from_mavacamten, na.rm = TRUE) >= 0,
  max(lvot_long_new$months_from_mavacamten, na.rm = TRUE) <= 15,
  min(lvef_long_new$months_from_mavacamten, na.rm = TRUE) >= 0,
  max(lvef_long_new$months_from_mavacamten, na.rm = TRUE) <= 15
)


# -----------------------------------------------------------------------------
# 4.1 Resting LVOT longitudinal model
# Primary longitudinal inference:
# overall genotype-by-time interaction P = 0.076
# -----------------------------------------------------------------------------

model_lvot_mixed_new <- lmerTest::lmer(
  resting_LVOT ~
    months_from_mavacamten * GeneticGroup3 +
    mava_start_age +
    Sex +
    (1 | HCM_DCC_PID),
  data = lvot_long_new,
  REML = FALSE
)

summary(model_lvot_mixed_new)

# Overall fixed-effect tests, including genotype-by-time interaction
anova(model_lvot_mixed_new)

# Group-specific monthly slopes for interpretation
lvot_slopes <- emmeans::emtrends(
  model_lvot_mixed_new,
  ~ GeneticGroup3,
  var = "months_from_mavacamten"
)

summary(
  lvot_slopes,
  infer = c(TRUE, TRUE)
)


# -----------------------------------------------------------------------------
# 4.2 LVEF longitudinal model
# Primary longitudinal inference:
# overall genotype-by-time interaction P = 0.835
# -----------------------------------------------------------------------------

model_lvef_mixed_new <- lmerTest::lmer(
  lvef ~
    months_from_mavacamten * GeneticGroup3 +
    mava_start_age +
    Sex +
    (1 | HCM_DCC_PID),
  data = lvef_long_new,
  REML = FALSE
)

summary(model_lvef_mixed_new)

anova(model_lvef_mixed_new)

lvef_slopes <- emmeans::emtrends(
  model_lvef_mixed_new,
  ~ GeneticGroup3,
  var = "months_from_mavacamten"
)

summary(
  lvef_slopes,
  infer = c(TRUE, TRUE)
)


# =============================================================================
# 5. TREATMENT PERSISTENCE AND DISCONTINUATION
#    Observation periods vary between patients; this is NOT a fixed
#    15-month landmark analysis.
# =============================================================================

persistence3 <- persistence3 %>%
  filter(
    GeneticGroup %in% c(
      "Genotype-negative",
      "Sarcomeric P/LP",
      "VUS"
    )
  ) %>%
  mutate(
    GeneticGroup = factor(
      GeneticGroup,
      levels = c(
        "Genotype-negative",
        "Sarcomeric P/LP",
        "VUS"
      )
    )
  )

stopifnot(
  nrow(persistence3) == 431,
  dplyr::n_distinct(persistence3$HCM_MedicationList_PID) == 431
)

# Treatment-status summary
persistence_summary <- persistence3 %>%
  count(
    GeneticGroup,
    treatment_status
  ) %>%
  group_by(GeneticGroup) %>%
  mutate(
    total = sum(n),
    percent = 100 * n / total
  ) %>%
  ungroup()

print(persistence_summary)

# Association between genotype and treatment status
persistence_tab <- table(
  persistence3$GeneticGroup,
  persistence3$treatment_status
)

print(persistence_tab)
chisq.test(persistence_tab)
fisher.test(persistence_tab)

# Time to discontinuation among patients who discontinued and had both
# medication start and end ages available.
persistence3 <- persistence3 %>%
  mutate(
    time_to_discontinuation_months =
      (Med_End_Age - Med_Start_Age) * 12
  )

time_to_discontinuation <- persistence3 %>%
  filter(
    treatment_status == "Discontinued",
    !is.na(time_to_discontinuation_months)
  ) %>%
  group_by(GeneticGroup) %>%
  summarise(
    n = n(),
    median_months = median(
      time_to_discontinuation_months,
      na.rm = TRUE
    ),
    Q1 = quantile(
      time_to_discontinuation_months,
      0.25,
      na.rm = TRUE
    ),
    Q3 = quantile(
      time_to_discontinuation_months,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(time_to_discontinuation)

# Check timing availability amongst discontinuers
persistence3 %>%
  filter(treatment_status == "Discontinued") %>%
  group_by(GeneticGroup) %>%
  summarise(
    discontinued = n(),
    time_available = sum(
      !is.na(Med_Start_Age) &
        !is.na(Med_End_Age)
    ),
    missing_time = sum(
      is.na(Med_Start_Age) |
        is.na(Med_End_Age)
    ),
    .groups = "drop"
  ) %>%
  print()


# =============================================================================
# 6. KEY REPRODUCIBILITY CHECKS
# =============================================================================

cat("\nExpected final analytic sample sizes\n")
cat("0–3 months: resting 234; Valsalva 231; LVEF 235; NYHA 142\n")
cat("9–15 months: resting 228; Valsalva 220; LVEF 228; NYHA 154\n")
cat("Gene-level 0–3 months: resting 178; Valsalva 177; LVEF 179; NYHA 107\n")
cat("Longitudinal resting LVOT: 1,521 observations / 364 patients\n")
cat("Longitudinal LVEF: 1,523 observations / 366 patients\n")
cat("Persistence: 431 patients\n\n")

cat("Expected overall P values\n")
cat("0–3 months: resting .658; Valsalva .515; LVEF .097; LVEF10 .826; NYHA .578\n")
cat("9–15 months: resting .728; Valsalva .751; LVEF .538; LVEF10 .938; NYHA .696\n")
cat("Gene level: resting .440; Valsalva .651; LVEF .070; LVEF10 .441; NYHA .256\n")
cat("Mixed-effects interactions: resting LVOT .076; LVEF .835\n")
cat("Persistence: Pearson chi-square .358; Fisher exact .362\n")

# END OF ESSENTIAL VERIFIED ANALYSIS SCRIPT
