# MAVACAMTEN HCM THESIS
# Essential Statistical Analysis Code
#
# This script documents the verified statistical analyses.
# Patient-level source data are not included.

library(dplyr)
library(emmeans)
library(car)
library(lmerTest)
library(lme4)
library(ggplot2)
library(survival)

# 1. Primary cohort -------------------------------------------------------
stopifnot(dplyr::n_distinct(mava_cohort$HCM_DCC_PID) == 605)

genotype_cohort <- mava_cohort |>
  filter(GeneticGroup3 %in% c(
    "Genotype-negative", "Sarcomeric P/LP", "VUS"
  ))

stopifnot(
  nrow(genotype_cohort) == 440,
  dplyr::n_distinct(genotype_cohort$HCM_DCC_PID) == 440
)

# 2. Primary 0-3 month resting LVOT --------------------------------------
model_rest_03m <- lm(
  LVOT_change_rest ~ GeneticGroup3 + LVOT_rest_baseline +
    FirstEncounter_Age + Sex,
  data = rest_lvot_primary
)
summary(model_rest_03m)
car::Anova(model_rest_03m, type = 2)
summary(emmeans(model_rest_03m, ~ GeneticGroup3), infer = c(TRUE, TRUE))

# 3. Primary 0-3 month Valsalva LVOT -------------------------------------
model_valsalva_3group <- lm(
  LVOT_change_valsalva ~ GeneticGroup3 +
    echo_LVOT_valsalva_gradient0 + FirstEncounter_Age + Sex,
  data = genetic_sensitivity
)
summary(model_valsalva_3group)
car::Anova(model_valsalva_3group, type = 2)
summary(emmeans(model_valsalva_3group, ~ GeneticGroup3),
        infer = c(TRUE, TRUE))

# 4. Primary 0-3 month LVEF change ---------------------------------------
model_lvef_change_3group <- lm(
  LVEF_change ~ GeneticGroup3 + echo_LVEF0 +
    FirstEncounter_Age + Sex,
  data = genetic_sensitivity
)
summary(model_lvef_change_3group)
car::Anova(model_lvef_change_3group, type = 2)
summary(emmeans(model_lvef_change_3group, ~ GeneticGroup3),
        infer = c(TRUE, TRUE))

# 5. LVEF decline >=10 percentage points -------------------------------
model_lvef_drop_3group <- glm(
  LVEF_drop10 ~ GeneticGroup3 + echo_LVEF0 +
    FirstEncounter_Age + Sex,
  family = binomial,
  data = genetic_sensitivity
)
summary(model_lvef_drop_3group)
car::Anova(model_lvef_drop_3group, type = 2)
exp(cbind(
  OR = coef(model_lvef_drop_3group),
  confint(model_lvef_drop_3group)
))

# 6. Primary 0-3 month NYHA change ---------------------------------------
model_nyha_change_0_3m <- lm(
  NYHA_Change ~ GeneticGroup + NYHA_Baseline +
    FirstEncounter_Age + Sex,
  data = nyha_analysis_0_3m
)
summary(model_nyha_change_0_3m)
car::Anova(model_nyha_change_0_3m, type = 2)
summary(emmeans(model_nyha_change_0_3m, ~ GeneticGroup),
        infer = c(TRUE, TRUE))

# 7. 9-15 month LVEF -----------------------------------------------------
model_lvef_9to15 <- lm(
  LVEF_change ~ GeneticGroup + baseline_LVEF +
    FirstEncounter_Age + Sex,
  data = lvef_model_data_9to15
)
summary(model_lvef_9to15)
car::Anova(model_lvef_9to15, type = 2)
summary(emmeans(model_lvef_9to15, ~ GeneticGroup),
        infer = c(TRUE, TRUE))

# 8. 9-15 month NYHA -----------------------------------------------------
model_nyha_9to15 <- lm(
  followup_NYHA ~ GeneticGroup + baseline_NYHA +
    FirstEncounter_Age + Sex,
  data = nyha_model_data_9to15
)
summary(model_nyha_9to15)
car::Anova(model_nyha_9to15, type = 2)
summary(emmeans(model_nyha_9to15, ~ GeneticGroup),
        infer = c(TRUE, TRUE))

# 9. 9-15 month NYHA change ----------------------------------------------
nyha_change_model_9to15 <- lm(
  NYHA_change ~ GeneticGroup + baseline_NYHA +
    FirstEncounter_Age + Sex,
  data = nyha_model_data_9to15
)
summary(nyha_change_model_9to15)
car::Anova(nyha_change_model_9to15, type = 2)
summary(emmeans(nyha_change_model_9to15, ~ GeneticGroup),
        infer = c(TRUE, TRUE))

# 10. Exploratory MYBPC3/MYH7 LVEF analysis -------------------------------
# Gene-level source cohort: MYBPC3 n=58, MYH7 n=43.
# Analytic subset in this model: genotype-negative n=163,
# MYBPC3 n=34, MYH7 n=16; total n=213.
model_lvef_gene_change <- lm(
  LVEF_change ~ GeneGroup + baseline_LVEF +
    FirstEncounter_Age + Sex,
  data = lvef_gene_03m
)
summary(model_lvef_gene_change)
car::Anova(model_lvef_gene_change, type = 2)

emmeans_lvef_gene <- emmeans(model_lvef_gene_change, ~ GeneGroup)
summary(emmeans_lvef_gene, infer = c(TRUE, TRUE))
pairs(emmeans_lvef_gene, adjust = "tukey")

# 11. Longitudinal mixed-effects LVEF -------------------------------
model_lvef_mixed <- lmerTest::lmer(
  lvef ~ months_from_start * GeneticGroup +
    FirstEncounter_Age + Sex + (1 | HCM_DCC_PID),
  data = lvef_mixed_0to15,
  REML = FALSE
)
summary(model_lvef_mixed)
anova(model_lvef_mixed)

# 12. Longitudinal mixed-effects resting LVOT -------------------------------
lvot_mixed_0to15 <- echo_rest_long |>
  filter(
    months_from_mavacamten >= 0,
    months_from_mavacamten <= 15
  ) |>
  left_join(
    mava_cohort |>
      select(
        HCM_DCC_PID,
        GeneticGroup3,
        FirstEncounter_Age,
        Sex
      ) |>
      distinct(),
    by = "HCM_DCC_PID"
  )

model_lvot_mixed <- lmerTest::lmer(
  resting_LVOT ~ months_from_mavacamten * GeneticGroup3 +
    FirstEncounter_Age + Sex + (1 | HCM_DCC_PID),
  data = lvot_mixed_0to15 |>
    filter(!is.na(GeneticGroup3)),
  REML = FALSE
)

summary(model_lvot_mixed)
anova(model_lvot_mixed)

lvot_slopes <- emmeans::emtrends(
  model_lvot_mixed,
  ~ GeneticGroup3,
  var = "months_from_mavacamten"
)

summary(lvot_slopes, infer = c(TRUE, TRUE))
pairs(lvot_slopes, adjust = "tukey")

lvot_pred <- emmeans::emmeans(
  model_lvot_mixed,
  ~ GeneticGroup3 | months_from_mavacamten,
  at = list(months_from_mavacamten = seq(0, 15, by = 1))
)

lvot_pred_df <- as.data.frame(lvot_pred)

ggplot(
  lvot_pred_df,
  aes(
    x = months_from_mavacamten,
    y = emmean,
    group = GeneticGroup3,
    colour = GeneticGroup3
  )
) +
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      fill = GeneticGroup3
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c(
      "Genotype-negative" = "#F8766D",
      "Sarcomeric P/LP" = "#00BA38",
      "VUS" = "#619CFF"
    ),
    name = "Genotype"
  ) +
  scale_fill_manual(
    values = c(
      "Genotype-negative" = "#F8766D",
      "Sarcomeric P/LP" = "#00BA38",
      "VUS" = "#619CFF"
    ),
    name = "Genotype"
  ) +
  scale_x_continuous(
    breaks = seq(0, 15, by = 3),
    limits = c(0, 15)
  ) +
  labs(
    x = "Months after Mavacamten initiation",
    y = "Adjusted resting LVOT gradient (mmHg)"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# 13. Treatment persistence: Kaplan–Meier analysis ------------------------
#
# km_data is created upstream from the medication records and last available
# clinic follow-up. It should contain:
#   HCM_DCC_PID, GeneticGroup3, Med_Start_Age, Med_End_Age,
#   treatment_status, event, followup_years
#
# event = 1 for recorded treatment discontinuation and 0 for patients who
# remained ongoing and were censored at their last available clinic visit.
# Time is measured in years from mavacamten initiation.

stopifnot(exists("km_data"))

km_analysis <- km_data |>
  filter(
    !is.na(GeneticGroup3),
    !is.na(Med_Start_Age),
    !is.na(followup_years),
    followup_years >= 0
  ) |>
  mutate(
    GeneticGroup3 = factor(
      GeneticGroup3,
      levels = c(
        "Genotype-negative",
        "Sarcomeric P/LP",
        "VUS"
      )
    )
  )

# Basic checks for the final KM analysis cohort
km_analysis |>
  summarise(
    N = n(),
    events = sum(event),
    censored = sum(event == 0),
    median_followup = median(followup_years),
    max_followup = max(followup_years)
  )

km_analysis |>
  count(GeneticGroup3, event)

# Kaplan–Meier model
km_fit <- survival::survfit(
  survival::Surv(followup_years, event) ~ GeneticGroup3,
  data = km_analysis
)

summary(km_fit)$table

# Log-rank comparison between genotype groups
km_logrank <- survival::survdiff(
  survival::Surv(followup_years, event) ~ GeneticGroup3,
  data = km_analysis
)

km_logrank

p_logrank <- 1 - pchisq(
  km_logrank$chisq,
  df = length(km_logrank$n) - 1
)

p_logrank

# Kaplan–Meier median estimates and confidence intervals
km_table <- as.data.frame(summary(km_fit)$table)
km_table$GeneticGroup <- sub(
  "GeneticGroup3=",
  "",
  rownames(km_table)
)
rownames(km_table) <- NULL

km_table <- km_table[, c(
  "GeneticGroup",
  "records",
  "events",
  "median",
  "X0.95LCL",
  "X0.95UCL"
)]

names(km_table) <- c(
  "GeneticGroup",
  "N",
  "Discontinuations",
  "Median_years",
  "Lower_95CI",
  "Upper_95CI"
)

km_table

# Kaplan–Meier figure
plot(
  km_fit,
  col = c("#0072B2", "#D55E00", "#009E73"),
  lwd = 2,
  xlim = c(0, 5),
  ylim = c(0, 1),
  xlab = "Time since mavacamten initiation (years)",
  ylab = "Probability of remaining on mavacamten",
  main = "Kaplan–Meier estimates of mavacamten treatment persistence according to genotype",
  mark.time = FALSE
)

legend(
  "right",
  legend = levels(km_analysis$GeneticGroup3),
  col = c("#0072B2", "#D55E00", "#009E73"),
  lwd = 2,
  bty = "n",
  title = "Genotype"
)

# 14. Reproducibility checks ----------------------------------------------
mava_cohort |>
  summarise(
    N = n(),
    unique_patients = n_distinct(HCM_DCC_PID)
  )

genotype_cohort |>

# 14. Table A1 baseline P values ---------------------------------------------
# Overall unadjusted comparisons across the three genotype groups.
# Continuous variables: Kruskal-Wallis test.
# Binary variables: chi-squared test, or Fisher's exact test when any cell
# count is <5. Empty categories are handled safely.
#
# table1_cohort should contain one row per patient in the three principal
# genotype groups and the treatment-aligned Table 1 variables.

stopifnot(exists("table1_cohort"))

table1_cohort <- table1_cohort |>
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

table1_p_value <- function(data, variable, type = c("continuous", "binary")) {

  type <- match.arg(type)

  x <- data[[variable]]
  g <- data$GeneticGroup

  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]

  if (type == "continuous") {
    if (length(unique(x)) < 2 || length(unique(g)) < 2) {
      return(NA_real_)
    }

    return(kruskal.test(x, g)$p.value)
  }

  tab <- table(g, x)

  tab <- tab[
    rowSums(tab) > 0,
    colSums(tab) > 0,
    drop = FALSE
  ]

  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NA_real_)
  }

  if (any(tab < 5)) {
    return(fisher.test(tab)$p.value)
  }

  chisq.test(tab)$p.value
}

format_table1_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf("%.3f", p)
}

table1_variables <- tibble::tribble(
  ~Characteristic, ~variable, ~type,
  "Age at mavacamten initiation, years", "mava_start_age", "continuous",
  "Female sex", "female", "binary",
  "Hypertension", "hypertension", "binary",
  "Atrial fibrillation", "atrial_fibrillation", "binary",
  "Family history of HCM", "family_history_HCM", "binary",
  "Family history of sudden cardiac death", "family_history_SCD", "binary",
  "Body mass index, kg/m²", "baseline_BMI", "continuous",
  "NYHA functional class II", "NYHA_II", "binary",
  "NYHA functional class III", "NYHA_III", "binary",
  "Implantable cardioverter-defibrillator", "ICD", "binary",
  "Septal reduction therapy", "septal_reduction_therapy", "binary",
  "β-blocker", "beta_blocker", "binary",
  "Calcium channel blocker", "calcium_channel_blocker", "binary",
  "Disopyramide", "disopyramide", "binary",
  "Diuretic", "diuretic", "binary",
  "LVEF, %", "baseline_LVEF", "continuous",
  "Maximum LV wall thickness, mm", "baseline_MaxLVT", "continuous",
  "Resting LVOT gradient, mmHg", "baseline_rest_LVOT", "continuous",
  "Valsalva LVOT gradient, mmHg", "baseline_valsalva_LVOT", "continuous",
  "Post-exercise LVOT gradient, mmHg", "baseline_exercise_LVOT", "continuous"
)

table1_p_values <- table1_variables |>
  rowwise() |>
  mutate(
    p_value = table1_p_value(
      table1_cohort,
      variable,
      type
    ),
    `P value` = format_table1_p(p_value)
  ) |>
  ungroup() |>
  select(
    Characteristic,
    `P value`
  )

table1_p_values

# P values are unadjusted overall comparisons across the three genotype
# groups and are intended for Table 1/Table A1.

# 15. Reproducibility checks -------------------------------------------------
mava_cohort |>
  summarise(
    N = n(),
    unique_patients = n_distinct(HCM_DCC_PID)
  )

genotype_cohort |>
  summarise(
    N = n(),
    unique_patients = n_distinct(HCM_DCC_PID)
  )

persistence3 |>
  summarise(
    N = n(),
    unique_patients = n_distinct(HCM_MedicationList_PID),
    discontinued = sum(treatment_status == "Discontinued"),
    ongoing = sum(treatment_status == "Ongoing")
  )

# END OF ESSENTIAL ANALYSIS SCRIPT
