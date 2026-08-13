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

# 13. Treatment persistence/discontinuation ------------------------------

table(
  is.na(persistence3$Med_End_Age),
  persistence3$treatment_status,
  useNA = "ifany"
)

persistence_summary <- persistence3 |>
  count(GeneticGroup, treatment_status) |>
  group_by(GeneticGroup) |>
  mutate(
    total = sum(n),
    percent = 100 * n / total
  )

persistence_summary

persistence_tab <- table(
  persistence3$GeneticGroup,
  persistence3$treatment_status
)
chisq.test(persistence_tab)
fisher.test(persistence_tab)

# 14. Reproducibility checks ----------------------------------------------
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
