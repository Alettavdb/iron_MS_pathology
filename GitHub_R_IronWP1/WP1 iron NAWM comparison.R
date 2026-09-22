#Aletta
#Script iron WP1 NAWM comparison

library(plyr)
library(lme4)
library(car)
library(readxl)
library(ggplot2)
library(emmeans)
library(ggridges)
library(ggsignif)
library(psych)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(emmeans)
library(openxlsx)
library(reshape2)
library(ggrepel)
library(broom)
library(scales)
library(dplyr)
library(patchwork)
library(geomtextpath)
library(chisq.posthoc.test)
library(glmmTMB)
library(ggnewscale)
library(rstatix)


setwd("T:/Lab Absinta/15_ALETTA/Iron/WP1")


#Are NAWM patterns are associated with donor characteristics? Fig. 1D

vars_numeric <- c(
  "age", "weight", "ph", "pmd_min",
  "DOD", "ARMSS", "Severity_score"
)


Donor_TBB <- read_excel("Data analysis/Datasets/Donor_data_TBB.xlsx")

NAWM_plot <- Donor_TBB %>%
  filter(NAWM_pattern %in% c("G-R-", "G-R+", "G+R-", "G+R+")) %>%
  select(
    NBB,
    NAWM_pattern,
    sex,
    all_of(vars_numeric)
  )

NAWM_donor <- NAWM_plot %>%
  group_by(NBB, NAWM_pattern) %>%
  summarise(
    sex = first(na.omit(sex)),
    across(all_of(vars_numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

library(rstatix)

stats_results <- NAWM_donor %>%
  pivot_longer(
    cols = all_of(vars_numeric),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  kruskal_test(value ~ NAWM_pattern) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(significant = p.adj < 0.05)

heatmap_data <- NAWM_donor %>%
  group_by(NAWM_pattern) %>%
  summarise(
    across(all_of(vars_numeric), ~ mean(.x, na.rm = TRUE))
  ) %>%
  column_to_rownames("NAWM_pattern") %>%
  scale()

sig_labels <- stats_results %>%
  mutate(label = ifelse(significant, "*", "")) %>%
  select(variable, label) %>%
  deframe()

sex_summary <- NAWM_donor %>%
  distinct(NBB, NAWM_pattern, sex) %>%
  group_by(NAWM_pattern) %>%
  summarise(
    prop_F = mean(sex == "F", na.rm = TRUE),
    .groups = "drop"
  )

sex_test <- NAWM_donor %>%
  distinct(NBB, NAWM_pattern, sex) %>%
  count(NAWM_pattern, sex) %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0)

fisher.test(as.matrix(sex_test[, -1]))

library(ComplexHeatmap)
library(circlize)

sex_annot <- HeatmapAnnotation(
  `Prop. female` = sex_summary$prop_F,
  col = list(
    `Prop. female` =
      colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  ),
  annotation_name_side = "left"
)

var_order <- vars_numeric
pattern_order <- c("G-R-", "G-R+", "G+R-", "G+R+")

Heatmap(
  t(heatmap_data)[var_order, pattern_order, drop = FALSE],
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_labels = paste0(
    var_order,
    sig_labels[var_order]
  ),
  row_title = "Donor characteristics",
  top_annotation = sex_annot,
  heatmap_legend_param = list(title = "Scaled mean")
)

#we now add count data

# Count of donors per NAWM pattern
pattern_counts <- NAWM_donor %>%
  count(NAWM_pattern) %>%
  arrange(factor(NAWM_pattern, levels = pattern_order))

# Top annotation: both prop female and counts
top_annot <- HeatmapAnnotation(
  `Donor count` = anno_barplot(pattern_counts$n,
                               gp = gpar(fill = "grey"),
                               border = FALSE),
  `Prop. female` = sex_summary$prop_F,
  col = list(
    `Prop. female` = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  ),
  annotation_name_side = "left"
)

Heatmap(
  t(heatmap_data)[var_order, pattern_order, drop = FALSE],
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_labels = paste0(var_order, sig_labels[var_order]),
  row_title = "Donor characteristics",
  top_annotation = top_annot,
  heatmap_legend_param = list(title = "Scaled mean")
)


#Okay now lets take out the donor characteristics and make a separate one where we cluster based on pathology


vars_numeric <- c(
  "Reactive_load", "Lesion_load", "Prop_act", "Prop_mixed", "Prop_inact", "Prop_remy", "MMAS"
)

NAWM_plot <- Donor_TBB %>%
  filter(NAWM_pattern %in% c("G-R-", "G-R+", "G+R-", "G+R+")) %>%
  select(
    NBB,
    NAWM_pattern,
    all_of(vars_numeric)
  )

NAWM_donor <- NAWM_plot %>%
  group_by(NBB, NAWM_pattern) %>%
  summarise(
    across(all_of(vars_numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

stats_results <- NAWM_donor %>%
  pivot_longer(
    cols = all_of(vars_numeric),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  kruskal_test(value ~ NAWM_pattern) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(significant = p.adj < 0.05)

heatmap_data <- NAWM_donor %>%
  group_by(NAWM_pattern) %>%
  summarise(
    across(all_of(vars_numeric), ~ mean(.x, na.rm = TRUE))
  ) %>%
  column_to_rownames("NAWM_pattern") %>%
  scale()

sig_labels <- stats_results %>%
  mutate(label = ifelse(significant, "*", "")) %>%
  select(variable, label) %>%
  deframe()

library(ComplexHeatmap)
library(circlize)

Heatmap(
  t(heatmap_data),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # CLUSTER BOTH
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  row_names_side = "left",
  
  row_labels = paste0(
    rownames(t(heatmap_data)),
    sig_labels[rownames(t(heatmap_data))]
  ),
  
  row_title = "Severity characteristics",
  heatmap_legend_param = list(title = "Scaled mean")
)

