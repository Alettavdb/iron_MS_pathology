#Aletta
#Script iron WP1 classic lesion comparison

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

#TBB Fig.1G
Results_combined <- read_excel("Data analysis/Datasets/Results_combined.xlsx")

summary_df <- Results_combined %>%
  group_by(Lesion_Region_Microglia) %>%
  summarise(
    mean_TBB = mean(TBB_per_mm2, na.rm = TRUE),
    sem_TBB = sd(TBB_per_mm2, na.rm = TRUE) / sqrt(n())
  )

order_levels <- c(
  "NAWM_NAWM_",
  "Reactive_lesion_",
  "Active_lesion_Non-foamy",
  "Active_peri_Non-foamy",
  "Active_lesion_Foamy",
  "Active_peri_Foamy",
  "Mixed_core_Non-foamy",
  "Mixed_border_Non-foamy",
  "Mixed_peri_Non-foamy",
  "Mixed_core_Foamy",
  "Mixed_border_Foamy",
  "Mixed_peri_Foamy",
  "Inactive_lesion_",
  "Inactive_peri_",
  "Remyelinated_lesion_",
  "Remyelinated_peri_"
)

ggplot() +
  # Mean ± SEM bars
  geom_col(data = summary_df, 
           aes(x = factor(Lesion_Region_Microglia, levels = order_levels),
               y = mean_TBB),
           width = 0.6, fill = "white", color = "black") +
  geom_errorbar(data = summary_df,
                aes(x = factor(Lesion_Region_Microglia, levels = order_levels),
                    ymin = mean_TBB - sem_TBB, ymax = mean_TBB + sem_TBB),
                width = 0.2, color = "black") +
  # Individual data points
  geom_jitter(data = Results_combined,
              aes(x = factor(Lesion_Region_Microglia, levels = order_levels),
                  y = TBB_per_mm2),
              width = 0.15, size = 0.7, color = "black", alpha = 0.6) +
  labs(x = "Lesion type", y = "TBB per mm2") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = c(1.5, 2.5, 4.5, 6.5, 9.5, 12.5, 14.5), 
             linetype = "dashed", color = "grey40")

model <- lmer(
  TBB_per_mm2 ~ Lesion_Region_Microglia + (1 | Donor),
  data = Results_combined
)

emm <- emmeans(model, ~ Lesion_Region_Microglia)

pairs(emm, adjust = "fdr")




#PLIN2 Suppl.Fig.1H

summary_plin2 <- Results_combined %>%
  group_by(Region_grouped_microglia) %>%
  summarise(
    mean_PLIN2 = mean(dens_HLApos_PLIN2pos, na.rm = TRUE),
    sem_PLIN2 = sd(dens_HLApos_PLIN2pos, na.rm = TRUE) / sqrt(n())
  )

order_levels <- c(
  "NAWM",
  "core_Non-foamy",
  "core_Foamy",
  "HLA_accu_Non-foamy",
  "HLA_accu_Foamy",
  "peri_Non-foamy",
  "peri_Foamy"
)

ggplot(Results_combined,
       aes(x = factor(Region_grouped_microglia, levels = order_levels),
           y = dens_HLApos_PLIN2pos)) +
  
  # Boxplot
  geom_boxplot(
    width = 0.6,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  
  # Individual data points
  geom_jitter(
    width = 0.15,
    size = 0.7,
    color = "black",
    alpha = 0.6
  ) +
  
  labs(
    x = "Lesion type",
    y = "PLIN2+ HLA+ cells/mm2"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

df2 <- Results_combined %>%
  filter(!is.na(Region_grouped_microglia),
         Region_grouped_microglia != "NA")

model <- lmer(
  dens_HLApos_PLIN2pos ~ Region_grouped_microglia + (1 | Donor),
  data = df2
)

emm <- emmeans(model, ~ Region_grouped_microglia)

pairs(emm, adjust = "fdr")

