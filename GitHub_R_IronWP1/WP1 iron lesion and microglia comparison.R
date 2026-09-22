#Aletta
#Script iron WP1 HypIL vs IsoIL Figures 2-3
#Plot first, then statistics


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

Results_combined <- read_excel("Data analysis/Datasets/Results_combined.xlsx")
Results_combined_active <- subset(Results_combined, Region_grouped == "HLA_accu")
Results_combined_select <- subset(Results_combined, Selection == "yes")




# Define desired order
order_levels <- c(
  "NAWM",
  "core_Isointense_Non-foamy",
  "core_Isointense_Foamy",
  "core_Hyperintense_Non-foamy",
  "core_Hyperintense_Foamy",
  "HLA_accu_Isointense_Non-foamy",
  "HLA_accu_Isointense_Foamy",
  "HLA_accu_Hyperintense_Non-foamy",
  "HLA_accu_Hyperintense_Foamy",
  "peri_Isointense_Non-foamy",
  "peri_Isointense_Foamy",
  "peri_Hyperintense_Non-foamy",
  "peri_Hyperintense_Foamy"
)

# Define custom colors for your groups
custom_colors <- c(
  "NAWM" = "#B0B0B0",
  "core_Isointense_Non-foamy" = "#74C0FC",
  "core_Isointense_Foamy" = "#1F78B4",
  "core_Hyperintense_Non-foamy" = "#FFB3C6",
  "core_Hyperintense_Foamy" = "#E31A1C",
  "HLA_accu_Isointense_Non-foamy" = "#74C0FC",
  "HLA_accu_Isointense_Foamy" = "#1F78B4",
  "HLA_accu_Hyperintense_Non-foamy" = "#FFB3C6",
  "HLA_accu_Hyperintense_Foamy" = "#E31A1C",
  "peri_Isointense_Non-foamy" = "#74C0FC",
  "peri_Isointense_Foamy" = "#1F78B4",
  "peri_Hyperintense_Non-foamy" = "#FFB3C6",
  "peri_Hyperintense_Foamy" = "#E31A1C"
)


cols_to_plot <- c(
  "Nonvasc_fibri_per_mm2",
  "Fibri_pary",
  "APP_per_mm2",
  "oxPL_bootstrap_norm")

# Loop over measures and plot separately
for (col in cols_to_plot) {
  p <- ggplot(Results_combined_select, 
              aes(x = factor(Region_grouped_iron_grouped_microglia_NAWM, levels = order_levels), 
                  y = .data[[col]], 
                  fill = Region_grouped_iron_grouped_microglia_NAWM)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "black") +
    scale_fill_manual(values = custom_colors) +
    labs(x = "Region Group", y = col, title = paste("Boxplot of", col)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = c(1.5, 5.5, 9.5), linetype = "dashed", color = "grey40")
  
  print(p)
}




#Now we analyse the cellnumbers

#--------------------------------------------------
# Variables to plot
#--------------------------------------------------

cols_to_plot <- c(
  "GFAP_per_mm2",
  "SOX10_per_mm2",
  "Micro_per_mm2"
)

titles <- c(
  "Astrocytes",
  "Oligodendrocytes",
  "Microglia"
)

#--------------------------------------------------
# Create plotting group
#--------------------------------------------------

Results_combined_select <- Results_combined_select %>%
  mutate(
    Region_grouped_iron_grouped_microglia_NAWM = case_when(
      Region_grouped == "NAWM" ~ "NAWM",
      TRUE ~ paste(
        Region_grouped,
        Iron_pos_neg_microglia,
        sep = "_"
      )
    )
  )

#--------------------------------------------------
# Define plotting order
#--------------------------------------------------

order_levels <- c(
  "NAWM",
  "core_Isointense_Non-foamy",
  "core_Isointense_Foamy",
  "core_Hyperintense_Non-foamy",
  "core_Hyperintense_Foamy",
  "HLA_accu_Isointense_Non-foamy",
  "HLA_accu_Isointense_Foamy",
  "HLA_accu_Hyperintense_Non-foamy",
  "HLA_accu_Hyperintense_Foamy",
  "peri_Isointense_Non-foamy",
  "peri_Isointense_Foamy",
  "peri_Hyperintense_Non-foamy",
  "peri_Hyperintense_Foamy"
)

#--------------------------------------------------
# Colours
#--------------------------------------------------

custom_colors <- c(
  "NAWM" = "#B0B0B0",
  "core_Isointense_Non-foamy" = "#74C0FC",
  "core_Isointense_Foamy" = "#1F78B4",
  "core_Hyperintense_Non-foamy" = "#FFB3C6",
  "core_Hyperintense_Foamy" = "#E31A1C",
  "HLA_accu_Isointense_Non-foamy" = "#74C0FC",
  "HLA_accu_Isointense_Foamy" = "#1F78B4",
  "HLA_accu_Hyperintense_Non-foamy" = "#FFB3C6",
  "HLA_accu_Hyperintense_Foamy" = "#E31A1C",
  "peri_Isointense_Non-foamy" = "#74C0FC",
  "peri_Isointense_Foamy" = "#1F78B4",
  "peri_Hyperintense_Non-foamy" = "#FFB3C6",
  "peri_Hyperintense_Foamy" = "#E31A1C"
)

#--------------------------------------------------
# Plot each cell type
#--------------------------------------------------

for (i in seq_along(cols_to_plot)) {
  
  col <- cols_to_plot[i]
  
  p <- ggplot(
    Results_combined_select,
    aes(
      x = factor(
        Region_grouped_iron_grouped_microglia_NAWM,
        levels = order_levels
      ),
      y = .data[[col]],
      fill = factor(
        Region_grouped_iron_grouped_microglia_NAWM,
        levels = order_levels
      )
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      alpha = 0.8,
      colour = "black"
    ) +
    
    geom_jitter(
      width = 0.2,
      alpha = 0.6,
      size = 1,
      colour = "black"
    ) +
    
    scale_fill_manual(values = custom_colors) +
    
    labs(
      x = "Region Group",
      y = expression("Cell density (cells/mm"^2*")"),
      title = titles[i]
    ) +
    
    theme_minimal(base_size = 10) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    
    geom_vline(
      xintercept = c(1.5, 5.5, 9.5),
      linetype = "dashed",
      colour = "grey40"
    )
  
  print(p)
  
}

#Now FTH+ cell density

cols_to_plot <- c(
  "FTH_GFAP_per_mm2",
  "FTH_SOX10_per_mm2",
  "nr_HLA_FTH_per_mm2"
)

titles <- c(
  "FTH+Astrocytes",
  "FTH+Oligodendrocytes",
  "FTH+Microglia"
)

#--------------------------------------------------
# Create plotting group
#--------------------------------------------------

Results_combined_select <- Results_combined_select %>%
  mutate(
    Region_grouped_iron_grouped_microglia_NAWM = case_when(
      Region_grouped == "NAWM" ~ "NAWM",
      TRUE ~ paste(
        Region_grouped,
        Iron_pos_neg_microglia,
        sep = "_"
      )
    )
  )

#--------------------------------------------------
# Define plotting order
#--------------------------------------------------

order_levels <- c(
  "NAWM",
  "core_Isointense_Non-foamy",
  "core_Isointense_Foamy",
  "core_Hyperintense_Non-foamy",
  "core_Hyperintense_Foamy",
  "HLA_accu_Isointense_Non-foamy",
  "HLA_accu_Isointense_Foamy",
  "HLA_accu_Hyperintense_Non-foamy",
  "HLA_accu_Hyperintense_Foamy",
  "peri_Isointense_Non-foamy",
  "peri_Isointense_Foamy",
  "peri_Hyperintense_Non-foamy",
  "peri_Hyperintense_Foamy"
)

#--------------------------------------------------
# Colours
#--------------------------------------------------

custom_colors <- c(
  "NAWM" = "#B0B0B0",
  "core_Isointense_Non-foamy" = "#74C0FC",
  "core_Isointense_Foamy" = "#1F78B4",
  "core_Hyperintense_Non-foamy" = "#FFB3C6",
  "core_Hyperintense_Foamy" = "#E31A1C",
  "HLA_accu_Isointense_Non-foamy" = "#74C0FC",
  "HLA_accu_Isointense_Foamy" = "#1F78B4",
  "HLA_accu_Hyperintense_Non-foamy" = "#FFB3C6",
  "HLA_accu_Hyperintense_Foamy" = "#E31A1C",
  "peri_Isointense_Non-foamy" = "#74C0FC",
  "peri_Isointense_Foamy" = "#1F78B4",
  "peri_Hyperintense_Non-foamy" = "#FFB3C6",
  "peri_Hyperintense_Foamy" = "#E31A1C"
)

#--------------------------------------------------
# Plot each cell type
#--------------------------------------------------

for (i in seq_along(cols_to_plot)) {
  
  col <- cols_to_plot[i]
  
  p <- ggplot(
    Results_combined_select,
    aes(
      x = factor(
        Region_grouped_iron_grouped_microglia_NAWM,
        levels = order_levels
      ),
      y = .data[[col]],
      fill = factor(
        Region_grouped_iron_grouped_microglia_NAWM,
        levels = order_levels
      )
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      alpha = 0.8,
      colour = "black"
    ) +
    
    geom_jitter(
      width = 0.2,
      alpha = 0.6,
      size = 1,
      colour = "black"
    ) +
    
    scale_fill_manual(values = custom_colors) +
    
    labs(
      x = "Region Group",
      y = expression("Cell density (cells/mm"^2*")"),
      title = titles[i]
    ) +
    
    theme_minimal(base_size = 10) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    
    geom_vline(
      xintercept = c(1.5, 5.5, 9.5),
      linetype = "dashed",
      colour = "grey40"
    )
  
  print(p)
  
}




#Distribution of FTH+ cells

df_long <- Results_combined_select %>%
  select(Region_grouped_iron_grouped_microglia_NAWM,
         Perc_FTH_that_is_iba1,
         Perc_FTH_that_is_GFAP, 
         perc_FTH_that_is_SOX10) %>%
  pivot_longer(cols = starts_with("Perc_FTH_that"), 
               names_to = "Marker", values_to = "Percentage")

# Compute means for each Marker * Region group
means <- df_long %>%
  group_by(Region_grouped_iron_grouped_microglia_NAWM, Marker) %>%
  summarise(mean_val = mean(Percentage, na.rm = TRUE), .groups = "drop")

# Plot with facets by region
ggplot(df_long, aes(x = Percentage, fill = Marker)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = means,
             aes(xintercept = mean_val, color = Marker),
             linetype = "dashed", size = 1) +
  facet_wrap(~ Region_grouped_iron_grouped_microglia_NAWM, scales = "free_y") +
  labs(
    title = "% FTH+ Cells Co-expressing Markers by Region",
    x = "% Positive Cells",
    y = "Density",
    fill = "Marker",
    color = "Mean Line"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "Perc_FTH_that_is_iba1" = "#1f77b4",
    "perc_FTH_that_is_SOX10" = "magenta",
    "Perc_FTH_that_is_GFAP" = "orange"
  )) +
  scale_color_manual(values = c(
    "Perc_FTH_that_is_iba1" = "#1f77b4",
    "perc_FTH_that_is_SOX10" = "magenta",
    "Perc_FTH_that_is_GFAP" = "orange"
  ))


df <- subset(Results_combined_select, Microglia_grouped_NAWM != "NAWM")

# Define pairs
pairs <- list(
  c("Perc_iba1_FTHneg_that_is_p22", "Perc_iba1_FTH_that_is_p22", "p22 (Iba1)"),
  c("Perc_HLA_FTHneg_that_is_iNOS_new", "Perc_HLA_FTH_that_is_iNOS_new", "iNOS (HLA)"),
  c("Perc_iba1_FTHneg_that_is_GAL3", "Perc_iba1_FTH_that_is_GAL3", "GAL3 (Iba1)"),
  c("Perc_FTHneg_HLA_that_is_TSPO", "Perc_FTH_HLA_that_is_TSPO", "TSPO (HLA)"),
  c("percHLAposFTHneg_thatis_PLIN2pos", "percHLAposFTHpos_thatis_PLIN2pos", "PLIN2pos (HLA)"))




for (p in pairs) {
  neg_col <- p[1]
  pos_col <- p[2]
  marker <- p[3]
  
  temp <- df %>%
    select(all_of(c(neg_col, pos_col, "Microglia_grouped_NAWM"))) %>%
    mutate(Sample = row_number()) %>%
    pivot_longer(cols = c(all_of(neg_col), all_of(pos_col)),
                 names_to = "FTH_status", values_to = "Percentage") %>%
    mutate(
      FTH_status = factor(FTH_status, 
                          levels = c(neg_col, pos_col),
                          labels = c("FTH−", "FTH+")),
      Group_FTH = paste(Microglia_grouped_NAWM, FTH_status, sep = "_")
    ) %>%
    mutate(
      Group_FTH = factor(
        Group_FTH,
        levels = c(
          "Non-foamy_FTH−", "Non-foamy_FTH+",
          "Foamy_FTH−", "Foamy_FTH+"
        )
      )
    )
  
  plot <- ggplot(temp, aes(x = Group_FTH, y = Percentage, group = interaction(Sample, Microglia_grouped_NAWM))) +
    geom_boxplot(aes(group = Group_FTH), width = 0.8, alpha = 0.2, color = "black") +
    geom_line(aes(color = Microglia_grouped_NAWM), alpha = 0.7) +
    geom_point(aes(fill = Microglia_grouped_NAWM), shape = 21, color = "black", size = 2, alpha = 0.6) +
    scale_fill_manual(values = c("Non-foamy" = "#1f77b4", "Foamy" = "#ff7f0e")) +
    scale_color_manual(values = c("Non-foamy" = "#1f77b4", "Foamy" = "#ff7f0e")) +
    scale_y_continuous(limits = c(0,100))+
    labs(
      title = marker,
      x = "",
      y = "% Positive cells",
      fill = "Microglia group",
      color = "Microglia group"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(plot)
}


Results_LAMP1_GAL3_summary <- read_excel("Data analysis/Datasets/Results_LAMP1_GAL3_summary.xlsx")

# Variables to plot
cols_to_plot_LAMP1 <- c(
  "Lysosomal_load",
  "Perc_damaged_lysosomes"
)

# Define the order (in the sequence you want them on the x-axis)
order_levels_LAMP1 <- c(
  "Isointense_non-foamy",
  "Isointense_foamy",
  "Hyperintense_non-foamy",
  "Hyperintense_foamy"
)

# Define matching colors (names MUST match your data exactly)
custom_colors_LAMP1 <- c(
  "Isointense_non-foamy" = "#74C0FC",  # light blue
  "Isointense_foamy"    = "#1F78B4",   # dark blue
  "Hyperintense_non-foamy" = "#FFB3C6", # light pink
  "Hyperintense_foamy"  = "#E31A1C"    # red
)

# Ensure the factor levels are set correctly
Results_LAMP1_GAL3_summary$Iron_microglia <- factor(
  Results_LAMP1_GAL3_summary$Iron_microglia,
  levels = order_levels_LAMP1
)

# Loop over your variables
for (col in cols_to_plot_LAMP1) {
  p <- ggplot(Results_LAMP1_GAL3_summary,
              aes(x = Iron_microglia,
                  y = .data[[col]],
                  fill = Iron_microglia)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "black") +
    scale_fill_manual(values = custom_colors_LAMP1) +
    labs(x = "Iron / Microglia Group",
         y = col,
         title = paste("Boxplot of", col)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
}












# And now statistics
############################################################
## Statistics
## Iron lesion manuscript

############################################################

#===========================================================

# Packages

#===========================================================


library(dplyr)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(broom)
library(openxlsx)
library(readxl)





#===========================================================
# Helper functions
#===========================================================

#-----------------------------------------------------------
# Order of lesion groups
#-----------------------------------------------------------

get_region_levels <- function() {
  
  c(
    "NAWM",
    "core_Isointense_Non-foamy",
    "core_Isointense_Foamy",
    "core_Hyperintense_Non-foamy",
    "core_Hyperintense_Foamy",
    "HLA_accu_Isointense_Non-foamy",
    "HLA_accu_Isointense_Foamy",
    "HLA_accu_Hyperintense_Non-foamy",
    "HLA_accu_Hyperintense_Foamy",
    "peri_Isointense_Non-foamy",
    "peri_Isointense_Foamy",
    "peri_Hyperintense_Non-foamy",
    "peri_Hyperintense_Foamy"
  )
  
}

#-----------------------------------------------------------
# Create plotting/statistics group
#-----------------------------------------------------------

create_region_group <- function(data){
  
  data %>%
    mutate(
      
      Plot_group = case_when(
        
        Region_grouped == "NAWM" ~ "NAWM",
        
        TRUE ~ paste(
          Region_grouped,
          Iron_pos_neg_microglia,
          sep = "_"
        )
        
      ),
      
      Plot_group = factor(
        Plot_group,
        levels = get_region_levels()
      )
      
    )
  
}

#-----------------------------------------------------------
# Check inputs
#-----------------------------------------------------------

check_region_inputs <- function(data,
                                variables,
                                donor){
  
  required <- c(
    "Region_grouped",
    "Iron_pos_neg_microglia",
    donor
  )
  
  missing <- required[!required %in% names(data)]
  
  if(length(missing) > 0){
    
    stop(
      paste(
        "Missing columns:",
        paste(missing, collapse = ", ")
      )
    )
    
  }
  
  bad_vars <- variables[!variables %in% names(data)]
  
  if(length(bad_vars) > 0){
    
    stop(
      paste(
        "These variables are not present:",
        paste(bad_vars, collapse = ", ")
      )
    )
    
  }
  
}

#===========================================================
# Count variable -> offset area lookup
#===========================================================

count_area_lookup <- c(
  
  #---------------------------------------------------------
  # Cell types
  #---------------------------------------------------------
  
  Num_GFAP = "Area_GFAP",
  Num_FTH_GFAP = "Area_GFAP",
  
  Num_SOX10 = "Area_SOX10",
  num_SOX10_FTH = "Area_SOX10",
  
  Num_micro_total = "Area_micro",
  Num_active_micro_total = "Area_micro",
  
  #---------------------------------------------------------
  # APP
  #---------------------------------------------------------
  
  APP_events = "APP_area",
  
  #---------------------------------------------------------
  # TBB
  #---------------------------------------------------------
  
  Num_TBB = "Area_TBB",
  
  #---------------------------------------------------------
  # iNOS / HLA
  #---------------------------------------------------------
  
  Num_HLA = "Area_iNOS",
  Num_iNOS = "Area_iNOS",
  Num_FTH = "Area_iNOS",
  
  Num_HLA_FTH = "Area_iNOS",
  Num_HLA_FTHneg = "Area_iNOS",
  
  Num_HLA_iNOS = "Area_iNOS",
  Num_HLA_iNOS_FTH = "Area_iNOS",
  
  Num_FTH_iNOS = "Area_iNOS",
  
  #---------------------------------------------------------
  # GAL3
  #---------------------------------------------------------
  
  Num_iba1_GAL3 = "Area_GAL3",
  
  Num_iba1_FTH = "Area_GAL3",
  Num_iba1_FTHneg = "Area_GAL3",
  
  #---------------------------------------------------------
  # TSPO
  #---------------------------------------------------------
  
  nr_HLA = "Area_TSPO_quant",
  nr_FTH = "Area_TSPO_quant",
  nr_TSPO = "Area_TSPO_quant",
  
  nr_HLA_FTH = "Area_TSPO_quant",
  nr_HLA_TSPO = "Area_TSPO_quant",
  nr_HLA_FTH_TSPO = "Area_TSPO_quant",
  
  #---------------------------------------------------------
  # p22
  #---------------------------------------------------------
  
  Nr_iba1 = "Area_p22",
  Nr_FTH = "Area_p22",
  Nr_p22 = "Area_p22",
  
  Nr_iba1_FTH = "Area_p22",
  Nr_iba1_p22 = "Area_p22",
  Nr_FTH_p22 = "Area_p22",
  Nr_iba1_FTH_p22 = "Area_p22",
  
  #---------------------------------------------------------
  # PLIN2
  #---------------------------------------------------------
  
  Nr_HLA = "Area_PLIN2",
  Nr_FTH = "Area_PLIN2",
  Nr_PLIN2 = "Area_PLIN2",
  
  Nr_HLA_FTH = "Area_PLIN2",
  Nr_HLA_PLIN2 = "Area_PLIN2",
  Nr_FTH_PLIN2 = "Area_PLIN2",
  Nr_HLA_FTH_PLIN2 = "Area_PLIN2",
  
  #---------------------------------------------------------
  # fibrinogen
  #---------------------------------------------------------
  
  FibriPos_UlexNeg = "Area_fibri"
  
)


#===========================================================
# State count variable -> offset area lookup
#===========================================================

phenotype_lookup <- list(
  
  p22 = list(
    pos_success = "Nr_iba1_FTH_p22",
    pos_total   = "Nr_iba1_FTH",
    neg_success = "Nr_iba1_FTHneg_p22",
    neg_total   = "Nr_iba1_FTHneg"
  ),
  
  GAL3 = list(
    pos_success = "Num_iba1_FTH_GAL3high",
    pos_total   = "Num_iba1_FTH",
    neg_success = "Num_iba1_FTHneg_GAL3high",
    neg_total   = "Num_iba1_FTHneg"
  ),
  
  iNOS = list(
    pos_success = "Num_HLA_iNOS_FTH_new",
    pos_total   = "Num_HLA_FTH_new",
    neg_success = "Num_HLA_iNOS_FTHneg_new",
    neg_total   = "Num_HLA_FTHneg_new"
  ),
  
  PLIN2 = list(
    pos_success = "HLAposFTHposPLIN2pos",
    pos_total = "HLAposFTHpos",
    neg_success = "HLAposFTHnegPLIN2pos",
    neg_total = "HLAposFTHneg"
  ),
  
  TSPO = list(
    pos_success = "nr_HLA_FTH_TSPO",
    pos_total   = "nr_HLA_FTH",
    neg_success = "nr_HLA_FTHneg_TSPO",
    neg_total   = "nr_HLA_FTHneg"
  )
  
)




#===========================================================
# Main mixed model function
#===========================================================

run_region_lmm <- function(data,
                           variables,
                           donor = "Donor_sample"){
  
  #---------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------
  
  check_region_inputs(
    data = data,
    variables = unname(variables),
    donor = donor
  )
  
  #---------------------------------------------------------
  # Create Plot_group
  #---------------------------------------------------------
  
  data <- create_region_group(data)
  
  #---------------------------------------------------------
  # Containers
  #---------------------------------------------------------
  
  model_list <- list()
  
  anova_table <- tibble()
  
  emmeans_table <- tibble()
  
  contrast_table <- tibble()
  
  #---------------------------------------------------------
  # Loop over outcomes
  #---------------------------------------------------------
  
  for(i in seq_along(variables)){
    
    outcome_label <- names(variables)[i]
    
    outcome <- unname(variables[i])
    
    message("Running: ", outcome_label)
    
    #-------------------------------------------------------
    # Build formula
    #-------------------------------------------------------
    
    form <- as.formula(
      
      paste0(
        outcome,
        " ~ Plot_group + (1|",
        donor,
        ")"
      )
      
    )
    
    #-------------------------------------------------------
    # Fit model
    #-------------------------------------------------------
    
    model <- lmer(
      form,
      data = data,
      REML = FALSE
    )
    
    model_list[[outcome_label]] <- model
    
    #-------------------------------------------------------
    # ANOVA
    #-------------------------------------------------------
    
    tmp_anova <- anova(model) |>
      tibble::rownames_to_column("Effect") |>
      mutate(
        Outcome = outcome_label,
        .before = 1
      )
    
    anova_table <- bind_rows(
      anova_table,
      tmp_anova
    )
    
    #-------------------------------------------------------
    # Estimated marginal means
    #-------------------------------------------------------
    
    emm <- emmeans(
      model,
      ~ Plot_group
    )
    
    tmp_emm <- as.data.frame(emm) |>
      mutate(
        Outcome = outcome_label,
        .before = 1
      )
    
    emmeans_table <- bind_rows(
      emmeans_table,
      tmp_emm
    )
    
    #-------------------------------------------------------
    # All pairwise contrasts
    #-------------------------------------------------------
    
    #-------------------------------------------------------
    # All pairwise contrasts (UNADJUSTED)
    #-------------------------------------------------------
    
    tmp_contrasts <- pairs(
      emm,
      adjust = "none"
    ) |>
      as.data.frame() |>
      mutate(
        Outcome = outcome_label,
        .before = 1
      )
    
    contrast_table <- bind_rows(
      contrast_table,
      tmp_contrasts
    )
    
  }
  
  #---------------------------------------------------------
  # Return results
  #---------------------------------------------------------
  
  results <- list(
    
    models = model_list,
    
    anova = anova_table,
    
    emmeans = emmeans_table,
    
    contrasts_all = contrast_table,
    
    contrasts = filter_region_contrasts(contrast_table)
    
  )
  
  return(results)
  
}


run_region_glmm <- function(data,
                            variables,
                            donor = "Donor_sample"){
  
  #---------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------
  
  check_region_inputs(
    data = data,
    variables = unname(variables),
    donor = donor
  )
  
  #---------------------------------------------------------
  # Create Plot_group
  #---------------------------------------------------------
  
  data <- create_region_group(data)
  
  #---------------------------------------------------------
  # Containers
  #---------------------------------------------------------
  
  model_list <- list()
  
  anova_table <- tibble()
  
  emmeans_table <- tibble()
  
  contrast_table <- tibble()
  
  #---------------------------------------------------------
  # Loop over outcomes
  #---------------------------------------------------------
  
  for(i in seq_along(variables)){
    
    outcome_label <- names(variables)[i]
    
    outcome <- unname(variables[i])
    
    message("Running: ", outcome_label)
    
    #-------------------------------------------------------
    # Build formula
    #-------------------------------------------------------
    
    area <- count_area_lookup[[outcome]]
    
    if(is.null(area)){
      stop(paste("No area defined for", outcome))
    }
    
    form <- as.formula(
      
      paste0(
        outcome,
        " ~ Plot_group + offset(log(",
        area,
        ")) + (1|",
        donor,
        ")"
      )
      
    )
    
    #-------------------------------------------------------
    # Fit model
    #-------------------------------------------------------
    
    model <- glmmTMB(
      form,
      data = data,
      family = nbinom2()
    )
    
    model_list[[outcome_label]] <- model
    
    #-------------------------------------------------------
    # Likelihood ratio test
    #-------------------------------------------------------
    
    null_model <- glmmTMB(
      update(
        form,
        . ~ . - Plot_group
      ),
      data = data,
      family = nbinom2()
    )
    
    lrt <- anova(
      null_model,
      model
    )
    
    tmp_anova <- tibble(
      
      Outcome = outcome_label,
      
      Effect = "Plot_group",
      
      Chisq = lrt$Chisq[2],
      
      Df = lrt$Df[2],
      
      `Pr(>Chisq)` = lrt$`Pr(>Chisq)`[2]
      
    )
    
    anova_table <- bind_rows(
      anova_table,
      tmp_anova
    )
    
    #-------------------------------------------------------
    # Estimated marginal means
    #-------------------------------------------------------
    
    emm <- emmeans(
      model,
      ~ Plot_group
    )
    
    tmp_emm <- summary(
      emm,
      type = "response"
    ) |>
      as.data.frame() |>
      mutate(
        Outcome = outcome_label,
        .before = 1
      )
    
    emmeans_table <- bind_rows(
      emmeans_table,
      tmp_emm
    )
    
    #-------------------------------------------------------
    # All pairwise contrasts
    #-------------------------------------------------------
    
    #-------------------------------------------------------
    # All pairwise contrasts (UNADJUSTED)
    #-------------------------------------------------------
    
    tmp_contrasts <- pairs(
      emm,
      adjust = "none"
    ) |>
      as.data.frame() |>
      mutate(
        Outcome = outcome_label,
        .before = 1
      )
    
    contrast_table <- bind_rows(
      contrast_table,
      tmp_contrasts
    )
    
  }
  
  #---------------------------------------------------------
  # Return results
  #---------------------------------------------------------
  
  results <- list(
    
    models = model_list,
    
    anova = anova_table,
    
    emmeans = emmeans_table,
    
    contrasts_all = contrast_table,
    
    contrasts = filter_region_contrasts(contrast_table)
    
  )
  
  return(results)
  
}

run_phenotype_glmm <- function(data,
                               markers,
                               donor = "Donor",
                               lesion = "Lesion_ID"){
  
  #---------------------------------------------------------
  # Containers
  #---------------------------------------------------------
  
  model_list <- list()
  
  anova_table <- tibble()
  
  emmeans_table <- tibble()
  
  contrast_table <- tibble()
  
  #---------------------------------------------------------
  # Remove NAWM
  #---------------------------------------------------------
  
  data <- subset(
    data,
    Microglia_grouped_NAWM != "NAWM"
  )
  
  #---------------------------------------------------------
  # Loop over markers
  #---------------------------------------------------------
  
  for(i in seq_along(markers)){
    
    marker_label <- markers[i]
    
    message("Running: ", marker_label)
    
    #-------------------------------------------------------
    # Get lookup
    #-------------------------------------------------------
    
    lookup <- phenotype_lookup[[marker_label]]
    
    #-------------------------------------------------------
    # Create long dataset
    #-------------------------------------------------------
    
    temp <- data |>
      
      mutate(
        
        FTH_pos_success = .data[[lookup$pos_success]],
        FTH_pos_total   = .data[[lookup$pos_total]],
        
        FTH_neg_success = .data[[lookup$neg_success]],
        FTH_neg_total   = .data[[lookup$neg_total]]
        
      ) |>
      
      tidyr::pivot_longer(
        
        cols = c(
          FTH_pos_success,
          FTH_pos_total,
          FTH_neg_success,
          FTH_neg_total
        ),
        
        names_to = c("FTH_status",".value"),
        
        names_pattern = "FTH_(pos|neg)_(.*)"
        
      ) |>
      
      mutate(
        
        FTH_status = factor(
          FTH_status,
          levels = c("neg","pos"),
          labels = c("FTH−","FTH+")
        ),
        
        Microglia_grouped_NAWM = factor(
          Microglia_grouped_NAWM,
          levels = c(
            "Non-foamy",
            "Foamy"
          )
        )
        
      )
    
    #-------------------------------------------------------
    # Fit model
    #-------------------------------------------------------
    
    model <- glmmTMB(
      
      cbind(
        success,
        total-success
      ) ~
        
        FTH_status *
        Microglia_grouped_NAWM +
        
        (1|Donor/Lesion_ID),
      
      family = betabinomial(),
      
      data = temp
      
    )
    
    model_list[[marker_label]] <- model
    
    #-------------------------------------------------------
    # Likelihood ratio test
    #-------------------------------------------------------
    
    null_model <- update(
      model,
      . ~ . -
        FTH_status *
        Microglia_grouped_NAWM
    )
    
    lrt <- anova(
      null_model,
      model
    )
    
    tmp_anova <- tibble(
      
      Marker = marker_label,
      
      Effect = "Interaction",
      
      Chisq = lrt$Chisq[2],
      
      Df = lrt$Df[2],
      
      `Pr(>Chisq)` = lrt$`Pr(>Chisq)`[2]
      
    )
    
    anova_table <- bind_rows(
      anova_table,
      tmp_anova
    )
    
    #-------------------------------------------------------
    # Estimated marginal means
    #-------------------------------------------------------
    
    emm <- emmeans(
      model,
      ~ FTH_status * Microglia_grouped_NAWM
    )
    
    tmp_emm <- summary(
      emm,
      type = "response"
    ) |>
      
      as.data.frame() |>
      
      mutate(
        Marker = marker_label,
        .before = 1
      )
    
    emmeans_table <- bind_rows(
      emmeans_table,
      tmp_emm
    )
    
    #-------------------------------------------------------
    # Pairwise contrasts
    #-------------------------------------------------------
    
    tmp_contrasts <- pairs(
      emm,
      adjust = "none"
    ) |>
      
      as.data.frame() |>
      
      mutate(
        Marker = marker_label,
        .before = 1
      )
    
    contrast_table <- bind_rows(
      contrast_table,
      tmp_contrasts
    )
    
  }
  
  #---------------------------------------------------------
  # Return results
  #---------------------------------------------------------
  
  results <- list(
    
    models = model_list,
    
    anova = anova_table,
    
    emmeans = emmeans_table,
    
    contrasts = contrast_table
    
  )
  
  
  return(results)
  
}

#===========================================================
# Keep planned contrasts and apply FDR correction
#===========================================================

filter_region_contrasts <- function(contrast_table){
  
  planned <- contrast_table %>%
    
    mutate(
      contrast = gsub("[()]", "", contrast)
    ) %>%
    
    filter(
      
      grepl("^NAWM - ", contrast) |
        
        (grepl("^core_", contrast) &
           grepl(" - core_", contrast)) |
        
        (grepl("^HLA_accu_", contrast) &
           grepl(" - HLA_accu_", contrast)) |
        
        (grepl("^peri_", contrast) &
           grepl(" - peri_", contrast))
      
    ) %>%
    
    group_by(Outcome) %>%
    
    mutate(
      p.value = p.adjust(
        p.value,
        method = "fdr"
      )
    ) %>%
    
    ungroup()
  
  planned
  
}



Results_combined <- read_excel("T:/Lab Absinta/15_ALETTA/Iron/WP1/Data analysis/Datasets/Results_combined.xlsx")
Results_combined_active <- subset(Results_combined, Region_grouped == "HLA_accu")
Results_combined_select <- subset(Results_combined, Selection == "yes")
Results_combined_select_active <- subset(Results_combined_select, Region_grouped == "HLA_accu"|Region_grouped == "NAWM")

stats <- run_region_lmm(
  data = Results_combined_select,
  variables = c(
    Fibrinogen_parenchyma = "Fibri_pary",
    oxPL = "oxPL_bootstrap_norm",
    GAL3_semi = "GAL3_semi_quant"
  )
)

print(stats$contrasts, n = Inf)

stats_cells <- run_region_glmm(
  data = Results_combined_select,
  variables = c(
    Fibrinogen = "FibriPos_UlexNeg",
    APP = "APP_events",
    Astrocytes = "Num_GFAP",
    Oligodendrocytes = "Num_SOX10",
    Microglia = "Nr_iba1",
    Astrocytes_FTH = "Num_FTH_GFAP",
    Oligos_FTH = "num_SOX10_FTH",
    Microglia_FTH = "Nr_iba1_FTH"
  )
)

print(stats_cells$contrasts, n = Inf)

stats_phenotype <- run_phenotype_glmm(
  data = Results_combined_select,
  markers = c(
    "p22",
    "GAL3",
    "iNOS",
    "TSPO",
    "PLIN2"
  )
)

print(stats_phenotype$contrasts, n = Inf)





Results_LAMP1_GAL3_summary <- read_excel("Data analysis/Datasets/Results_LAMP1_GAL3_summary.xlsx")

# Variables to plot
cols_to_plot_LAMP1 <- c(
  "Lysosomal_load",
  "Perc_damaged_lysosomes"
)

# Define the order (in the sequence you want them on the x-axis)
order_levels_LAMP1 <- c(
  "Isointense_non-foamy",
  "Isointense_foamy",
  "Hyperintense_non-foamy",
  "Hyperintense_foamy"
)

# Define matching colors (names MUST match your data exactly)
custom_colors_LAMP1 <- c(
  "Isointense_non-foamy" = "#74C0FC",  # light blue
  "Isointense_foamy"    = "#1F78B4",   # dark blue
  "Hyperintense_non-foamy" = "#FFB3C6", # light pink
  "Hyperintense_foamy"  = "#E31A1C"    # red
)

# Ensure the factor levels are set correctly
Results_LAMP1_GAL3_summary$Iron_microglia <- factor(
  Results_LAMP1_GAL3_summary$Iron_microglia,
  levels = order_levels_LAMP1
)

# Loop over your variables
for (col in cols_to_plot_LAMP1) {
  p <- ggplot(Results_LAMP1_GAL3_summary,
              aes(x = Iron_microglia,
                  y = .data[[col]],
                  fill = Iron_microglia)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "black") +
    scale_fill_manual(values = custom_colors_LAMP1) +
    labs(x = "Iron / Microglia Group",
         y = col,
         title = paste("Boxplot of", col)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
}


LAMP1_model <- lmer(
  Lysosomal_load ~ Iron_microglia + (1 | Donor),
  data = Results_LAMP1_GAL3_summary
)
emmeans(LAMP1_model, pairwise ~ Iron_microglia, adjust = "fdr")

GAL3_model <- lmer(
  Perc_damaged_lysosomes ~ Iron_microglia + (1 | Donor), 
  data = Results_LAMP1_GAL3_summary)
emmeans(GAL3_model, pairwise ~ Iron_microglia, adjust = "fdr")

