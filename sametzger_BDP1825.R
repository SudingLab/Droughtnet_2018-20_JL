#' ---
#' title: "Updates to Larson et al. -- Analyses for 2018-2025 Boulder Drought Project"
#' author: "Sam Metzger"
#' genesis date: "04 march 2026"
#' ---
#'

#' 
#' **Welcome!**
#'
#' This R script contains code to conduct analyses and create figures adapted from 
#' in the following manuscript in preparation:
#' 
#' Larson, J., B. Anacker, T. Merchant, K. Suding. (In Prep.) Experimental 
#' grazing and rainfall treatments reveal key contingencies underlying 
#' rangeland resistance
#' 
#' This code continues analysis from the first three years (2018-2020) to the present (2025)
#' of an ongoing field experiment manipulating rainfall (dry, wet, control) 
#' and grazing (growing season, dormant season, ungrazed) in a grassland (Boulder, CO, USA).
#' 
#' This ongoing project is conducted in collaboration with the *University of Colorado* 
#' and *City of Boulder Open Space and Mountain Parks (Boulder, CO, USA)*.
#' 

#' 
#' **Load packages, functions, and options**
#' 


# Packages
#+ results=FALSE, message=FALSE, warning=FALSE
#install.packages("...")
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(performance)
library(scales)
library(codyn)
library(patchwork)
library(magick)
library(ggeffects)
library(emmeans)
library(ggExtra)
library(car)
library(GGally)
library(ggcorrplot)
library(ggpubr)
library(cowplot)
#library(Hmisc)
#library(gridExtra)

# Load source code for Oridicenter() function
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')
# Ordihull() source code
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:orglhull?do=export_code&codeblock=1')

# Load function to calculate standard errors
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#' 
#' **Load data**
#'
#'set path
wdpath <- "/Users/sametzger/Library/CloudStorage/OneDrive-SharedLibraries-UCB-O365/Thomas Merchant - Plant_Map_Team/"

#'
#'   **Plant Community Composition**
#'              


#'  **Overview**
#'  
#'  This script explores ...
#'  
#'  
#'  

#' 
#' **Read in data**
#' 

comp_path <-  paste0(wdpath, "Data/plant_area_clipped_18_25.csv")
comp_raw <- read.csv(comp_path)

#' View dataset structure:
str(comp_raw)               

#' Remove non-plant objects, unknown species that could not be uniquely identified,
#' and early-phenology species that could not be reliably detected 
exclude_species <- c("cow", "rock",   # non-plant objects, 
                     "dr", "lo", "castilleja_sp", "mertensia_sp",     # early-phenology species
                     "uk_g", "uk_f", "uk f", # unknown grasses/forbs
                     "uk", "ukf", "ukg", "unk", "f", "uk1")  
comp <- comp_raw %>%
  mutate(across(c(year, block, plot, species_clean), as.factor)) %>%
  filter(!species_clean %in% exclude_species) %>% 
  filter(!is.na(species_clean))

#'check this has worked 
unique(comp$species_clean)

###
#' Recode to pool species that were challenging to differentiate but functionally-similar
levels(comp$species_clean)[levels(comp$species_clean)=="ch"] <- "cg" # Chondrosum (syn. Bouteloua) gracilis and C. hirsuta
levels(comp$species_clean)[levels(comp$species_clean)=="pl"] <- "td" # Podospermum lacinatum and Tragopogon dubius

###

#' Incorporate 'present' species?
#' species marked as present but with no cover data
#' dont have prcessed yet for 2023-2025

#' Create Wide matrix of plot- and species-level aerial cover with
#' summed aerial cover for each species within each 0.5m2 plot and year.
comp_wide <- comp %>%
  filter(!species_clean %in% exclude_species) %>%
  group_by(year, block, plot, species_clean) %>%
  summarise(cover = sum(area_clipped, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species_clean, values_from = cover, values_fill = 0)

#write.csv(comp_wide, 'Data/plant_comp_wide_18_25.csv')

#' Revised long matrix of total species cover
comp_long <- comp_wide %>%
  pivot_longer(cols = -c(year, block, plot), names_to = "species clean", values_to = "cover")

#write.csv(comp_long, 'Data/plant_comp_long_18_25.csv')

#' Create info dataframe sorted by year, block, plot and 
#' year and block as factors + reference levels for fixed effects
comp_wide_info <- comp %>%
  distinct(year, block, plot, rain_trt, graze_trt) %>%
  arrange(year, block, plot) %>%
  mutate(block = as.factor(block),
         year = factor(year, levels = c("2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")),
         rain_trt = factor(rain_trt, levels = c("D", "W", "C")),
         graze_trt = factor(graze_trt, levels = c("FG", "SG", "C")),
         b_p = as.factor(paste(block, plot, sep = "_")))

#' Create a species-only matrix 
comp_spp <- comp_wide %>%
  ungroup() %>%
  select(-year, -block, -plot)

#' sqrt transform
comp_spp_sqrt <- sqrt(comp_spp)

#' Convert absolute abundances to relative abundances
spp_rel      <- decostand(comp_spp, "total")        # relative abundance
spp_rel_sqrt <- decostand(sqrt(comp_spp), "total")  # sqrt + relative

#' Function for relative abundances 
summarise_rel_abun <- function(df, col_name, info = comp_wide_info) {
  df %>%
    bind_cols(info %>% select(year)) %>%  # attach year for grouping
    group_by(year) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-year, names_to = "species", values_to = col_name) %>%
    mutate(across(all_of(col_name), ~ . * 100))
}

#' Apply once for each transformation
spp_mean_rel <- summarise_rel_abun(spp_rel, "rel_abun")
comp_mean_rel_sqrt <- summarise_rel_abun(spp_rel_sqrt, "sqrt_rel_abun")

#' long-format version of community data with sqrt-transformed relative abundances, 
#' with site/treatment info attached.
comp_dat_long <- comp_wide %>%
  # transform species matrix BEFORE pivoting
  mutate(across(-c(year, block, plot), ~ sqrt(.))) %>%  # sqrt transform
  {bind_cols(select(., year, block, plot), decostand(select(., -year, -block, -plot), "total")
  )} %>% # decostand on whole matrix
  # recombine treatment info
  left_join(comp_wide_info, by = c("year", "block", "plot")) %>%
  mutate(b_p = as.factor(paste(block, plot, sep = "_")),
    year = factor(year, levels = c("2018","2019","2020","2021", "2022","2023","2024","2025")),
    year_num = as.integer(as.character(year))) %>%
  pivot_longer(cols = -c(year, year_num, block, plot, rain_trt, graze_trt, b_p),
               names_to = "species", values_to = "rel_cov") %>%
  mutate(rel_cov = replace_na(rel_cov, 0))

#' Verify transformation worked correctly
comp_dat_long %>%
  filter(plot == 1, block == 1, year == "2018") %>%
  filter(rel_cov > 0) %>%
  select(species, rel_cov) %>%
  summarise(total = sum(rel_cov))  # should sum to ~1

#' 
#' *Gathering Weather Data From NLR Station*
#' **Load data**
#'

#' read in files
#' filter out NA and trace codes
#' make columns for experimental plots (wet and dry treatments)
precipData <- read.csv(paste0(wdpath, "Data/BoulderDailyWeather.csv")) %>% 
  rename(year = Year, month = Month, day = Day) %>% 
  filter(precip != -998.00, precip != -999.00) %>% 
  mutate(month = as.numeric(unlist(month)),
         precip_cm = as.numeric(unlist(precip)) * 2.54) %>% 
  mutate(drought_cm = if_else(month %in% 5:9, precip_cm * 0.67, NA_real_),
         wet_cm = if_else(month %in% 5:9, precip_cm * 1.67, NA_real_)) 

#' match calculated experimental rain treatments to rain_trt codes
#' focus on experimental period May - September
precipData_long <- precipData %>%
  select(-precip) %>% # Drops the original raw column
  pivot_longer(cols = c(precip_cm, drought_cm, wet_cm),
               names_to = "rain_trt", values_to = "precip") %>%
  mutate(rain_trt = case_match(rain_trt, "precip_cm" ~ "C", "drought_cm" ~ "D", "wet_cm" ~ "W", .default = rain_trt)) %>% 
  filter(month %in% 5:9) 

#' sum total of precip per growing year
precip_annual <- precipData_long %>%
  filter(year %in% 2018:2025) %>% 
  filter(month %in% 5:9) %>% 
  group_by(year, rain_trt) %>%
  summarise(total_precip = sum(precip, na.rm = TRUE), .groups = "drop") %>%
  mutate(rain_trt = factor(rain_trt, levels = c("D", "C", "W")))


#'
#' *Estimate composition metrics*
#' *1. Richness, Evenness, & Shannon Diversity* 
#' 

#' *Step 1: calculate diversity metrics*
structure_sqrt <- community_structure(comp_dat_long, time.var = "year", replicate.var = "b_p", abundance.var = "rel_cov")
div_sqrt <- community_diversity(comp_dat_long, time.var = "year", replicate.var = "b_p", abundance.var = "rel_cov")

#' Step 2: combine + join treatments
structure_raw_sqrt <- structure_sqrt %>%
  mutate(shan = div_sqrt$Shannon) %>%
  rename(rich = richness, even = Evar) %>%
  left_join(comp_wide_info, by = c("b_p", "year"))

#' Step 3: scale for analysis
structure_all_sqrt <- structure_raw_sqrt %>%
  ungroup() %>%
  mutate(across(where(is.numeric), scale),
         b_p = as.factor(b_p))

#' Step 4: format for plotting 
structure_sum_sqrt <- structure_raw_sqrt %>%
  select(-b_p, -block) %>%
  pivot_longer(c(rich, even, shan), names_to = "metric", values_to = "value") %>%
  group_by(rain_trt, graze_trt, year, metric) %>%
  summarise(mean = mean(value), se = se(value), .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year)),
         rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
         graze_trt = fct_recode(graze_trt, "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
         fct_relevel("Ungrazed", "Fall Graze", "Spring Graze"),
         metric = factor(metric, levels = c("rich", "even", "shan"), labels = c("Richness", "Evenness", "Diversity")))

#' Format dataframes for analysis
#' Create several sets of unique treatment contrasts for models
#' 
#' *artefacts from Julie - probably not using this stepwise contrast method in 2026*
#' Function for creating contrasts
#set_contrasts <- function(df, year_ref, rain_ref, graze_ref) {
 # df %>% mutate(year = fct_relevel(year, year_ref),
  #              rain_trt  = fct_relevel(rain_trt, rain_ref),
   #             graze_trt = fct_relevel(graze_trt, graze_ref))
#  }

#' Control rain as reference ("does drought/wet differ from ambient?")
  #structure_con <- set_contrasts(structure_all_sqrt, "2018", "C", "C")
#' Drought as reference ("does control/wet differ from drought?")
  #structure_dry <- set_contrasts(structure_all_sqrt, "2018", "D", "C")
#' Wet as reference ("does control/drought differ from wet?")
  #structure_wet <- set_contrasts(structure_all_sqrt, "2018", "W", "C")
  
#' join precipitation into composition dataframe
structure_precip <- structure_raw_sqrt %>%
  mutate(year_num = as.numeric(as.character(year))) %>%
  left_join(precip_annual, by = c("rain_trt", "year_num" = "year")) %>% 
  mutate(rain_trt = factor(rain_trt, levels = c("D", "C", "W")))

#'
#' *2. Functional Groups* 
#' 

#' Starting dataframe
species_info <- comp_raw%>%
  dplyr::select(species_clean, genus, spp, origin, growth_form, life_history, photo_path, fun_grp)%>%
  distinct(species_clean, .keep_all = T)

#' Step 1: join functional group info
comp_dat_long_fungrp <- comp_dat_long %>%
    left_join(species_info %>% select(species_clean, fun_grp), 
              by = c("species" = "species_clean")) %>%
    filter(!fun_grp %in% c("G_NA", "F_NA", "S", NA))
  
#' Step 2: summarize plot-level cover by functional group - all years at once
  fungrp_df <- comp_dat_long_fungrp %>%
    group_by(year, b_p, block, rain_trt, graze_trt, fun_grp) %>%
    summarise(fungrp_cov = sum(rel_cov, na.rm = TRUE), .groups = "drop")
  
#' Step 3: save 2018 baseline covers - all functional groups at once
  fungrp_baseline <- fungrp_df %>%
    filter(year == "2018") %>%
    ungroup() %>%
    mutate(across(where(is.numeric), scale)) %>%
    select(b_p, fun_grp, fungrp_cov) %>%
    pivot_wider(names_from = fun_grp, values_from = fungrp_cov)
  
#' Step 4: format for figures
  fungrp_sum <- fungrp_df %>%
    select(-b_p, -block) %>%
    group_by(year, rain_trt, graze_trt, fun_grp) %>%
    summarise(mean = mean(fungrp_cov), se = se(fungrp_cov), .groups = "drop") %>%
    mutate(year = as.numeric(as.character(year)),
           rain_trt = factor(rain_trt, levels = c("C", "D", "W")),
           graze_trt = fct_recode(graze_trt, "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
           fct_relevel("Ungrazed", "Fall Graze", "Spring Graze")) %>%
    rename(metric = fun_grp)
  
#' Step 5: contrast sets using the helper function from before - applied to fungrp_df directly
  fungrp_con <- set_contrasts(fungrp_df, "2018", "C", "C")  # control references
  fungrp_dry <- set_contrasts(fungrp_df, "2018", "D", "C")  # drought reference
  fungrp_wet <- set_contrasts(fungrp_df, "2018", "W", "C")  # wet reference
  
#' Step 6: filter and scale by functional group using a function
  scale_fungrp <- function(df, grp) {
    df %>%
      filter(fun_grp == grp) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), scale))
  }
  
#' Apply cleanly 
  g_c4 <- scale_fungrp(fungrp_con, "G_C4")
  g_c3 <- scale_fungrp(fungrp_con, "G_C3")
  ann  <- scale_fungrp(fungrp_con, "G_A-BI")
  f    <- scale_fungrp(fungrp_con, "F")
  
#' 
#' *Combined Community Metrics Dataframe*
#' 
#' Step 1: combine diversity + functional group summaries
  community_sum_combined <- bind_rows(structure_sum_sqrt, fungrp_sum) %>%
    mutate(metric = fct_recode(as.character(metric), "Forb" = "F") %>%
    factor(levels = c("Diversity", "Richness", "Evenness", "G_C4", "G_C3", "G_A-BI", "Forb")),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
    graze_trt = fct_recode(as.character(graze_trt),
            "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
    factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    year = as.numeric(as.character(year)),
    treat = as.factor(paste(graze_trt, rain_trt, sep = "-")))
  
#' Step 2: wide mean and SE dataframes
  community_sum_mean <- community_sum_combined %>%
    select(-se) %>%
    pivot_wider(names_from = metric, values_from = mean) %>%
    arrange(graze_trt, year)
  
  community_sum_se <- community_sum_combined %>%
    select(-mean) %>%
    pivot_wider(names_from = metric, values_from = se) %>%
    arrange(graze_trt, year)
  
#' 
#' Step 3: staggered x-axis positions (generalized to any years)
#' 
  year_seq <- community_sum_mean %>%
    distinct(year, graze_trt) %>%
    mutate(year_seq = case_when(
      graze_trt == "Ungrazed" ~ year - 0.15,
      graze_trt == "Fall Graze" ~ year,
      graze_trt == "Spring Graze" ~ year + 0.15 ))
  
  community_sum_mean <- community_sum_mean %>%
    left_join(year_seq, by = c("year", "graze_trt"))

#'
#' **ANALYSIS: Species composition**
#'   
  
#' PERMANOVA
spp_perm_all_sqrt <- adonis2(spp_rel_sqrt ~ rain_trt * graze_trt * year, 
                     strata = comp_wide_info$block, data = comp_wide_info, 
                     permutations = 999, method = "bray", by = "terms")
spp_perm_all_sqrt
  
#' NMDS
nmds_allyrs_sqrt <- metaMDS(spp_rel_sqrt, k = 3, trymax = 100, distance = "bray")
nmds_allyrs_sqrt
stressplot(nmds_allyrs_sqrt)
  
#' Plot scores 
plot_scores_all <- data.frame(nmds_allyrs_sqrt$points) %>%
  bind_cols(comp_wide_info) %>%
  mutate(b_p = paste(block, plot, sep = "_"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
    graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze"   = "FG", "Spring Graze" = "SG") %>%
    factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    year = factor(year, levels = c("2018","2019","2020","2021", "2022","2023","2024","2025"))) 

#' Mean relative abundance per species per year
spp_mean_rel_sqrt <- comp_dat_long %>%
  group_by(year, species) %>%
  summarise(sqrt_rel_abun = mean(rel_cov, na.rm = TRUE), .groups = "drop")

#' Species mean relative abundance across all years
spp_mean_rel_sqrt_allyrs <- spp_mean_rel_sqrt %>%
  group_by(species) %>%
  summarise(sqrt_rel_abun = mean(sqrt_rel_abun), .groups = "drop")

#' Species scores
nmds_species <- data.frame(nmds_allyrs_sqrt$species) %>%
  mutate(species = row.names(.)) %>%
  left_join(species_info, by = c("species" = "species_clean")) %>%  
  left_join(spp_mean_rel_sqrt_allyrs, by = "species")
  
#' Trajectory dataframe 
  plot_scores_trajectories <- plot_scores_all %>%
    filter(year %in% c(2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025)) %>%
    select(block, plot, graze_trt, rain_trt, year, MDS1, MDS2) %>%
    pivot_wider(names_from = year, values_from = c(MDS1, MDS2), names_sep = "_") 
  
#' Formatting helper - apply once to both dataframes
  format_scores <- function(df, exclude_year = NULL) {
    df %>%
      {if (!is.null(exclude_year)) filter(., year != exclude_year) else .} %>%
      mutate(
        graze_trt = fct_recode(as.character(graze_trt), 
                    "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
        factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
        rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
        year = factor(year, levels = c("2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025")))
  }
  
#' Apply formatting
#' *why are we excluding 2019?*
#plot_scores_no2019     <- format_scores(plot_scores_all, exclude_year = 2019) 
#plot_scores_2018_20    <- format_scores(plot_scores_trajectories)


#' *NMDS FIGURE - Figure S6*
#' 
#' Composition over 3 years, from a single NMDS (3 axes) with all years
#'    and vectors showing change over time 

# Extract PERMANOVA p-values dynamically
p_rain <- spp_perm_all_sqrt["rain_trt", "Pr(>F)"]
p_int <- spp_perm_all_sqrt["rain_trt:graze_trt", "Pr(>F)"]
p_year <- spp_perm_all_sqrt["year", "Pr(>F)"]
perm_label <- paste0("PerMANOVA:\n  Rainfall p=", p_rain,
                     "\n  Rainfall:Graze p=", p_int,
                     "\n  Year p=", p_year)

#' Figure 1: trajectory plot
#' Create one segment per consecutive year pair automatically
all_years <- c(2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025)

arrow_segments <- map_dfr(seq_along(all_years[-length(all_years)]), function(i) {
  yr1 <- all_years[i]
  yr2 <- all_years[i + 1]
  plot_scores_trajectories %>%
    mutate(rain_trt = factor(rain_trt, levels = c("D", "C", "W"))) %>%
    transmute(block, plot, graze_trt, rain_trt,
              MDS1_start = .data[[paste0("MDS1_", yr1)]],
              MDS2_start = .data[[paste0("MDS2_", yr1)]],
              MDS1_end   = .data[[paste0("MDS1_", yr2)]],
              MDS2_end   = .data[[paste0("MDS2_", yr2)]],
              year_from  = yr1, year_to = yr2)
}) %>%
  filter(!is.na(MDS1_start) & !is.na(MDS1_end) &   # remove missing segments
           !is.na(MDS2_start) & !is.na(MDS2_end))

# Geom_segment for all years
geom_segment(aes(x = MDS1_start, xend = MDS1_end,
                 y = MDS2_start, yend = MDS2_end,
                 color = rain_trt, lty = graze_trt),
             linewidth = 0.6,
             arrow = arrow(length = unit(0.2, "cm")),
             data = arrow_segments)


n_years <- length(all_years)
alpha_vals <- seq(0.15, 1, length.out = n_years)
names(alpha_vals) <- as.character(all_years)

fig_nmds_vec <- ggplot() +
  geom_point(aes(x = MDS1, y = MDS2, col = rain_trt, shape = graze_trt, alpha = as.factor(year)),
           size = 2.5, data = plot_scores_all) +
  # Consecutive year arrows - built dynamically
  geom_segment(aes(x = MDS1_start, xend = MDS1_end,
                   y = MDS2_start, yend = MDS2_end,
                   color = rain_trt, lty = graze_trt),
               linewidth = 0.6, arrow = arrow(length = unit(0.2, "cm")), data = arrow_segments) +
  annotate("text", x = -0.85, y = 0.85, label = "2018-2025", fontface = 2, size = 4, hjust = 0) +
  scale_color_manual(values = c("#d7191c", "#c994c7", "#2c7bb6")) +
  scale_alpha_manual(values = alpha_vals) +  # dynamic alpha for all years
  labs(x = "NMDS 1 - Species Composition",
       y = "NMDS 2 - Species Composition",
       col = "Rain trt.", shape = "Graze trt.",
       lty = "Graze trt.", alpha = "Year") +
  facet_grid(graze_trt~rain_trt) + 
  theme_bw()

#' view and save
fig_nmds_vec
ggsave(fig_nmds_vec, filename = paste0(wdpath, "Figures26/fig_nmds_vec.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")


# Figure 2: species scores
fig_nmds_species <- ggplot(data = nmds_species) +
  geom_point(aes(x = MDS1, y = MDS2, col = rain_trt, shape = graze_trt, alpha = as.factor(year)),  
             size = 2.5, data = plot_scores_all) +  
  geom_text(aes(x = MDS1, y = MDS2, size = sqrt_rel_abun, label = species)) +
  annotate("text", x = -1.5, y = 1.1, hjust = 0, fontface = 2, size = 3,
           label = "Species weighted by\nmean rel. abundance\n(sqrt-transformed)\n2018-2025") + 
  scale_color_manual(values = c("#d7191c", "#c994c7", "#2c7bb6")) +
  scale_alpha_manual(values = alpha_vals) + 
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "NMDS 1 - Species Composition", y = "NMDS 2 - Species Composition",
       col = "Rain trt.", shape = "Graze trt.") +
  theme_bw() +
  theme(legend.position = "none")

#' view and save
fig_nmds_species
ggsave(fig_nmds_species, filename = paste0(wdpath, "Figures26/fig_nmds_species.jpg"),
       device = "jpeg", height = 6.5, width = 5.5, units = "in")

#' **ANALYSIS: Linear Mixed Models**
#'  
#'  Compositional Responses:
#'  
#'  --Shannon Diversity
#'     (Supporting - Evenness, Richness)
#'  --C4 grasses
#'  --C3 grasses
#'  --Annual-Biennials
#'  --Forbs
#'  

#'  *Shannon Diversity*
#'  

#' Step 1: distribution check 
structure_all_sqrt %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt),
                           "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  
  ggplot() +
  geom_point(aes(x = rain_trt, y = shan, color = rain_trt),
             alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = shan, fill = rain_trt),
               alpha = 0.4, outlier.shape = 8, outlier.colour = "#7b3294") +
  scale_fill_manual(values = c("#d7191c", "#ffffbf", "#2c7bb6")) +
  scale_color_manual(values = c("#d7191c", "#ffffbf", "#2c7bb6")) +
  labs(y = "Shannon Diversity") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' run linear mixed models for shannon diversity 
#' using rain_trt as a category
div_mod_cat <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_all_sqrt)
  div_cat_means <- emmeans(div_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
  div_cat_contrasts <- div_cat_means %>% 
    contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")
  
#' using total precip cm 
div_mod_precip <- lmer(shan ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip)
summary(div_mod_precip)

#' Check residuals
plot(div_mod_cat)
qqnorm(residuals(div_mod_cat))

#' Model summary
anova(div_mod_cat)
r2(div_mod_cat)
summary(div_mod_cat)

#' Step 4: Shannon diversity figure - generalized for any years
all_years <- sort(unique(structure_sum_sqrt$year))

structure_sum_div <- structure_sum_sqrt %>%
  filter(metric == "Diversity") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
         rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
         year      = as.numeric(as.character(year)))

fig_div <- ggplot(structure_sum_div, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                linewidth = 1, width = 0, alpha = 0.45,
                position = position_dodge(0.3)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +         
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet"), name = "Rain Treatment") +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  labs(x = "Year", y = "Shannon Diversity", col = "Rain Treatment", shape = "Graze Treatment") +
  ggtitle("Shannon Diversity", subtitle = "Rain x Graze Treatments") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_div
ggsave(fig_div, filename = paste0(wdpath, "Figures26/fig_div.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' LMM predicted values from coninuous rainfall model
div_precip_pred <- ggpredict(div_mod_precip, terms = c("total_precip [all]", "graze_trt"))

pred_div_df <- as.data.frame(div_precip_pred) %>%
  mutate(graze_trt = fct_recode(as.character(group), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"))

slope_df <- slope_df %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt),"Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
         label = paste0("Slope: ", round(total_precip.trend, 3), "\np = ", round(as.data.frame(test(slopes))$p.value, 3)))

#' extract slope and p values for lmm
div_trends <- emtrends(div_mod_precip, ~ graze_trt, var = "total_precip")

stats_div_df <- as.data.frame(test(div_trends)) %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ "ns"),
    label = paste0("beta == ", round(total_precip.trend, 3), "~('", stars, "')"))

#' extract R2 values from model
model_r2 <- r2(div_mod_precip)
m_r2 <- round(model_r2$R2_marginal, 3)
c_r2 <- round(model_r2$R2_conditional, 3)

r2_label <- paste0("Marginal ~ R^2 == ", m_r2, " * ',' ~ Conditional ~ R^2 == ", c_r2)

#' plot LMM trend line and raw data 
fig_div_precip <- ggplot() +
  # Raw data points
  geom_point(data = structure_precip %>%  
               mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG")),
             aes(x = total_precip, y = shan, col = rain_trt, shape = graze_trt), alpha = 0.8, size = 2.5) +
  geom_ribbon(data = pred_div_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = graze_trt), alpha = 0.15) +
  geom_line(data = pred_div_df, aes(x = x, y = predicted, group = graze_trt), color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet"), name = "Rain Treatment") +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  scale_fill_manual(values = c("Ungrazed" = "darkgreen", "Fall Graze" = "orange", "Spring Graze" = "purple"), 
                    labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  labs(x = "Total Growing Season Precipitation (cm)", y = "Shannon Diversity", col = "Rain Treatment", shape = "Graze Treatment") +
  ggtitle("Shannon Diversity", subtitle = "Total Precipitation (cm)") +
  facet_wrap(~graze_trt) +
  geom_text(data = stats_div_df, aes(x = 60, y = 3.1, label = label), parse = TRUE, size = 4) +
  annotate("text", x = Inf, y = -Inf, label = r2_label, parse = TRUE, hjust = 1.1, vjust = -1.5, size = 2.5, fontface = "italic") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

#' view and save
fig_div_precip
ggsave(fig_div_precip, filename = paste0(wdpath, "Figures26/fig_div_precip.jpg"),
       device = "jpeg", height = 4.5, width = 7, units = "in")

#' Combined figure: Shannon Diversity organized two ways 
fig_combined_div <- (fig_div / fig_div_precip) +
  plot_annotation(title = "Shannon Diversity", subtitle = "By rain treatment and total precipitation (cm)", tag_levels = "A") + 
  plot_layout(guides = "collect", axis_titles = "collect")

#' view and save
fig_combined_div
ggsave(fig_combined_div, filename = paste0(wdpath, "Figures26/fig_combined_div.jpg"),
       device = "jpeg", height = 10, width = 10, units = "in")

#'  
#'  * Richness *
#' 

#' Distribution check 
structure_all_sqrt %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
         rain_trt = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = rich, color = rain_trt),
             alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = rich, fill = rain_trt),
               alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "Richness") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' run linear mixed models for richness 
#' using rain_trt as a category
rich_mod_cat <- lmer(rich ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_all_sqrt)
  rich_cat_means <- emmeans(rich_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
  rich_cat_contrasts <- rich_cat_means %>% 
  contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")
  
#' using total precip cm 
rich_mod_precip <- lmer(rich ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip) 
summary(rich_mod_precip)

#' Check residuals
plot(rich_mod_cat)
qqnorm(residuals(rich_mod_cat))

#' Model summaries
anova(rich_mod_cat)
r2(rich_mod_cat)
summary(rich_mod_cat)

#' Richness Figure - generalized for all years
structure_sum_rich <- structure_sum_sqrt %>%
  filter(metric == "Richness") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
    year = as.numeric(as.character(year)))

fig_rich <- ggplot(structure_sum_rich, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                linewidth = 1, width = 0, alpha = 0.45,         
                position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Species Richness", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_rich
ggsave(fig_rich, filename = paste0(wdpath, "Figures26/fig_rich.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' Richness and total precip figure
#' LMM predicted values from coninuous rainfall model
rich_precip_pred <- ggpredict(rich_mod_precip, terms = c("total_precip [all]", "graze_trt"))

pred_rich_df <- as.data.frame(rich_precip_pred) %>%
  mutate(graze_trt = fct_recode(as.character(group), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"))

stats_rich_df <- emtrends(rich_mod_precip, ~ graze_trt, var = "total_precip") %>%
  test() %>%           
  as.data.frame() %>%  
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
         stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05  ~ "*", TRUE ~ "ns"),
         label = paste0("beta == ", round(total_precip.trend, 3), "~('", stars, "')"))

#' extract R2 values from model
model_r2 <- r2(rich_mod_precip)
m_r2 <- round(model_r2$R2_marginal, 3)
c_r2 <- round(model_r2$R2_conditional, 3)

r2_label <- paste0("Marginal ~ R^2 == ", m_r2, " * ',' ~ Conditional ~ R^2 == ", c_r2)

#' plot LMM trend line and raw data 
fig_rich_precip <- ggplot() +
  # Raw data points
  geom_point(data = structure_precip %>% 
               mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG")),
             aes(x = total_precip, y = rich, col = rain_trt, shape = graze_trt), alpha = 0.8, size = 2.5) +
  geom_ribbon(data = pred_rich_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = graze_trt), alpha = 0.15) +
  geom_line(data = pred_rich_df, aes(x = x, y = predicted, group = graze_trt), color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet"), name = "Rain Treatment") +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  scale_fill_manual(values = c("Ungrazed" = "darkgreen", "Fall Graze" = "orange", "Spring Graze" = "purple"), 
                    labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  labs(x = "Total Growing Season Precipitation (cm)", y = "Richness", col = "Rain Treatment", shape = "Graze Treatment") +
  ggtitle("Richness", subtitle = "Total Precipitation (cm)") +
  #coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~graze_trt) +
  geom_text(data = stats_even_df, aes(x = 60, y = 30, label = label), parse = TRUE, size = 4) +
  annotate("text", x = Inf, y = -Inf, label = r2_label, parse = TRUE, hjust = 1.1, vjust = -1.5, size = 2.5, fontface = "italic") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

#' view and save
fig_rich_precip
ggsave(fig_rich_precip, filename = paste0(wdpath, "Figures26/fig_rich_precip.jpg"),
       device = "jpeg", height = 4.5, width = 7, units = "in")

#' Combined figure: Richness organized two ways 
fig_combined_rich <- (fig_rich / fig_rich_precip) +
  plot_annotation(title = "Richness", subtitle = "By rain treatment and total precipitation (cm)", tag_levels = "A") + 
  plot_layout(guides = "collect", axis_titles = "collect")

#' view and save
fig_combined_rich
ggsave(fig_combined_rich, filename = paste0(wdpath, "Figures26/fig_combined_rich.jpg"),
       device = "jpeg", height = 10, width = 10, units = "in")

#'  
#'  * Evenness *
#'  

#' Distribution check 
structure_raw_sqrt %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = even, color = rain_trt), alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = even, fill = rain_trt), alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "Evenness") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' run linear mixed models for evenness 
#' using rain_trt as a category
even_mod_cat <- lmer(even ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_all_sqrt)
  even_cat_means <- emmeans(even_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
  even_cat_contrasts <- even_cat_means %>% 
    contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")

#' using total precip cm
even_mod_precip <- lmer(even ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip) 
summary(even_mod_precip)

#' Check residuals
plot(even_mod_cat)
qqnorm(residuals(even_mod_cat))

#' Model summaries
anova(even_mod_cat)
r2(even_mod_cat)
summary(even_mod_cat)

#' Evenness Figure
structure_sum_even <- structure_sum_sqrt %>%
  filter(metric == "Evenness") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt = factor(rain_trt,  levels = c("D", "C", "W")), year = as.numeric(as.character(year)))

fig_even <- ggplot(structure_sum_even, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  # Raw data points
  geom_point(data = structure_raw_sqrt %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
      factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
      rain_trt = factor(rain_trt, levels = c("D", "C", "W")), year = as.numeric(as.character(year))),
             aes(x = year, y = even, col = rain_trt, shape = graze_trt),
             alpha = 0.15, size = 1, position = position_dodge(0.3), show.legend = FALSE) +
  # Means
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.3)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Species Evenness", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))
  #coord_cartesian(ylim = c(0, 0.75))
  
#' view and save
fig_even
ggsave(fig_even, filename = paste0(wdpath, "Figures26/fig_even.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' Evenness and total precip figure
#' LMM predicted values from coninuous rainfall model
even_precip_pred <- ggpredict(even_mod_precip, terms = c("total_precip [all]", "graze_trt"))

pred_even_df <- as.data.frame(even_precip_pred) %>%
  mutate(graze_trt = fct_recode(as.character(group), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"))

stats_even_df <- emtrends(even_mod_precip, ~ graze_trt, var = "total_precip") %>%
  test() %>%           
  as.data.frame() %>%  
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
         stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05  ~ "*", TRUE ~ "ns"),
         label = paste0("beta == ", round(total_precip.trend, 3), "~('", stars, "')"))

#' extract R2 values from model
model_r2 <- r2(even_mod_precip)
m_r2 <- round(model_r2$R2_marginal, 3)
c_r2 <- round(model_r2$R2_conditional, 3)

r2_label <- paste0("Marginal ~ R^2 == ", m_r2, " * ',' ~ Conditional ~ R^2 == ", c_r2)

#' plot LMM trend line and raw data 
fig_even_precip <- ggplot() +
  # Raw data points
  geom_point(data = structure_precip %>% 
               mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG")),
             aes(x = total_precip, y = even, col = rain_trt, shape = graze_trt), alpha = 0.8, size = 2.5) +
  geom_ribbon(data = pred_even_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = graze_trt), alpha = 0.15) +
  geom_line(data = pred_even_df, aes(x = x, y = predicted, group = graze_trt), color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet"), name = "Rain Treatment") +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  scale_fill_manual(values = c("Ungrazed" = "darkgreen", "Fall Graze" = "orange", "Spring Graze" = "purple"), 
                    labels = c("Ungrazed", "Fall Graze", "Spring Graze"), name = "Graze Treatment") +
  labs(x = "Total Growing Season Precipitation (cm)", y = "Evenness", col = "Rain Treatment", shape = "Graze Treatment") +
  ggtitle("Evenness", subtitle = "Total Precipitation (cm)") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~graze_trt) +
  geom_text(data = stats_even_df, aes(x = 60, y = 0.95, label = label), parse = TRUE, size = 4) +
  annotate("text", x = Inf, y = -Inf, label = r2_label, parse = TRUE, hjust = 1.1, vjust = -1.5, size = 2.5, fontface = "italic") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

#' view and save
fig_even_precip
ggsave(fig_even_precip, filename = paste0(wdpath, "Figures26/fig_even_precip.jpg"),
       device = "jpeg", height = 4.5, width = 7, units = "in")

#' Combined figure: Evenness organized two ways 
fig_combined_even <- (fig_even / fig_even_precip) +
  plot_annotation(title = "Evenness", subtitle = "By rain treatment and total precipitation (cm)", tag_levels = "A") + 
  plot_layout(guides = "collect", axis_titles = "collect")

#' view and save
fig_combined_even
ggsave(fig_combined_even, filename = paste0(wdpath, "Figures26/fig_combined_even.jpg"),
       device = "jpeg", height = 10, width = 10, units = "in")

#' 
#' *Combined figure for diversity, richness, and evenness*

#' Combined figure: diversity richness, evenness, 
fig_combined_dre <- (fig_div / fig_rich / fig_even) +
  plot_annotation(title = "Community Structure Metrics", tag_levels = "A") + 
  plot_layout(guides = "collect", axis_titles = "collect")

#' view and save
fig_combined_dre
ggsave(fig_combined_dre, filename = paste0(wdpath, "Figures26/fig_combined_dre.jpg"),
       device = "jpeg", height = 10, width = 10, units = "in")


#'  
#'  * C4 grasses *
#'

#' Distribution check - inline, no separate object
scale_fungrp(fungrp_con, "G_C4") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = fungrp_cov, color = rain_trt), alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = fungrp_cov, fill = rain_trt),
               alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "C4 Perennial Grass") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' Linear mixed models
#' using rain_trt as a category
C4_mod_cat <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "G_C4"))
  C4_cat_means <- emmeans(C4_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
  C4_cat_contrasts <- C4_cat_means %>% 
  contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")

#' using total precip cm 
C4_mod_precip <- lmer(fungrp_cov ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip) 

#' Check residuals
plot(C4_mod_cat)
qqnorm(residuals(C4_mod_cat))

#' Model summaries
anova(C4_mod_cat)
r2(C4_mod_cat)
summary(C4_mod_cat)

#' C4 Grass Figure - generalized for all years
c4_sum <- fungrp_sum %>%
  filter(metric == "G_C4") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
    year = as.numeric(as.character(year)))

fig_c4 <- ggplot(c4_sum, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  # Raw data points behind means
  geom_point(data = scale_fungrp(fungrp_con, "G_C4") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
  factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")), rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
  year = as.numeric(as.character(year))), aes(x = year, y = fungrp_cov, col = rain_trt, shape = graze_trt),
             alpha = 0.15, size = 1, position = position_dodge(0.35)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "C4 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  coord_cartesian(ylim = c(-0, 1)) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_c4
ggsave(fig_c4, filename = paste0(wdpath, "Figures26/fig_c4.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' 
#' *C3 perennial grasses*
#'

#' Distribution check - inline
scale_fungrp(fungrp_con, "G_C3") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = fungrp_cov, color = rain_trt), alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = fungrp_cov, fill = rain_trt), alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "C3 Perennial Grass") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' Linear mixed models
#' using rain_trt as a category
C3_mod_cat <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "G_C4"))
C3_cat_means <- emmeans(C3_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
C3_cat_contrasts <- C3_cat_means %>% 
  contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")

#' using total precip cm 
C3_mod_precip <- lmer(fungrp_cov ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip) 

#' Check residuals
plot(C3_mod_cat)
qqnorm(residuals(C3_mod_cat))

#' Model summaries
anova(C3_mod_cat)
r2(C3_mod_cat)
summary(C3_mod_cat)

# Summary dataframe for figure
c3_sum <- fungrp_sum %>%
  filter(metric == "G_C3") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt = factor(rain_trt, levels = c("D", "C", "W")), year = as.numeric(as.character(year)))

# Figure
fig_c3 <- ggplot(c3_sum, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(data = scale_fungrp(fungrp_con, "G_C3") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
  factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")), rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
  year = as.numeric(as.character(year))), aes(x = year, y = fungrp_cov, col = rain_trt, shape = graze_trt),
    alpha = 0.15, size = 1, position = position_dodge(0.35)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "C3 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_c3
ggsave(fig_c3, filename = paste0(wdpath, "Figures26/fig_c3.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' 
#' *Annual-Biennials*
#' 

#' Distribution check - inline
scale_fungrp(fungrp_con, "G_A-BI") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = fungrp_cov, color = rain_trt), alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = fungrp_cov, fill = rain_trt), alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "Annual-Biennials") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' Outlier check - compare full vs outlier-removed model
ann_mod_full <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|b_p), data = scale_fungrp(fungrp_con, "G_A-BI"))
ann_mod_out  <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|b_p), data = scale_fungrp(fungrp_con, "G_A-BI") %>% filter(b_p != "3_4"))

#' Compare full vs outlier-removed
anova(ann_mod_full)
anova(ann_mod_out)
plot(ann_mod_full);  plot(ann_mod_out)       # fitted vs residuals
qqnorm(residuals(ann_mod_full)); qqnorm(residuals(ann_mod_out))  # QQ plots

#' Conclusion: remove outlier 3_4 - improves diagnostics, minimally affects results

# Helper to remove outlier from any fungrp contrast set
remove_outlier <- function(df, grp, outlier = "3_4") {
  scale_fungrp(df, grp) %>% filter(b_p != outlier)
}

# Models with outlier removed - note: block removed due to singularity
ann_mod_con <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|b_p), data = remove_outlier(fungrp_con, "G_A-BI"))
ann_mod_dry <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|b_p), data = remove_outlier(fungrp_dry, "G_A-BI"))
ann_mod_wet <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|b_p), data = remove_outlier(fungrp_wet, "G_A-BI"))

#' Check residuals
plot(ann_mod_con)
qqnorm(residuals(ann_mod_con))

#' Model summaries
Anova(ann_mod_con)
anova(ann_mod_con)
r2(ann_mod_con)
summary(ann_mod_con)

#' Extract fixed effects
ann_fixef_all <- bind_rows(extract_fixef(ann_mod_con, "ann"),
  extract_fixef(ann_mod_dry, "ann"),
  extract_fixef(ann_mod_wet, "ann")) %>%
  distinct() %>%
  arrange(term)
ann_fixef_all

#' Summary dataframe - remove outlier before summarising
ann_sum <- fungrp_sum %>%
  filter(metric == "G_A-BI") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
    year = as.numeric(as.character(year)))

#' Note: fungrp_sum should be recalculated without outlier for figure accuracy
ann_sum_out <- comp_dat_long_fungrp %>%
  filter(fun_grp == "G_A-BI", b_p != "3_4") %>%  # remove outlier
  group_by(year, rain_trt, graze_trt, fun_grp) %>%
  summarise(mean = mean(rel_cov, na.rm = TRUE), se = se(rel_cov), .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year)), rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
    graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
      factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")))

#' Annual - biannual Figure
fig_ann <- ggplot(ann_sum_out, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(data = remove_outlier(fungrp_con, "G_A-BI") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
    factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
    year = as.numeric(as.character(year))), aes(x = year, y = fungrp_cov, col = rain_trt, shape = graze_trt),
      alpha = 0.15, size = 1, position = position_dodge(0.35)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Annual-Biennials\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_ann
ggsave(fig_ann, filename = paste0(wdpath, "Figures26/fig_ann.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")

#' 
#' *Forbs*
#' 

#' Distribution check
scale_fungrp(fungrp_con, "F") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG"),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W"))) %>%
  ggplot() +
  geom_point(aes(x = rain_trt, y = fungrp_cov, color = rain_trt), alpha = 0.4, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(x = rain_trt, y = fungrp_cov, fill = rain_trt), alpha = 0.4, outlier.shape = 8, outlier.colour = "darkred") +
  scale_fill_manual(values  = c("goldenrod2", "forestgreen", "navyblue")) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "navyblue")) +
  labs(y = "Perennial Forbs") +
  facet_grid(year ~ graze_trt, scales = "free") +
  theme_bw()

#' Models
f_mod_con <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "F"))
f_mod_dry <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_dry, "F"))
f_mod_wet <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_wet, "F"))

#' Linear mixed models
#' using rain_trt as a category
f_mod_cat <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "G_C4"))
f_cat_means <- emmeans(f_mod_cat, specs = c("rain_trt", "graze_trt", "year")) 
f_cat_contrasts <- f_cat_means %>% 
  contrast(method = 'pairwise', by = c("graze_trt", "year"), adjust = "none")

#' using total precip cm 
f_mod_precip <- lmer(fungrp_cov ~ total_precip * graze_trt * year + (1|block) + (1|b_p), data = structure_precip) 

#' Check residuals
plot(f_mod_cat)
qqnorm(residuals(f_mod_cat))

#' Model summaries
Anova(f_mod_cat)
anova(f_mod_cat)
r2(f_mod_cat)
summary(f_mod_cat)

# Summary dataframe
f_sum <- fungrp_sum %>%
  filter(metric == "F") %>%
  mutate(graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W")), year = as.numeric(as.character(year)))

#' Forbs Figure
fig_forb <- ggplot(f_sum, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(data = scale_fungrp(fungrp_con, "F") %>%
  mutate(graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
  factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
     rain_trt  = factor(rain_trt, levels = c("D", "C", "W")), year = as.numeric(as.character(year))),
     aes(x = year, y = fungrp_cov, col = rain_trt, shape = graze_trt), alpha = 0.15, size = 1, position = position_dodge(0.35)) +
     geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Perennial Forbs\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_forb
ggsave(fig_forb, filename = paste0(wdpath, "Figures26/fig_forb.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")


#' 
#' *Combined figure for C4, C3, annuals, forbs*

#' Combined figure: diversity richness, evenness, 
fig_combined_43af <- (fig_c4 + fig_c3) / (fig_ann + fig_forb) +
  plot_annotation(title = "Community Composition", tag_levels = "A") + 
  plot_layout(guides = "collect", axis_titles = "collect")

#' view and save
fig_combined_43af
ggsave(fig_combined_43af, filename = paste0(wdpath, "Figures26/fig_combined_43af.jpg"),
       device = "jpeg", height = 10, width = 10, units = "in")

#' 
#' ** Supporting data figure: Absolute relative covers for each functional groups**
#' 
fg_abs_full_allyears <- fungrp_df %>%
  pivot_wider(names_from = fun_grp, 
              values_from = fungrp_cov,
              values_fill = 0) %>%
  rename_with(~ paste0(tolower(gsub("G_", "", .)), "_abs"), 
              .cols = -c(year, b_p, block, rain_trt, graze_trt)) %>%
  #removes outlier
  mutate(across(contains("ann"), ~ if_else(b_p == "3_4", NA_real_, .)))

fg_abs_sum_dat <- fg_abs_full_allyears %>%
  #mutate(ann_abs = if_else(b_p == "3_4", NA_real_, ann_abs)) %>%
  pivot_longer(cols = c4_abs:f_abs, names_to = "fg", values_to = "abs_cov") %>%
  group_by(fg, year, graze_trt, rain_trt) %>%
  summarise(abs_cov_mean = mean(abs_cov, na.rm = TRUE),
            abs_cov_se   = se(abs_cov, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year)), rain_trt = factor(rain_trt, levels = c("D", "C", "W")),
    graze_trt = fct_recode(as.character(graze_trt), "Ungrazed" = "C", "Fall Graze" = "FG", "Spring Graze" = "SG") %>%
    factor(levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    fg = fct_recode(fg, "Annual-Biennials" = "ann_abs", "C3 Perennial Grasses" = "c3_abs",
                    "C4 Perennial Grasses" = "c4_abs", "Perennial Forbs" = "f_abs") %>%
    factor(levels = c("Annual-Biennials", "Perennial Forbs", "C3 Perennial Grasses", "C4 Perennial Grasses")))

fig_fg_abs <- ggplot(fg_abs_sum_dat, 
                     aes(x = year, y = abs_cov_mean, 
                         col = rain_trt, shape = graze_trt)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = abs_cov_mean - abs_cov_se, 
                    ymax = abs_cov_mean + abs_cov_se),
                linewidth = 1, width = 0, alpha = 0.45,    # linewidth not size
                position = position_dodge(0.35)) +
  scale_x_continuous(breaks = all_years,                   # dynamic
                     limits = c(min(all_years) - 0.25, 
                                max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6")) +
  scale_fill_manual(values  = c("#d7191c", "tan", "#2c7bb6")) +
  labs(x = "Year", y = "Absolute Cover (%)", 
       col = "Rain trt.", shape = "Graze trt.") +
  facet_grid(fg ~ graze_trt, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_fg_abs

#' 
#' *species richness X precip figure*
precipRich <- precip_annual %>%
  #left_join(treatments, by = c("block", "plot")) %>% 
  left_join(structure_sum_rich, by = c("year", "rain_trt"))

fig_precipRich <- ggplot(precipRich, aes(x = as.factor(year))) +
  geom_bar(aes(y = total_precip, fill = rain_trt, group = rain_trt),
           stat = "identity", width = 0.35, position = position_dodge(0.4), alpha = 0.3) +
  geom_point(aes(y = mean * max(total_precip, na.rm = TRUE) / 30, 
                 col = rain_trt, shape = graze_trt), size = 2.5, alpha = 0.8, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = (mean - se) * max(total_precip, na.rm = TRUE) / 30,  
                    ymax = (mean + se) * max(total_precip, na.rm = TRUE) / 30,
                    group = rain_trt), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.4)) +  
  scale_y_continuous(name = "Total Precipitation (cm)", limits = c(0, NA),  # force y to start at 0
    sec.axis = sec_axis(~ . / max(precipRich$total_precip, na.rm = TRUE) * 30, name = "Species Richness")) +
  scale_fill_manual(values  = c("#d7191c", "tan", "#2c7bb6"), labels  = c("Drought", "Control", "Wet")) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", fill = "Rain Treatment", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_precipRich

precipRich2 <- precipRich %>% 
  group_by(rain_trt, graze_trt)

fig_precipRich2 <- ggplot(precipRich2, aes(x = total_precip, y = mean, col = rain_trt, shape = graze_trt)) + 
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(limits = c(0, max(precipRich2$total_precip) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Precipitation (cm)", y = "Species Richness", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10))

fig_precipRich2

fig_combined_precipRich <- fig_precipRich / fig_precipRich2 +
  plot_annotation(title = "Precipitation and Species Richness", tag_levels = "A") + 
  plot_layout(guides = "collect")

fig_combined_precipRich

#' 
#' *Species Diversity X precip figure*
precipDiv <- precip_annual %>%
  #left_join(treatments, by = c("block", "plot")) %>% 
  left_join(structure_sum_div, by = c("year", "rain_trt"))

fig_precipDiv <- ggplot(precipDiv, aes(x = as.factor(year))) +
  geom_bar(aes(y = total_precip, fill = rain_trt, group = rain_trt),
           stat = "identity", width = 0.35, position = position_dodge(0.4), alpha = 0.3) +
  #geom_point(aes(y = mean * max(total_precip, na.rm = TRUE) / 3, 
  #col = rain_trt, shape = graze_trt), size = 2.5, alpha = 0.8, position = position_dodge(0.4)) +
  geom_point(aes(y = mean, col = rain_trt, shape = graze_trt), size = 2.5, alpha = 0.8, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.4)) +
  #geom_errorbar(aes(ymin = (mean - se) * max(total_precip, na.rm = TRUE) / 3,  
                    #ymax = (mean + se) * max(total_precip, na.rm = TRUE) / 3,
                    #group = rain_trt), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.4)) +  
  scale_y_continuous(name = "Shannon Diversity", 
                     sec.axis = sec_axis(limits = c(0, ~ . /max(precipDiv$mean, na.rm = T)) * 100, name = "Total Precipitation (cm)")) +
  #scale_y_continuous(name = "Total Precipitation (cm)", limits = c(0, NA),  
                     #sec.axis = sec_axis(~ . / max(precipDiv$total_precip, na.rm = TRUE) * 3, name = "Shannon Diversity")) +
  scale_fill_manual(values  = c("#d7191c", "tan", "#2c7bb6"), labels  = c("Drought", "Control", "Wet")) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", fill = "Rain Treatment", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

# 1. Calculate a scaling factor
# This ensures the max precip aligns with the top of your diversity scale
scaleFactor_precipDiv <- max(precipDiv$total_precip, na.rm = TRUE) / max(precipDiv$mean, na.rm = TRUE)

fig_precipDiv <- ggplot(precipDiv, aes(x = as.factor(year))) +
  geom_bar(aes(y = total_precip / scaleFactor_precipDiv, fill = rain_trt),
           stat = "identity", width = 0.35, position = position_dodge(0.4), alpha = 0.3) +
  geom_point(aes(y = mean, col = rain_trt, shape = graze_trt), 
             size = 2.5, alpha = 0.8, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = rain_trt), 
                linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.4)) +
    scale_y_continuous(name = "Shannon Diversity", limits = c(0, NA),  
    sec.axis = sec_axis(~ . * scaleFactor_precipDiv, name = "Total Precipitation (cm)")) +
  scale_fill_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", fill = "Rain Treatment", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_precipDiv

precipDiv2 <- precipDiv %>% 
  group_by(rain_trt, graze_trt)

fig_precipDiv2 <- ggplot(precipDiv2, aes(x = total_precip, y = mean, col = rain_trt, shape = graze_trt)) + 
  geom_point(aes(alpha = year), size = 2.5, position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(limits = c(0, max(precipDiv2$total_precip) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Precipitation (cm)", y = "Shannon Diversity", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10))

fig_precipDiv2

fig_combined_precipDiv <- fig_precipDiv / fig_precipDiv2 +
  plot_annotation(title = "Precipitation and Shannon Diversity", tag_levels = "A") + 
  plot_layout(guides = "collect")

fig_combined_precipDiv

#' 
#' *precipitation X c3 cover Figure*
precipC3 <- precip_annual %>%
  #left_join(treatments, by = c("block", "plot")) %>% 
  left_join(c3_sum, by = c("year", "rain_trt"))

fig_precipC3 <- ggplot(precipC3, aes(x = total_precip, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_text(aes(label = year), size = 3, vjust = -1)  +
  geom_point(aes(size = year), position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(limits = c(0, max(precipDiv2$total_precip) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Total Precipitation (cm)", y = "C3 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  #coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_precipC3

#' 
#' *precipitation X c4 cover Figure*
precipC4 <- precip_annual %>%
  #left_join(treatments, by = c("block", "plot")) %>% 
  left_join(c4_sum, by = c("year", "rain_trt"))

fig_precipC4 <- ggplot(precipC4, aes(x = total_precip, y = mean, col = rain_trt, shape = graze_trt, group = interaction(rain_trt, graze_trt))) +
  geom_text(aes(label = year), size = 3, vjust = -1)  +
  geom_point(aes(size = year), position = position_dodge(0.35)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.35)) +
  scale_x_continuous(limits = c(0, max(precipDiv2$total_precip) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "tan", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Total Precipitation (cm)", y = "C4 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  #coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_precipC4

fig_combined_precipC3C4 <- fig_precipC3 / fig_precipC4 +
  plot_annotation(title = "Precipitation and Relative Proportions of C3 and C4 Perennial Grasses", tag_levels = "A") + 
  plot_layout(guides = "collect")

fig_combined_precipC3C4

#####
#'    
#'    
#'    
#'    [Part 2] : *FORAGE SERVICES RESPONSES*
#'                   to rainfall and grazing interactions
#'  
#'  
#'  [A] : Compile raw dataframes
#'  
#'  [B] : Prep model and figure dataframes
#'  
#'  [C] : Explore rainfall and grazing interactions on forage response
#'  
#'  [D] : Identify key indices of plant-driven functions and services
#'  
#'  
#'  
#'###         

##
#' 
#' 
#' [A]: **RAW DATA PREPARATION**
#' 
#'                Compile dataframes to assess rainfall and grazing effects on
#'                     plant communities and services
#' 
#'
#'##

#'
#'  **PLANT BIOMASS and BARE GROUND - Prep data**
#'     

#' 
#' *Data preparation*
#' 

#' Load dataframe with plant biomass (production), bare ground, and counts for some annual species
plant <- read.csv(paste0(wdpath, "Plant_Plot_Metrics_2018-23.csv"))

#' View structure of data and remove sample_ID column
str(plant)
