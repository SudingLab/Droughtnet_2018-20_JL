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
  filter(!species_clean %in% exclude_species)

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

#' Revised long matrix of total species cover
comp_long <- comp_wide %>%
  pivot_longer(cols = -c(year, block, plot), names_to = "species clean", values_to = "cover")

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
#' *Estimate composition metrics*
#' *1. Richness, Evenness, & Shannon Diversity* 
#' 

#' Step 1: calculate diversity metrics
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
#' Function for creating contrasts
set_contrasts <- function(df, year_ref, rain_ref, graze_ref) {
  df %>% mutate(year = fct_relevel(year, year_ref),
                rain_trt  = fct_relevel(rain_trt, rain_ref),
                graze_trt = fct_relevel(graze_trt, graze_ref))
  }

#' Control rain as reference ("does drought/wet differ from ambient?")
  structure_con <- set_contrasts(structure_all_sqrt, "2018", "C", "C")
#' Drought as reference ("does control/wet differ from drought?")
  structure_dry <- set_contrasts(structure_all_sqrt, "2018", "D", "C")
#' Wet as reference ("does control/drought differ from wet?")
  structure_wet <- set_contrasts(structure_all_sqrt, "2018", "W", "C")
  
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
  
#' Step 6: filter and scale by functional group using a function instead of 16 objects
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
  
#' If you need different contrasts for a specific model, just filter inline:
  scale_fungrp(fungrp_dry, "G_C4")  # drought reference for C4 grasses
  scale_fungrp(fungrp_wet, "G_C3")  # wet reference for C3 grasses 
  
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
  
#' ###########
#' Step 4: filter inline when needed
  # e.g. for C3 grasses:
  #community_sum_mean %>% 
    #filter(metric == "G_C3") %>%
    #left_join(community_sum_se %>% select(year, rain_trt, graze_trt, se = G_C3),
     #         by = c("year", "rain_trt", "graze_trt"))
  # OR keep combined and use facet_wrap(~metric) in ggplot 

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
})

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

#' Step 2: run models using contrast sets from set_contrasts() function
div_mod_con <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_con)
div_mod_dry <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_dry)
div_mod_wet <- lmer(shan ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_wet)

#' Check residuals
plot(div_mod_con)
qqnorm(residuals(div_mod_con))

#' Model summary
anova(div_mod_con)
r2(div_mod_con)
summary(div_mod_con)

#' Step 3: extract fixed effects 
extract_fixef <- function(model, response_name) {
  fix <- data.frame(round(fixef(model), 4)) %>%
    rename(fixef = 1) %>%
    rownames_to_column("term")
  ci <- data.frame(round(confint(model), 4)) %>%
    rownames_to_column("term")
  left_join(fix, ci, by = "term") %>%
    mutate(response = response_name)
}

#' Apply to each contrast set and combine
div_fixef_all <- bind_rows(
  extract_fixef(div_mod_con, "div"),
  extract_fixef(div_mod_dry, "div"),
  extract_fixef(div_mod_wet, "div")
) %>%
  distinct() %>%
  arrange(term)

#' Step 4: Shannon diversity figure - generalized for any years
all_years <- sort(unique(structure_sum_sqrt$year))

structure_sum_div <- structure_sum_sqrt %>%
  filter(metric == "Diversity") %>%
  mutate(
    graze_trt = factor(graze_trt, levels = c("Ungrazed", "Fall Graze", "Spring Graze")),
    rain_trt  = factor(rain_trt, levels = c("D", "C", "W")),
    year      = as.numeric(as.character(year))
  )

fig_div <- ggplot(structure_sum_div, aes(x = year, y = mean, col = rain_trt, shape = graze_trt)) +
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                linewidth = 1, width = 0, alpha = 0.45,
                position = position_dodge(0.3)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +         
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6")) +
  scale_fill_manual(values = c("#d7191c", "#8856a7", "#2c7bb6")) +
  labs(x = "Year", y = "Shannon Diversity",
       col = "Rain trt.", shape = "Graze trt.") +
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

#' Model - using named contrast sets
rich_mod_con <- lmer(rich ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_con)
rich_mod_dry <- lmer(rich ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_dry)
rich_mod_wet <- lmer(rich ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_wet)

#' Check residuals
plot(rich_mod_con)
qqnorm(residuals(rich_mod_con))

#' Model summaries
anova(rich_mod_con)
r2(rich_mod_con)
summary(rich_mod_con)

#' Extract fixed effects
rich_fixef_all <- bind_rows(
  extract_fixef(rich_mod_con, "rich"),
  extract_fixef(rich_mod_dry, "rich"),
  extract_fixef(rich_mod_wet, "rich")) %>%
  distinct() %>%
  arrange(term)

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
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
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

#'  
#'  * Evenness *
#'  

#' Distribution check - inline
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

#' Models
even_mod_con <- lmer(even ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_con)
even_mod_dry <- lmer(even ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_dry)
even_mod_wet <- lmer(even ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = structure_wet)

#' Check residuals
plot(even_mod_con)
qqnorm(residuals(even_mod_con))

#' Model summaries
anova(even_mod_con)
r2(even_mod_con)
summary(even_mod_con)

#' Fixed effects
even_fixef_all <- bind_rows(extract_fixef(even_mod_con, "even"),
  extract_fixef(even_mod_dry, "even"),
  extract_fixef(even_mod_wet, "even")) %>%
  distinct() %>%
  arrange(term)

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
             alpha = 0.15, size = 1, position = position_dodge(0.3)) +
  # Means
  geom_point(size = 2.5, alpha = 0.8, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), linewidth = 1, width = 0, alpha = 0.45, position = position_dodge(0.3)) +
  scale_x_continuous(breaks = all_years, limits = c(min(all_years) - 0.25, max(all_years) + 0.25)) +
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
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

#' Models - using named contrast sets
c4_mod_con <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "G_C4"))
c4_mod_dry <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_dry, "G_C4"))
c4_mod_wet <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_wet, "G_C4"))

#' Check residuals
plot(c4_mod_con)
qqnorm(residuals(c4_mod_con))

#' Model summaries
anova(c4_mod_con)
r2(c4_mod_con)
summary(c4_mod_con)

#' Extract fixed effects using reusable function
c4_fixef_all <- bind_rows(extract_fixef(c4_mod_con, "c4"),
  extract_fixef(c4_mod_dry, "c4"),
  extract_fixef(c4_mod_wet, "c4")) %>%
  distinct() %>%
  arrange(term)
c4_fixef_all

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
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "C4 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "none",
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

#' Models
c3_mod_con <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_con, "G_C3"))
c3_mod_dry <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_dry, "G_C3"))
c3_mod_wet <- lmer(fungrp_cov ~ rain_trt * graze_trt * year + (1|block) + (1|b_p), data = scale_fungrp(fungrp_wet, "G_C3"))

#' Check residuals
plot(c3_mod_con)
qqnorm(residuals(c3_mod_con))

#' Model summaries
anova(c3_mod_con)
r2(c3_mod_con)
summary(c3_mod_con)

#' Extract fixed effects
c3_fixef_all <- bind_rows(extract_fixef(c3_mod_con, "c3"),
  extract_fixef(c3_mod_dry, "c3"),extract_fixef(c3_mod_wet, "c3")) %>%
  distinct() %>%
  arrange(term)
c3_fixef_all

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
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "C3 Perennial Grasses\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
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
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Annual-Biennials\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
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

#' Check residuals
plot(f_mod_con)
qqnorm(residuals(f_mod_con))

#' Model summaries
Anova(f_mod_con)
anova(f_mod_con)
r2(f_mod_con)
summary(f_mod_con)

#' Extract fixed effects
f_fixef_all <- bind_rows(extract_fixef(f_mod_con, "f"),
  extract_fixef(f_mod_dry, "f"), extract_fixef(f_mod_wet, "f")) %>%
  distinct() %>%
  arrange(term)
f_fixef_all

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
  scale_color_manual(values = c("#d7191c", "#8856a7", "#2c7bb6"), labels = c("Drought", "Control", "Wet")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Ungrazed", "Fall Graze", "Spring Graze")) +
  labs(x = "Year", y = "Perennial Forbs\n(rel. cover proportion)", col = "Rain Treatment", shape = "Graze Treatment") +
  facet_wrap(~graze_trt) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1))

#' view and save
fig_forb
ggsave(fig_forb, filename = paste0(wdpath, "Figures26/fig_forb.jpg"),
       device = "jpeg", height = 4.64, width = 6.5, units = "in")
















