#' 
#' using map shapefiles to quantify bare / litter quantities
#' for the Boulder Drought Project (DroughtNet) years 2018 - 2025
#' SAM + TKM
#' 

#' load packages
library(sf)
library(tidyverse)

#' 
#' **Load data**
#'
#'set path
lbpath <- "/Users/sametzger/Library/CloudStorage/OneDrive-SharedLibraries-UCB-O365/Thomas Merchant - Plant_Map_Team/"

#' read in files
shp_files <- list.files(path = paste0(lbpath, "Data/shapefiles"), pattern = "*.shp", full.names = TRUE)
maps_all <- map_dfr(shp_files, st_read)

#' set crs and borders
#' Read all borders and extract block/plot from filename
#' Get both .shp and .gpkg border files
border_files_shp  <- list.files(path = "/Users/sametzger/Library/CloudStorage/OneDrive-SharedLibraries-UCB-O365/Thomas Merchant - Plant_Map_Team/Digitizing 2023 Maps/2023 Done", pattern = "half_grid.*\\.shp$", recursive = TRUE, full.names = TRUE) 
border_files_gpkg <- list.files(path = "/Users/sametzger/Library/CloudStorage/OneDrive-SharedLibraries-UCB-O365/Thomas Merchant - Plant_Map_Team/Digitizing 2023 Maps/2023 Done", pattern = "half_grid.*\\.gpkg$", recursive = TRUE, full.names = TRUE)

border_files <- c(border_files_shp, border_files_gpkg)
crs <- st_crs(maps_all)

borders <- lapply(border_files, function(f) {
  parts <- str_match(basename(f), "half_grid_(\\d+)_(\\d+)")
  if (is.na(parts[2]) || is.na(parts[3])) return(NULL)
  
#' Skip corrupted files: 1_12, 3_4, 4_1, 4_9, and 4_10
  tryCatch({st_read(f, quiet = TRUE) %>%
      st_union() %>%
      st_as_sf() %>%
      rename(geometry = x) %>%
      mutate(block = parts[2], plot = parts[3])
  }, error = function(e) {
    message("Skipping corrupted file: ", basename(f))
    return(NULL)
  })
}) %>%
  Filter(Negate(is.null), .) %>%
  do.call(rbind, .) %>%
  st_set_crs(crs) %>%
  distinct(block, plot, .keep_all = TRUE)

#' add half_grids for missing, corrupted files
template_1 <- borders %>% filter(block == "1", plot == "10")  # for 1_12
template_3 <- borders %>% filter(block == "3", plot == "3")  # for 3_4
template_4 <- borders %>% filter(block == "4", plot == "3")  # for 4_1, 4_2, 4_9, 4_10

#' Create missing rows
missing_borders <- bind_rows(
  template_1 %>% mutate(plot = "12"),
  template_3 %>% mutate(plot = "4"),
  template_4 %>% mutate(plot = "1"),
  template_4 %>% mutate(plot = "9"),
  template_4 %>% mutate(plot = "10"))

#' Add to borders
borders <- bind_rows(borders, missing_borders) %>%
  arrange(block, plot)

#' check that borders for all plots are there - should be 72
nrow(borders)

#' Create dataframe for calculating litter
litter_df <- maps_all %>%
  #exclude quads where litter was mapped as a species
  filter(!(year == 2023 & quad %in% c("1_1", "2_1", "4_1", "4_8")),
         !(year == 2024 & quad == "6_5")) %>% 
  group_by(block, plot, year) %>%
  group_map(~ {
    plot_block <- as.character(.y$block)
    plot_num   <- as.character(.y$plot)
    plot_year  <- .y$year
    
    message("Processing block ", plot_block, " plot ", plot_num, " year ", plot_year)
    
    border <- borders %>% filter(block == plot_block, plot == plot_num)
    if (nrow(border) == 0) {
      message("No border found - skipping")
      return(NULL)
    }
    tryCatch({
      total_area <- as.numeric(st_area(border))
      
      # Helper function to get clipped area for a subset of species
      get_area <- function(data) {
        if (nrow(data) == 0) return(0)
        unioned <- st_union(data)
        clipped <- st_intersection(unioned, st_geometry(border))
        if (length(clipped) == 0 || all(st_is_empty(clipped))) 0 else
          as.numeric(st_area(clipped))
      }
      
      #' Split into categories
      bare_data  <- .x %>% filter(species == "bare")
      other_data <- .x %>% filter(species %in% c("rock", "cow"))
      veg_data   <- .x %>% filter(!species %in% c("bare", "rock", "cow"))
      all_data   <- .x  # 
      
      # Calculate areas
      bare_area  <- get_area(bare_data)
      other_area <- get_area(other_data)
      veg_area   <- get_area(veg_data)
      
      # Total occupied 
      occupied_area <- get_area(all_data)
      litter_area   <- total_area - occupied_area
      
      data.frame(block = plot_block, plot = plot_num, year = plot_year, 
                 total_area = total_area, occupied_area = occupied_area, litter_area = litter_area, 
                 litter_pct = litter_area / total_area * 100, bare_area = bare_area, 
                 bare_pct = bare_area / total_area * 100, other_area = other_area, 
                 other_pct = other_area / total_area * 100, veg_area = veg_area, 
                 veg_pct = veg_area / total_area * 100)
    }, error = function(e) {
      message("Error in block ", plot_block, " plot ", plot_num,
              " year ", plot_year, ": ", e$message)
      return(NULL)
    })
  }, .keep = TRUE) %>%
  bind_rows()

#' Calculate for excluded plots separately
excluded_df <- maps_all %>%
  filter((year == 2023 & quad %in% c("1_1", "2_1", "4_1", "4_8")) |
          (year == 2024 & quad == "6_5")) %>%
  filter(species != "litter") %>%  # remove litter polygons
  group_by(block, plot, year) %>%
  group_map(~ {
    plot_block <- as.character(.y$block)
    plot_num   <- as.character(.y$plot)
    plot_year  <- .y$year
    
    message("Processing block ", plot_block, " plot ", plot_num, " year ", plot_year)
    
    border <- borders %>% filter(block == plot_block, plot == plot_num)
    if (nrow(border) == 0) {
      message("No border found - skipping")
      return(NULL)
    }
    tryCatch({
      total_area <- as.numeric(st_area(border))
      
      # Helper function to get clipped area for a subset of species
      get_area <- function(data) {
        if (nrow(data) == 0) return(0)
        unioned <- st_union(data)
        clipped <- st_intersection(unioned, st_geometry(border))
        if (length(clipped) == 0 || all(st_is_empty(clipped))) 0 else
          as.numeric(st_area(clipped))
      }
      
      #' Split into categories
      bare_data  <- .x %>% filter(species == "bare")
      other_data <- .x %>% filter(species %in% c("rock", "cow"))
      veg_data   <- .x %>% filter(!species %in% c("bare", "rock", "cow"))
      all_data   <- .x  # 
      
      # Calculate areas
      bare_area  <- get_area(bare_data)
      other_area <- get_area(other_data)
      veg_area   <- get_area(veg_data)
      
      # Total occupied 
      occupied_area <- get_area(all_data)
      litter_area   <- total_area - occupied_area
      
      data.frame(block = plot_block, plot = plot_num, year = plot_year, 
                 total_area = total_area, occupied_area = occupied_area, litter_area = litter_area, 
                 litter_pct = litter_area / total_area * 100, bare_area = bare_area, 
                 bare_pct = bare_area / total_area * 100, other_area = other_area, 
                 other_pct = other_area / total_area * 100, veg_area = veg_area, 
                 veg_pct = veg_area / total_area * 100)
    }, error = function(e) {
      message("Error in block ", plot_block, " plot ", plot_num,
              " year ", plot_year, ": ", e$message)
      return(NULL)
    })
  }, .keep = TRUE) %>%
  bind_rows()

#' Combine
litter_df <- bind_rows(litter_df, excluded_df) %>%
  arrange(block, plot, year)

  
  
  