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
lbpath <- "/Users/sametzger/Library/CloudStorage/OneDrive-SharedLibraries-UCB-O365/Thomas Merchant - Plant_Map_Team/lb_maps"

#' read in files
gpkg_files <- list.files(path = lbpath, pattern = "*.gpkg", full.names = TRUE)
maps_all <- map_dfr(gpkg_files, st_read)

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
template_1 <- borders %>% filter(block == "1", plot == "1")  # for 1_12
template_3 <- borders %>% filter(block == "3", plot == "3")  # for 3_4
template_4 <- borders %>% filter(block == "4", plot == "2")  # for 4_1, 4_9, 4_10

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
  group_by(block, plot, year) %>%
  group_map(~ {
    plot_block <- as.character(.y$block)
    plot_num   <- as.character(.y$plot)
    plot_year  <- .y$year
    
    # Debug message for each plot
    message("Processing block ", plot_block, " plot ", plot_num, " year ", plot_year)
    
    border <- borders %>%
      filter(block == plot_block, plot == plot_num)
    
    if (nrow(border) == 0) {
      message("No border found - skipping")
      return(NULL)
    }
    
    tryCatch({
      total_area <- as.numeric(st_area(border))
      occupied <- st_union(.x)
      intersection <- st_intersection(occupied, st_geometry(border))
      occupied_area <- if (length(intersection) == 0 || all(st_is_empty(intersection))) 0 else 
                             as.numeric(st_area(intersection))
      
      data.frame(block = plot_block, plot = plot_num, year = plot_year, total_area = total_area, 
                 occupied_area = occupied_area, litter_area = total_area - occupied_area, 
                 litter_pct = (total_area - occupied_area) / total_area * 100)
    }, error = function(e) {
      message("Error in block ", plot_block, " plot ", plot_num, 
              " year ", plot_year, ": ", e$message)
      return(NULL)
    })
  }, .keep = TRUE) %>%
  bind_rows()

oneTwelve <- maps_all %>% 
  filter(block == 1 & plot == 12 & year == 2023) %>% 
  ggplot() + geom_sf(fill = "red")  
  
  
  