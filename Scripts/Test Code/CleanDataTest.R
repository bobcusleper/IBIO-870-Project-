# Data Cleaning Script
# Cleans Level 0 Data to Level 1 Data
# Script Created 11/20/2022 by Aidan Morales
# Modified 12/03/2022 by Aidan Morales

# Removes existing objects from the workspace
rm(list = ls())

# Clears the console
cat("\014")

# Loads required libraries
library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(rstudioapi)

# Sets the working directory to the location of this script
setwd(dirname(getActiveDocumentContext()$path))

# Reads in the raw biomass and plot data
tree_raw <- read_csv("Data/Level_0/TreeData_raw.csv", show_col_types = FALSE)
plot_raw <- read_csv("Data/Level_0/PlotData_raw.csv", show_col_types = FALSE)

# Checks forests for complete data sets of 16+ trees
check <- tree_raw %>%
  group_by(loc, spp) %>%
  summarize(n = n()) %>%
  filter(n >= 16)
check

# Selects forests & trees of interest
tree <- tree_raw %>%
  filter((spp == "WP" & loc == "ASGA") |
    (spp == "BF" & loc == "FBIC") |
    (spp == "RM" & loc == "KF") |
    (spp == "HE" & loc == "LC") |
    (spp == "WO" & loc == "MC") |
    (spp == "TP" & loc == "NW") |
    (spp == "BE" & loc == "RF") |
    (spp == "SM" & loc == "RF"))

# Initial cleaning of the tree data
tree_clean <- tree %>%
  select(
    project, scientific.name, spp.FIA, spp,
    tree, loc, easting, northing, dbh.in, DWtotWOleaf.lb
  ) %>%
  drop_na() %>%
  distinct(easting, northing, .keep_all = TRUE) %>%
  mutate(loc = case_when(
    loc == "ASGA" ~ "Allegan State Game Area",
    loc == "FBIC" ~ "Forest Biomass Inovation Center",
    loc == "KF" ~ "Kellogg Forest",
    loc == "LC" ~ "Lake City",
    loc == "RF" ~ "Fred Russ Forest",
    loc == "MC" ~ "MacCready Reserve",
    loc == "NW" ~ "Newton Woods"
  ))

# # Selects forest location with the largest sample size
# max_sample_loc <- tree_clean %>%
#   group_by(loc) %>%
#   summarize(count = n()) %>%
#   slice_max(count) %>%
#   select(loc) %>%
#   pull()
#
# # Filters the tree data to the forest of interest
# tree_clean <- tree_clean %>%
#   filter(loc == max_sample_loc) %>%
#   filter(!spp.FIA == "YB") %>%
#   filter(!(spp.FIA == "LZ" & tree == 2))

# Gets the tree id info to filter the plot data
tree_id <- tree_clean %>%
  select(project, loc, spp.FIA, spp, tree) %>%
  mutate(index = 1)

# Filters and cleans the plot data according to FIA plot sizes
plot_clean <- left_join(plot_raw, tree_id, by = c("project", "spp", "tree")) %>%
  mutate(dist.ft = case_when(is.na(dist.ft) & tree.no == 1 ~ 0, TRUE ~ dist.ft)) %>%
  filter(index == 1) %>%
  filter(dist.ft <= 24.0) %>%
  select(loc, spp = spp.FIA, tree, tree.no, spp1, dbh.in) %>% 
  filter(spp1 %in% unique(tree_clean$spp.FIA))

plot_clean <- left_join(plot_clean, plot_clean %>%
  distinct(spp, tree) %>%
  mutate(plot = 1:n()), by = c("spp", "tree")) %>%
  select(loc, plot, tree.no, spp1, dbh.in)

# Simplifies the tree data
tree_clean <- tree_clean %>%
  mutate(scientific.name = case_when(
    spp == "WP" ~ "Pinus strobus",
    spp == "BF" ~ "Abies balsamea",
    spp == "HE" ~ "Tsuga canadensis",
    spp == "RM" ~ "Acer rubrum",
    spp == "SM" ~ "Acer saccharum",
    spp == "WO" ~ "Quercus alba",
    spp == "TP" ~ "Liriodendron tulipifera",
    spp == "BE" ~ "Fagus grandifolia",
  )) %>% 
  select(loc, scientific.name, spp = spp.FIA, tree, easting, northing, dbh.in, dwtot.lb = DWtotWOleaf.lb) %>%
  arrange(scientific.name, tree)

# Converts UTM 16 coordinates to lat lon
tree_utm_16 <- tree_clean %>%
  filter(easting > 300000)

coords <- tree_utm_16 %>%
  vect(geom = c("easting", "northing"), crs = "+proj=utm +zone=16") %>%
  project("+proj=longlat") %>%
  crds(df = TRUE) %>%
  rename("lat" = y, "lon" = x)

# Joins the new coordinates
tree_utm_16 <- bind_cols(tree_utm_16, coords) %>%
  relocate(loc, scientific.name, spp, tree, easting, northing, lat, lon, dbh.in, dwtot.lb)

# # Converts UTM 17 coordinates to lat lon
# tree_utm_17 <- tree_clean %>%
#   filter(easting < 300000)

# coords <- tree_utm_17 %>%
#   vect(geom = c("easting", "northing"), crs = "+proj=utm +zone=17") %>%
#   project("+proj=longlat") %>%
#   crds(df = TRUE) %>%
#   rename("lat" = y, "lon" = x)

# # Joins the new coordinates
# tree_utm_17 <- bind_cols(tree_utm_17, coords) %>%
#   relocate(loc, spp, tree, easting, northing, lat, lon, dbh.in, dwtot.lb)

# Joins tree data back together and adds 
tree_clean <- tree_utm_16 #bind_rows(tree_utm_16, tree_utm_17)

# Exports the cleaned data
write_csv(tree_clean, "Data/Level_1/TreeData_clean.csv", progress = FALSE)
write_csv(plot_clean, "Data/Level_1/PlotData_clean.csv", progress = FALSE)
