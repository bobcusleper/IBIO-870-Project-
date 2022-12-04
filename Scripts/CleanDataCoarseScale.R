# Data Cleaning Script
# Cleans Level 0 Data to Level 1 Data
# Script Created 11/20/2022 by Aidan Morales
# Modified 12/04/2022 by Aidan Morales

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
tree_raw <- read_csv("../Data/Level_0/TreeData_raw.csv", show_col_types = FALSE)
plot_raw <- read_csv("../Data/Level_0/PlotData_raw.csv", show_col_types = FALSE)

# Initial cleaning of the tree data
tree_clean <- tree_raw %>%
  select(
    project, scientific.name, spp.FIA, spp,
    tree, loc, easting, northing, dbh.in, DWtotWOleaf.lb
  ) %>%
  drop_na() %>%
  distinct(easting, northing, .keep_all = TRUE)

# Gets the tree id info to filter the plot data
tree_id <- tree_clean %>%
  select(project, scientific.name, loc, spp.FIA, spp, tree) %>%
  mutate(index = 1)

# Filters and cleans the plot data according to FIA plot sizes
plot_clean <- left_join(plot_raw, tree_id, by = c("project", "spp", "tree", "scientific.name")) %>%
  mutate(dist.ft = case_when(is.na(dist.ft) & tree.no == 1 ~ 0, TRUE ~ dist.ft)) %>%
  filter(index == 1) %>%
  filter(dist.ft <= 24.0) %>%
  select(loc, scientific.name, spp = spp.FIA, tree, tree.no, spp1, dbh.in) %>%
  filter(spp1 %in% unique(tree_clean$spp.FIA))

plot_clean <- left_join(plot_clean, plot_clean %>%
  distinct(spp, tree) %>%
  mutate(plot = 1:n()), by = c("spp", "tree")) %>%
  select(loc, plot, tree.no, spp1, dbh.in)

# Simplifies the tree data
tree_clean <- tree_clean %>%
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

# Converts UTM 17 coordinates to lat lon
tree_utm_17 <- tree_clean %>%
  filter(easting < 300000)

coords <- tree_utm_17 %>%
  vect(geom = c("easting", "northing"), crs = "+proj=utm +zone=17") %>%
  project("+proj=longlat") %>%
  crds(df = TRUE) %>%
  rename("lat" = y, "lon" = x)

# Joins the new coordinates
tree_utm_17 <- bind_cols(tree_utm_17, coords) %>%
  relocate(loc, spp, tree, easting, northing, lat, lon, dbh.in, dwtot.lb)

# Joins the new coordinates together
tree_clean <- bind_rows(tree_utm_16, tree_utm_17)

# Exports the cleaned data
write_csv(tree_clean, "../Data/Level_1/TreeDataCoarse_clean.csv", progress = FALSE)
write_csv(plot_clean, "../Data/Level_1/PlotDataCoarse_clean.csv", progress = FALSE)
