# Data Analysis Script
# Script Created 11/20/2022 by Aidan Morales
# Modified 10203/2022 by Aidan Morales

# Removes existing objects from the workspace
rm(list = ls())

# Clears the console
cat("\014")

# Loads required libraries
library(sf)
library(ggsn)
library(dplyr)
library(readr)
library(cowplot)
library(ggplot2)
library(rstudioapi)

# Sets the working directory to the location of this script
setwd(dirname(getActiveDocumentContext()$path))

# Reads in the raw biomass and plot data
tree <- read_csv("Data/Level_1/TreeData_clean.csv", show_col_types = FALSE)
plot <- read_csv("Data/Level_1/PlotData_clean.csv", show_col_types = FALSE)

################################################################################
####### 1. Map of Study Sites ##################################################
################################################################################

# Kellogg Forest Boundary
kf_boundary <- st_read("Data/Shapefiles/KF Boundary/Kellogg_Forest_Boundary.shp") %>%
  st_cast(to = "POLYGON")

# Michigan Counties
mi_counties <- st_read("Data/Shapefiles/Counties_(v17a)/Counties_(v17a).shp")

# Kalamazoo County
kzoo_county <- mi_counties %>%
  filter(NAME == "Kalamazoo")

# Michigan Map
p1 <- ggplot() +
  geom_sf(data = mi_counties, fill = "white", color = "black") +
  geom_sf(data = kzoo_county, color = "purple", fill = "purple", show.legend = FALSE) +
  geom_sf(data = kf_boundary, color = "black", fill = "black", show.legend = FALSE) +
  theme_classic() +
  labs(title = "Michigan") +
  theme(plot.title = element_text(vjust = -8, hjust = .75))

# Kalamazoo County Map
p2 <- ggplot() +
  geom_sf(data = kzoo_county, fill = "white") +
  geom_sf(data = kf_boundary, color = "purple", fill = "purple") +
  theme_classic() +
  labs(
    title = "Kalamazoo County",
    x = "",
    y = ""
  ) +
  theme(plot.title = element_text(vjust = -3, hjust = .15))

# Kellogg Forest Map
ylabs <- seq(42.350, 42.375, 0.005)
xlabs <- seq(85.370, 85.350, -0.005)

p3 <- ggplot() +
  geom_sf(data = kf_boundary, aes(color = "shape_col"), fill = "white") +
  geom_point(
    data = tree, aes(x = lon, y = lat, color = "point_col"),
    shape = 16, show.legend = TRUE
  ) +
  theme_classic() +
  labs(
    title = "MSU W.K. Kellogg Experimental Forest",
    x = "",
    y = ""
  ) +
  scale_x_continuous(breaks = waiver(), labels = paste0(xlabs, "°W")) +
  scale_y_continuous(breaks = waiver(), labels = paste0(ylabs, "°N")) +
  scale_color_manual(
    values = c("shape_col" = "purple", "point_col" = "black"),
    labels = c("Harvested Trees", "Property Boundary"),
    name = ""
  ) +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_legend(override.aes = list(linetype = c(0, 1))),
    color = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA)))
  ) +
  north(
    x.min = min(tree$lon), x.max = max(tree$lon),
    y.min = min(tree$lat), y.max = max(tree$lat),
    location = "topleft",
    symbol = 12, scale = .15
  ) +
  scalebar(
    x.min = min(tree$lon), x.max = max(tree$lon),
    y.min = min(42.35), y.max = max(tree$lat),
    dist = .25, dist_unit = "km", st.size = 3,
    transform = TRUE, model = "WGS84", location = "bottomright"
  )

# Combines maps together
map <- ggdraw() +
  draw_plot(p1, x = -0.03, y = .5, width = .5, height = .5) +
  draw_plot(p3, x = 0.15, y = 0, width = .9, height = .9) +
  draw_plot(p2, x = -0.03, y = 0, width = .5, height = .5)

dir.create("Figures", showWarnings = FALSE)
ggsave(plot = map, "Figures/StudySites.svg", width = 30, height = 20, units = "cm")