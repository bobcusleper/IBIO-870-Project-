# Data Analysis Script
# Script Created 11/20/2022 by Aidan Morales
# Modified 12/04/2022 by Aidan Morales

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
library(ggpubr)
library(gt)
library(ape)
library(nlme)

# Sets the working directory to the location of this script
setwd(dirname(getActiveDocumentContext()$path))

# Reads in the raw biomass and plot data
tree <- read_csv("../Data/Level_1/TreeDataCoarse_clean.csv", show_col_types = FALSE)
plot <- read_csv("../Data/Level_1/PlotDataCoarse_clean.csv", show_col_types = FALSE)

################################################################################
####### 1. Map of Study Sites ##################################################
################################################################################

# Michigan Counties
mi_counties <- st_read("../Data/Shapefiles/Counties_(v17a)/Counties_(v17a).shp")

# Study Site Map
map <- ggplot() +
  geom_sf(data = mi_counties, fill = "white", color = "black") +
  geom_point(
    data = tree, aes(x = lon, y = lat, color = loc),
    shape = 16, size = 3,
    show.legend = TRUE
  ) +
  theme_classic() +
  labs(title = "Michigan Study Sites", x = "", y = "", color = "Forest Code")

# Exports Map
dir.create("../Figures", showWarnings = FALSE)
ggsave(plot = map, "../Figures/StudySites.svg", width = 15, height = 15, units = "cm")

################################################################################
####### 2. Allometry Plots #####################################################
################################################################################

# Defines Power Function
pwr_fun <- function(x, a, b) {
  a * x^b
}

# Estimates power function coefficients through log transformation
lm <- lm(log(tree$dwtot.lb) ~ log(tree$dbh.in))
a <- exp(coef(lm)[1])
b <- coef(lm)[2]

# Plot labels
label <- vector(mode = "character", length = 3)
label[1] <- "y = a * x^b"
label[2] <- paste0("a = ", round(a, 3))
label[3] <- paste0("b = ", round(b, 3))

# Prediction Lines
x_predicted <- seq(min(tree$dbh.in), max(tree$dbh.in), 0.001)
y_predicted <- pwr_fun(x_predicted, a, b)
y_predicted <- y_predicted[y_predicted <= max(tree$dwtot.lb)]
x_predicted <- x_predicted[1:length(y_predicted)]

# Plots Allometry
p1 <- ggplot() +
  geom_point(data = tree, aes(x = dbh.in, y = dwtot.lb, color = spp), size = 2) +

  # Plots the predicted lines
  geom_line(aes(x = x_predicted, y = y_predicted), linewidth = 1) +
  labs(
    subtitle = "Spatially Unorrected",
    x = "DBH (in)",
    y = "Total Woody Biomass (lb)",
    color = "FIA SpCode"
  ) +

  # Plots the parameter values
  geom_label(
    aes(
      label = paste(label[1], label[2], label[3], sep = "\n"),
      x = min(x_predicted),
      y = max(y_predicted)
    ),
    label.size = .5,
    label.padding = unit(1, "lines"),
    label.r = unit(0, "lines"),
    hjust = "inward",
    vjust = "inward"
  ) +
  theme_classic() +
  theme(legend.position = "right")

################################################################################
####### 3. Spatial Autocorrelation Test ########################################
################################################################################

# Tests for a normal distribution
hist(tree$dwtot.lb)
shapiro.test(tree$dwtot.lb)

hist(log(tree$dwtot.lb))
shapiro.test(log(tree$dwtot.lb))

qqnorm(log(tree$dwtot.lb))
qqline(log(tree$dwtot.lb), col = "red")

# Log transformation to be more normally distributed
tree$log.dwtot.lb <- log(tree$dwtot.lb)

# Compute distance matrix
data.dist <- as.matrix(dist(cbind(tree$easting, tree$northing)))

# Create an inverse distance matrix
w <- 1 / data.dist
diag(w) <- 0

# Moran's I
moran <- Moran.I(tree$log.dwtot.lb, w, scaled = TRUE, na.rm = TRUE, alternative = "two.sided")

# Mantel Test
dwtot.dist <- as.matrix(dist(cbind(tree$dwtot.lb, tree$dwtot.lb)))
mantel <- ape::mantel.test(m1 = data.dist, m2 = dwtot.dist, nperm = 999, alternative = "two.sided")

################################################################################
####### 4. Autoregressive Model  ###############################################
################################################################################

dummy <- rep(1, nrow(tree))
lme <- lme(fixed = log(dwtot.lb) ~ log(dbh.in), data = tree, random = ~ 1 | dummy, method = "ML")
lme <- update(lme, correlation = corGaus(1, form = ~ easting + northing), method = "ML")

a2 <- pull(exp(coef(lme)[1]))
b2 <- pull(coef(lme)[2])

# Plot labels
label2 <- vector(mode = "character", length = 3)
label2[1] <- paste0("a = ", round(a2, 3))
label2[2] <- paste0("b = ", round(b2, 3))

# Prediction Lines
x_predicted <- seq(min(tree$dbh.in), max(tree$dbh.in), 0.001)
y_predicted <- pwr_fun(x_predicted, a2, b2)
y_predicted <- y_predicted[y_predicted <= max(tree$dwtot.lb)]
x_predicted <- x_predicted[1:length(y_predicted)]

# Plots Allometry
p2 <- ggplot() +
  geom_point(data = tree, aes(x = dbh.in, y = dwtot.lb, color = spp), size = 2) +

  # Plots the predicted lines
  geom_line(aes(x = x_predicted, y = y_predicted), linewidth = 1) +
  labs(
    subtitle = "Spatially Corrected",
    x = "DBH (in)",
    y = "Total Woody Biomass (lb)",
    color = "FIA SpCode"
  ) +

  # Plots the parameter values
  geom_label(
    aes(
      label = paste(label2[1], label2[2], sep = "\n"),
      x = min(x_predicted),
      y = max(y_predicted)
    ),
    label.size = .5,
    label.padding = unit(1, "lines"),
    label.r = unit(0, "lines"),
    hjust = "inward",
    vjust = "inward"
  ) +
  theme_classic() +
  theme(legend.position = "right")

# Combines both the uncorrected and corrected allometry plots
p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right") %>% 
  annotate_figure(
    top = text_grob("Coarse Scale Aboveground Biomass Allometry",
                    color = "black", face = "bold", size = 14
    ),
    bottom = text_grob(bquote("Formula: y = " * a * x^b * ""))
  )

# Exports the combined plot
ggsave(plot = p3, "../Figures/AllometryCoarse.svg", width = 40, height = 20, units = "cm")

################################################################################
####### 5. Biomass Calculations ################################################
################################################################################

# Calculates Per Acre Biomass
plot_size <- 24^2 * pi / 43560

plot_mass1 <- plot %>%
  mutate(mass.est.kg = pwr_fun(dbh.in, !!a, !!b)) %>%
  group_by(plot, spp1) %>%
  summarize(mass.est.ac.kg = sum(mass.est.kg, na.rm = TRUE) / plot_size) %>%
  group_by(spp1) %>%
  summarize("Uncorrected" = round(sum(mass.est.ac.kg) / 135))

plot_mass2 <- plot %>%
  mutate(mass.est.kg = pwr_fun(dbh.in, !!a2, !!b2)) %>%
  group_by(plot, spp1) %>%
  summarize(mass.est.ac.kg = sum(mass.est.kg, na.rm = TRUE) / plot_size) %>%
  group_by(spp1) %>%
  summarize("Corrected" = round(sum(mass.est.ac.kg) / 135))

# Welch’s t-test
ttest <- t.test(plot_mass1$Uncorrected, plot_mass2$Corrected)

plot_mass <- left_join(plot_mass1, plot_mass2, by = "spp1") %>%
  rename("FIA SpCode" = spp1) %>%
  summarize(
    "Uncorrected" = sum(Uncorrected),
    "Corrected" = sum(Corrected)
  ) %>%
  mutate("T" = round(ttest$statistic, 3), "DF" = round(ttest$parameter), "P-value" = round(ttest$p.value, 3))

################################################################################
####### 6. Tables ##############################################################
################################################################################

spatial_table <- bind_rows(
  tibble("Test" = "Moran's I", "Result" = moran$observed, "P-value" = moran$p.value),
  tibble("Test" = "Mantel", "Result" = mantel$z.stat, "P-value" = mantel$p)
) %>%
  gt() %>%
  tab_header(
    title = "Table 1: Spatial Autocorrelation Statistics"
  ) %>%
  tab_options(heading.align = "left") %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  opt_table_font(font = "Open Sans")

# Exports the spatial statistics table
gtsave(spatial_table, "Table1.png", "../Figures/")

biomass_table <- plot_mass %>%
  gt() %>%
  tab_header(
    title = "Table 2: AGB (kg) / Acre",
    subtitle = "Biomass Correction for Spatial Autocorrelation"
  ) %>%
  tab_options(heading.align = "left") %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  opt_table_font(font = "Open Sans")

# Exports the biomass table
gtsave(biomass_table, "Table2.png", "../Figures/")
