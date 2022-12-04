# Data Analysis Script
# Script Created 11/20/2022 by Aidan Morales
# Modified 12/03/2022 by Aidan Morales

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
library(doParallel)
library(doSNOW)
library(foreach)
library(ggpubr)

# Sets the working directory to the location of this script
setwd(dirname(getActiveDocumentContext()$path))

# Reads in the raw biomass and plot data
tree <- read_csv("Data/Level_1/TreeData_clean.csv", show_col_types = FALSE)
plot <- read_csv("Data/Level_1/PlotData_clean.csv", show_col_types = FALSE)

################################################################################
####### 1. Map of Study Sites ##################################################
################################################################################

# Michigan Counties
mi_counties <- st_read("Data/Shapefiles/Counties_(v17a)/Counties_(v17a).shp")

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

map

# Exports Map
dir.create("Figures", showWarnings = FALSE)
ggsave(plot = map, "Figures/StudySites.svg", width = 30, height = 20, units = "cm")

################################################################################
####### 2. Allometry Plots #####################################################
################################################################################

# Defines Power Function
pwr_fun <- function(x, a, b) {
  a * x^b
}

# Distinct species to loop over
species <- unique(tree$spp)

# Plot Colors
n <- length(species)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
set.seed(1)
palette <- sample(color, n, replace = FALSE)

# Uses parallel loop to ensure each graph gets its own worker
# and does not only plot the last graph in the loop
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Initializes the progress bar for the loop
progress <- function(n) setTxtProgressBar(txtProgressBar(max = length(species), style = 3), n)
opts <- list(progress = progress)

# Loops through each tree in files and extracts the data
results <- foreach(i = 1:length(species), .inorder = TRUE, .options.snow = opts) %dopar% {
  library(dplyr)
  library(ggplot2)

  # Filters the tree species
  temp <- tree %>%
    filter(spp == species[i])

  # Estimates power function coefficients through log transformation
  lm <- lm(log(temp$dwtot.lb) ~ log(temp$dbh.in))
  a <- exp(coef(lm)[1])
  b <- coef(lm)[2]

  # Plot labels
  label <- vector(mode = "character", length = 3)
  label[1] <- "y = a * x^b"
  label[2] <- paste0("a = ", round(a, 3))
  label[3] <- paste0("b = ", round(b, 3))

  # Prediction Lines
  x_predicted <- seq(min(temp$dbh.in), max(temp$dbh.in), 0.001)
  y_predicted <- pwr_fun(x_predicted, a, b)
  y_predicted <- y_predicted[y_predicted <= max(temp$dwtot.lb)]
  x_predicted <- x_predicted[1:length(y_predicted)]

  # Plots Allometry
  p <- ggplot() +
    geom_point(data = temp, aes(x = dbh.in, y = dwtot.lb), color = palette[i], size = 2) +

    # Plots the predicted lines
    geom_line(aes(x = x_predicted, y = y_predicted), linewidth = 1) +
    labs(
      # title = "Above Ground Biomass Allometry",
      subtitle = temp$scientific.name[i],
      x = "DBH (in)",
      y = "Total Woody Biomass (lb)"
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
    theme(legend.position = "none")

  return(list(p, a, b))
}

# Closes the clusters
stopCluster(cl)

# Vector to store results
a <- vector()
b <- vector()
plots <- list()

for (i in 1:length(results)) {
  plots[[i]] <- results[[i]][[1]]
  a[i] <- results[[i]][[2]]
  b[i] <- results[[i]][[3]]
}

# Table with results
data <- tibble(species, a, b)

# Plots results
ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
  plots[[5]], plots[[6]], plots[[7]], plots[[8]],
  ncol = 4, nrow = 2
) %>%
  annotate_figure(top = text_grob("Allometry - Uncorrected Spatial Autocorrelation",
    color = "black", face = "bold", size = 14
  ))

################################################################################
####### 3. Biomass Calculations ################################################
################################################################################

# Calculates Per Acre Biomass
plot_size <- 24^2 * pi / 43560

plot_mass <- plot %>%
  left_join(data %>% rename("spp1" = species), by = "spp1") %>%
  mutate(mass.est.kg = pwr_fun(dbh.in, a, b)) %>%
  group_by(plot, spp1) %>%
  summarize(mass.est.ac.kg = sum(mass.est.kg, na.rm = TRUE) / plot_size) %>%
  group_by(spp1) %>%
  summarize(mass.est.ac.kg = sum(mass.est.ac.kg) / 135)

biomass_table <- plot_mass %>% 
  rename("FIA SpCode" = spp1, "AGB (kg) / acre" = 2) %>% 
  gt() %>% 
  tab_header(title = md("Table 1: Per Acre Aboveground Biomass")) %>% 
  cols_align(
    align = "auto",
    columns = everything()
  )

################################################################################
####### 4. Spatial Autocorrelation Test ########################################
################################################################################

library(ape)
library(nlme)
  
moran.tests <- list()
mantel.tests <- list()
lmes <- list()

for(i in 1:length(species)){

  # Filters the tree species
  temp <- tree %>%
    filter(spp == species[i])

  # hist(temp$dwtot.lb)
  # shapiro.test(temp$dwtot.lb)
  #
  # hist(log(temp$dwtot.lb))
  # shapiro.test(log(temp$dwtot.lb))
  #
  # qqnorm(log(temp$dwtot.lb))
  # qqline(log(temp$dwtot.lb), col = "red")

  # Log transformation to be more normally distributed
  temp$log.dwtot.lb <- log(temp$dwtot.lb)

  # Compute distance matrix
  data.dist <- as.matrix(dist(cbind(temp$easting, temp$northing)))

  # Create an inverse distance matrix
  w <- 1 / data.dist
  diag(w) <- 0

  # Histogram breaks & counts
  #hist(data.dist)$breaks

  # Moran's I
  moran <- Moran.I(temp$log.dwtot.lb, w, scaled = TRUE, na.rm = TRUE, alternative = "two.sided")

  # Mantel Test
  dwtot.dist <- as.matrix(dist(cbind(temp$dwtot.lb, temp$dwtot.lb)))
  mantel <- ape::mantel.test(m1 = data.dist, m2 = dwtot.dist, nperm = 999, alternative = "two.sided")

  ################################################################################
  ####### 5. Autoregressive Model  ###############################################
  ################################################################################
  
  dummy <- rep(1, nrow(temp))
  lme <- lme(fixed = log(dwtot.lb) ~ log(dbh.in), data = temp, random = ~ 1 | dummy, method = "ML")
  lme <- update(lme, correlation = corGaus(1, form = ~ easting + northing), method = "ML")

  # Saves results
  moran.tests[[i]] <- moran
  mantel.tests[[i]] <- mantel
  lmes[[i]] <- lme
}
  
# Uses parallel loop to ensure each graph gets its own worker
# and does not only plot the last graph in the loop
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Initializes the progress bar for the loop
progress <- function(n) setTxtProgressBar(txtProgressBar(max = length(species), style = 3), n)
opts <- list(progress = progress)
  
# Loops through each tree in files and extracts the data
results2 <- foreach(i = 1:length(species), .inorder = TRUE, .options.snow = opts) %dopar% {
  library(dplyr)
  library(ggplot2)
  
  # Filters the tree species
  temp <- tree %>%
    filter(spp == species[i])

  a2 <- pull(exp(coef(lmes[[i]])[1]))
  b2 <- pull(coef(lmes[[i]])[2])

  # Plot labels
  label <- vector(mode = "character", length = 3)
  label[1] <- "y = a * x^b"
  label[2] <- paste0("a = ", round(a2, 3))
  label[3] <- paste0("b = ", round(b2, 3))

  # Prediction Lines
  x_predicted <- seq(min(temp$dbh.in), max(temp$dbh.in), 0.001)
  y_predicted <- pwr_fun(x_predicted, a2, b2)
  y_predicted <- y_predicted[y_predicted <= max(temp$dwtot.lb)]
  x_predicted <- x_predicted[1:length(y_predicted)]

  # Plots Allometry
  p <- ggplot() +
    geom_point(data = temp, aes(x = dbh.in, y = dwtot.lb), color = palette[i], size = 2) +

    # Plots the predicted lines
    geom_line(aes(x = x_predicted, y = y_predicted), linewidth = 1) +
    labs(
      # title = "Above Ground Biomass Allometry",
      subtitle = temp$scientific.name[i],
      x = "DBH (in)",
      y = "Total Woody Biomass (lb)"
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
    theme(legend.position = "none")

  return(list(p, a2, b2))
}

# Closes the clusters
stopCluster(cl)

# Vector to store results
a2 <- vector()
b2 <- vector()
plots2 <- list()

for (i in 1:length(results)) {
  plots2[[i]] <- results[[i]][[1]]
  a2[i] <- results[[i]][[2]]
  b2[i] <- results[[i]][[3]]
}

# Table with results
data2 <- tibble(species, a2, b2)

# Plots results
ggarrange(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]],
          plots2[[5]], plots2[[6]], plots2[[7]], plots2[[8]],
          ncol = 4, nrow = 2
) %>%
  annotate_figure(top = text_grob("Allometry - Corrected Spatial Autocorrelation",
                                  color = "black", face = "bold", size = 14
  ))

# Moran Table
mpi <- vector()
mpe <- vector()
mpsd <- vector()

for (i in 1:length(species)) {
  mpi[i] <- moran.tests[[i]]$p.value
  mpe[i] <- moran.tests[[i]]$expected
  mpsd[i] <- moran.tests[[i]]$sd
}

library(gt)

tibble(species, mpi, mpe, mpsd) %>% 
  arrange(mpi) %>% 
  rename("FIA SpCode" = species, "P-value" = mpi, "Expected" = mpe, "Standard Deviation" = mpsd) %>% 
  gt() %>% 
  tab_header(title = md("Moran's I Results Per Species")) %>% 
  cols_align(
    align = "auto",
    columns = everything()
  )

# Mantel Table
mtp <- vector()
mtz <- vector()

for (i in 1:length(species)) {
  mtp[i] <- mantel.tests[[i]]$p
  mtz[i] <- mantel.tests[[i]]$z.stat
}

tibble(species, mtp, mtz) %>% 
  arrange(mtp) %>% 
  rename("FIA SpCode" = species, "P-value" = mtp, "Z-Score" = mtz) %>% 
  gt() %>% 
  tab_header(title = md("Mantel Test Per Species")) %>% 
  cols_align(
    align = "auto",
    columns = everything()
  )

