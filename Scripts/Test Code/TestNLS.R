nls_model <- nls(dwtot.lb ~ a * dbh.in^b,
  start = list(b = 1, a = 2),
  algorithm = "port",
  data = tree
)

b2 <- environment(nls_model[["m"]][["incr"]])[["internalPars"]][1]
a2 <- environment(nls_model[["m"]][["incr"]])[["internalPars"]][2]

tree %>%
  ggplot(aes(x = dbh.in, y = dwtot.lb, color = spp)) +
  geom_point(size = 2) +
  geom_smooth(
    method = "nls",
    formula = y ~ a * x^b,
    method.args = list(start = c(a = 1, b = 2)),
    se = FALSE,
    color = "black"
  ) +
  labs(
    title = "Above Ground Biomass Allometry",
    x = "DBH (in)",
    y = "Total Woody Biomass (lb)"
  ) +
  theme_classic() +
  facet_wrap(~spp)