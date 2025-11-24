# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)


# Setting parameters ------------------------------------------------------

# set parameters as a single species model
params_initial <- newSingleSpeciesParams()

# parameters with a reduced resource level
params_rr <- setResource(params_initial,
                      resource_dynamics = "resource_semichemostat",
                      resource_level = 0.1)


# parameters with initial biomass doubled and initial resources halved
params_double <- params_rr
initialN(params_double) <- 2*initialN(params_double)
initialNResource(params_double) <- initialNResource(params_double)/2

# Fishing gear ------------------------------------------------------------

# set up fishing gear using sigmoid_weight() selectivity function

gear_params(params_double) <- data.frame(
  gear = "gear",
  species = "Target species",
  catchability = 0.2,
  sel_func = "sigmoid_weight",
  sigmoidal_weight = 30,
  sigmoidal_sigma = 5)

params_double <- setFishing(params_double, gear_params = gear_params)

gear_params(params_double)

# Simulation -------------------------------------------------

# simulate biomass density when initial biomass is doubled and resource level is reduced
# increase effort for increased fishing
sim_double <- project(params_double, t_max = 15, effort = 1)
animateSpectra(sim_double, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# shows biomass level oscillating --> predator prey relationship

# extract yield over time dependent on the fishing gear
getYieldGear(sim_double)
plotYieldGear(sim_double)


# Figures -----------------------------------------------------------------


# figure showing flux over increasing biomass density

# Extract weight and count number
N <- finalN(sim_double)["Target species", , drop = TRUE]
w <- w(params_double)

# calculate individual growth rate
E_growth <- getEGrowth(params_double)["Target species", , drop = TRUE]
gr <- w * E_growth
flux <- gr * N


plot_data <- data.frame(
  Weight = w,
  Flux = flux
)

# plot figure
flux_plot <- ggplot(plot_data, 
       aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,105)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux over increasing weight of fish") +
  theme_classic()

ggsave("figures/flux_plot.png",
       plot = flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)
