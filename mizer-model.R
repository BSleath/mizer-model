# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)


# Setting parameters ------------------------------------------------------

# Set parameters as a single species model
params_initial <- newSingleSpeciesParams()

# Parameters with a reduced resource level
params_rr <- setResource(params_initial,
                      resource_dynamics = "resource_semichemostat")


# Parameters with initial biomass doubled and initial resources halved
params_double <- params_rr
params_double <- setResource(params_double,
                             resource_level = 0.1)
initialN(params_double) <- 2*initialN(params_double)
initialNResource(params_double) <- initialNResource(params_double)/2

# Fishing gear ------------------------------------------------------------

# Set up fishing gear using sigmoid_weight() selectivity function

gear_params(params_double) <- data.frame(
  gear = "gear",
  species = "Target species",
  catchability = 0.3,
  sel_func = "sigmoid_weight",
  sigmoidal_weight = 15,
  sigmoidal_sigma = 5)

params_double <- setFishing(params_double, gear_params = gear_params)

gear_params(params_double)

# Simulation -------------------------------------------------

# Simulate biomass density when initial biomass is doubled and initial resources
# are reduced
sim_double <- project(params_double, t_max = 15, effort = 1)
animateSpectra(sim_double, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# Shows biomass level oscillating --> predator prey relationship

# Extract yield over time dependent on the fishing gear
getYieldGear(sim_double)
plotYieldGear(sim_double)


# Figures -----------------------------------------------------------------


# Figure showing flux over increasing biomass density

# Extract weight and count number
N <- finalN(sim_double)["Target species", , drop = TRUE]
w <- w(params_double)

# Calculate individual growth rate
E_growth <- getEGrowth(params_double)["Target species", , drop = TRUE]
gr <- w * E_growth
flux <- gr * N


flux_data <- data.frame(
  Weight = w,
  Flux = flux
)

# Plot figure
flux_plot <- ggplot(flux_data, 
       aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,105)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux over increasing weight of fish") +
  theme_classic()

flux_plot

ggsave("figures/flux_plot.png",
       plot = flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)


# Plot flux on a log-log axis
flux_plot <- ggplot(flux_data, 
                    aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_x_log10() + 
  scale_y_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux over increasing weight of fish") +
  theme_classic()

flux_plot

ggsave("figures/flux_plot_with_fishing.png",
       plot = flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)



# Changing parameters -----------------------------------------------------
# Make local reductions in resource replenishment rate and resource capacity

# Calculate flux
N <- finalN(sim_double)["Target species", , drop = TRUE]
w <- w(params_double)

E_growth <- getEGrowth(params_double)["Target species", , drop = TRUE]
gr <- w * E_growth
flux <- gr * N


# Plot flux before local reductions (without fishing gear)
initial_flux_data <- data.frame(Weight = w, 
                             Flux = flux)
initial_flux_plot <- ggplot(initial_flux_data,
                         aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux over increasing weight of fish") +
  theme_classic()
initial_flux_plot


# Reduce rr
rr <- resource_rate(params_double)
w_full <- w_full(params_double)
w_full[215:240]
rr[215:240] <- rr[215:240] / 100000000


# Plot resource rate against weight to check reduction
rr_data <- data.frame(
  Weight = params_double@w_full,
  rr = rr
)

rr_plot <- ggplot(rr_data, 
                  aes(x = Weight, y = rr)) +
  geom_smooth() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rr_plot


# Reduce resource capacity
rc <- resource_capacity(params_double)
w_full[215:240]
rc[215:240] <- rc[215:240] / 10000

# Plot resource capacity over weights to see reduction
rc_data <- data.frame(
  Weight = params_double@w_full,
  rc = rc
)

rc_plot <- ggplot(rc_data, 
                  aes(x = Weight, y = rc)) +
  geom_smooth() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rc_plot

params_reduced <- setResource(params_double,
                              resource_capacity = rc,
                              resource_rate = rr,
                              balance =  FALSE)

gear_params(params_reduced)


# Narrow predation kernel
pred_kernel <- getPredKernel(params_reduced)

pred_kernel_reduced <- pred_kernel[, , 215, drop = FALSE]

ggplot(melt(pred_kernel_reduced)) +
  geom_line(aes(x = w_pred, y = value)) +
  scale_x_log10()

select(species_params(params_reduced), beta, sigma)
given_species_params(params_reduced)$sigma <- 0.2
given_species_params(params_reduced)$beta <- 100

getPredKernel(params_reduced)[, , 180, drop = FALSE] %>% 
  melt() %>% 
  ggplot() +
  geom_line(aes(x = w_pred, y = value)) +
  scale_x_log10()


# lower maximum intake rate
species_params(params_reduced)$h
given_species_params(params_reduced)$h <- 30


# Calculate new flux
sim_reduced <- project(params_reduced, t_max = 15, effort = 1)

N_reduced <- finalN(sim_reduced)["Target species", , drop = TRUE]
w <- w(params_reduced)

E_growth_reduced <- getEGrowth(params_reduced)["Target species", , drop = TRUE]
grr <- w * E_growth_reduced
flux_reduced <- grr * N_reduced


reduced_flux_data <- data.frame(Weight = w, 
                             Flux = flux_reduced)

reduced_flux_plot <- ggplot(reduced_flux_data,
                         aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux with locally reduced resource replenishment rate and resource capacity") +
  theme_classic()

initial_flux_plot
reduced_flux_plot

# Save plots
ggsave("figures/flux_plot_initial.png",
       plot = initial_flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

ggsave("figures/flux_plot_reduced.png",
       plot = reduced_flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)
