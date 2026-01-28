# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)


# Setting parameters ------------------------------------------------------

# Set parameters as a single species model
params_single <- newSingleSpeciesParams()

# Parameters with a reduced resource level
params_initial <- setResource(params_single,
                      resource_dynamics = "resource_semichemostat")


# Parameters with initial biomass doubled and initial resources halved
params_double <- params_initial
params_double <- setResource(params_double,
                             resource_level = 0.1)

initialN(params_double) <- 2*initialN(params_double)
initialNResource(params_double) <- initialNResource(params_double)/2

# Set Reproduction parameters
#getReproductionProportion(params_double)
#repro_prop(params_double) <- repro_prop(params_double)/3


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
sim_double <- project(params_double, t_max = 20, effort = 1)
animateSpectra(sim_double, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# Shows biomass level oscillating --> predator prey relationship

# Extract yield over time dependent on the fishing gear
getYieldGear(sim_double)
plotYieldGear(sim_double)



# Changing parameters -----------------------------------------------------
# Make local reductions in resource replenishment rate and resource capacity



# Calculate flux
N <- finalN(sim_double)["Target species", , drop = TRUE]
w <- w(params_double)

E_growth <- getEGrowth(params_double)["Target species", , drop = TRUE]
gr <- w * E_growth
flux <- gr * N


# Plot flux before local reductions
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

ggsave("figures/flux_plot_initial.png",
       plot = initial_flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)



plotSpectra(params_double, power = 2)



# Reduce rr
rr <- resource_rate(params_double)
w_full <- w_full(params_double)
w_full[222:242]
rr[222:242] <- rr[222:242] / 100000000

# Plot resource rate against weight to check reduction
rr_data <- data.frame(
  Weight = params_double@w_full,
  rr = rr)

rr_plot <- ggplot(rr_data, 
                  aes(x = Weight, y = rr)) +
  geom_smooth() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rr_plot


# Reduce resource capacity
rc <- resource_capacity(params_double)
w_full[222:242]
rc[222:242] <- rc[222:242] / 100000000

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
                              balance = FALSE)



# Narrow predation kernel
pred_kernel <- getPredKernel(params_reduced)
pred_kernel

pred_kernel_reduced <- pred_kernel[, 90:98, , drop = FALSE]

ggplot(melt(pred_kernel_reduced)) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(limits = c(1e-4, 100))

select(species_params(params_reduced), beta, sigma)

params <- params_reduced

given_species_params(params)$sigma <- 0.2
given_species_params(params)$beta <- 100

getPredKernel(params)[, 90:98, , drop = FALSE] %>% 
  melt() %>% 
  ggplot() +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(limits = c(1e-4, 100))


plotSpectra(params, power = 2)


# Calculate new flux
sim_reduced <- project(params, t_max = 15, effort = 1)

N_reduced <- finalN(sim_reduced)["Target species", , drop = TRUE]
w <- w(params_reduced)

E_growth_reduced <- getEGrowth(params)["Target species", , drop = TRUE]
grr <- w * E_growth_reduced
flux_reduced <- grr * N_reduced


reduced_flux_data <- data.frame(Weight = w, 
                             Flux = flux_reduced)

# Plot flux on log-log axis

reduced_flux_log_plot <- ggplot(reduced_flux_data,
                         aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux with locally reduced resource replenishment rate and resource capacity") +
  theme_classic()

reduced_flux_log_plot

# Plot with log y-axis

reduced_flux_plot <- ggplot(reduced_flux_data,
                            aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)"),
    title = "Flux with locally reduced resource replenishment rate 
    and resource capacity") +
  theme_classic()

reduced_flux_plot

# Save plots

ggsave("figures/flux_plot_log_reduced.png",
       plot = reduced_flux_log_plot,
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


## Plot flux graph over a limited weight range
# Flux between weight of maturity and maximum weight

limited_flux_data <- data.frame(Weight = w[89:101], 
                                Flux = flux_reduced[89:101])

# Plot on log-log axis

limited_flux_log_plot <- ggplot(limited_flux_data,
                                aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10(limits = c(25, 100)) +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)")) +
  theme_classic()
limited_flux_log_plot


# Plot on log y-axis

limited_flux_plot <- ggplot(limited_flux_data,
                                aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  scale_x_continuous(limits = c(25, 100)) +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)")) +
  theme_classic()
limited_flux_plot


ggsave("figures/limited_flux_log_plot",
       plot = limited_flux_log_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

ggsave("figures/limited_flux_plot",
       plot = limited_flux_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)


## Yield graph
yield <- getYield(sim_reduced)
time <- getTimes(sim_reduced)

yield_data <- data.frame(Time = time, 
                         Yield = yield)

yield_graph <- ggplot(yield_data, 
                     aes(x = Time, y = Target.species)) + 
  geom_line() +
  #scale_y_log10() + 
  #scale_x_log10() +
  labs(x = paste0("Year"),
       y = paste0("Yield (g/year)") ) +
  theme_classic()
  yield_graph


ggsave("figures/yield_plot",
         plot = yield_graph,
         device = "png",
         width = 8,
         height = 6,
         units = "in",
         dpi = 300)
  


## growth rate plot
growth_data <- data.frame(growth = E_growth_reduced,
                          weight = w)
growth_plot <- ggplot(growth_data,
       aes(x = weight,
           y = growth)) +
  geom_smooth() +
  #scale_y_log10() +
  #scale_x_log10() +
  labs(x = paste0("Weight (g)"),
       y = paste0("Growth Rate (g/year)"),
       title = "Growth rate of fish across increasing size classes") +
  theme_classic()
growth_plot

ggsave("figures/growth_rate_plot.png",
       plot = growth_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)


## plot growth rate for larger fish
# plotted from the weight of maturity to maximum weight

limited_growth_data <- data.frame(growth = E_growth_reduced[89:101],
                          weight = w[89:101])
limited_growth_plot <- ggplot(limited_growth_data,
                      aes(x = weight,
                          y = growth)) +
  geom_smooth() +
  #scale_x_log10() +
  #scale_y_log10() +
  labs(x = paste0("Weight (g)"),
       y = paste0("Growth Rate (g/year)")) +
  theme_classic()
limited_growth_plot


ggsave("figures/limited_growth_rate_plot.png",
       plot = limited_growth_plot,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)


## reproduction
getRDD(params)

