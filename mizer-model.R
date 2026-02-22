# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)
library(plotly)
library(dplyr)


# Setting parameters ------------------------------------------------------

# Set parameters as a single species model
params_single <- newSingleSpeciesParams()

# Parameters with a reduced resource level
params_initial <- setResource(params_single,
                      resource_dynamics = "resource_semichemostat")


# Parameters with initial biomass doubled and initial resources halved
params_double <- params_initial
params_double <- setResource(params_double,
                             resource_level = 0.05)

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
  catchability = 0.8,
  sel_func = "sigmoid_weight",
  sigmoidal_weight = 15,
  sigmoidal_sigma = 5)

params_double <- setFishing(params_double, gear_params = gear_params)

gear_params(params_double)


# Simulation -------------------------------------------------

# Simulate biomass density when initial biomass is doubled and initial resources
# are reduced
sim_double <- project(params_double, t_max = 20, effort = 1)
animateSpectra(sim_double, total = FALSE, power = 2, 
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
w_full[192:232]
rr[192:232] <- rr[192:232] / 2

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
w_full[202:222]
rc[192:232] <- rc[192:232] / 5

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

pred_kernel_reduced <- pred_kernel[, 89, , drop = FALSE]

ggplot(melt(pred_kernel_reduced)) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(limits = c(1e-4, 100))

select(species_params(params_reduced), beta, sigma)

params <- params_reduced

given_species_params(params)$sigma <- 0.8
given_species_params(params)$beta <- 1000

getPredKernel(params)[, 89, , drop = FALSE] %>% 
  melt() %>% 
  ggplot() +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(limits = c(1e-4, 100))




# Calculate new flux
sim_reduced <- project(params, t_max = 30, effort = 1)

N_reduced <- finalN(sim_reduced)["Target species", , drop = TRUE]
w <- w(params_reduced)

E_growth_reduced <- getEGrowth(params)["Target species", , drop = TRUE]
grr <- w * E_growth_reduced
flux_reduced <- grr * N_reduced


reduced_flux_data <- data.frame(Weight = w, 
                             Flux = flux_reduced)
reduced_flux_data$species = "Target Species"

# Plot flux on log-log axis

flux_plot <- ggplot(reduced_flux_data, aes(x = Weight, y = Flux))

plot_ly(reduced_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~species) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                              title_text = "Flux (g/year)"),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)"))

reduced_flux_log_plot <- ggplot(reduced_flux_data,
                         aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)")) +
  theme_classic()

reduced_flux_log_plot

# Plot with log y-axis

reduced_flux_plot <- ggplot(reduced_flux_data,
                            aes(x = Weight, y = Flux)) +
  geom_smooth() +
  scale_y_log10() +
  labs(
    x = paste0("Weight (g)"),
    y = paste0("Flux (g/year)")) +
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
limited_flux_data$Species = "Target Species"

plot_ly(limited_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~Species) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (g/year)"),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)"))

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
growth_data$Species = "Target Species"

plot_ly(growth_data) |> 
  add_lines(x = ~weight, y = ~growth, color = ~Species) |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
                      title_text = "Growth Rate (g/year)"),
         xaxis = list(#type = "log", exponentformat = "power", 
                      title_text = "Weight (g)"))

growth_plot <- ggplot(growth_data,
       aes(x = weight,
           y = growth)) +
  geom_smooth() +
  #scale_y_log10() +
  #scale_x_log10() +
  labs(x = paste0("Weight (g)"),
       y = paste0("Growth Rate (g/year)")) +
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
limited_growth_data$species = "Target Species"

plot_ly(limited_growth_data) |> 
  add_lines(x = ~weight, y = ~growth, color = ~species) |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
    title_text = "Growth Rate (g/year)"),
    xaxis = list(#type = "log", exponentformat = "power", 
      title_text = "Weight (g)"))


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


plotSpectra(sim_reduced, power = 2, wlim = c(1e-8, NA), ylim = c(1e-5, NA),
            time_range = 30)


# animation 
nf <- melt(sim_reduced@n)
n_ppf <- melt(sim_reduced@n_pp)
n_ppf$sp <- "Resource"
nf <- rbind(nf, n_ppf)

plot_ly(nf) %>%
  # show only part of plankton spectrum
  filter(w > 10^-5) %>% 
  # start at time 20
  #filter(time >= 26) %>% 
  # calculate biomass density with respect to log size
  mutate(b = value * w^2) %>% 
  # Plot lines
  add_lines(
    x = ~w, y = ~b,
    color = ~sp,
    frame = ~time,
    line = list(simplify = FALSE)
  ) %>% 
  # Use logarithmic axes
  layout(params, xaxis = list(type = "log", exponentformat = "power",
                         title_text = "body mass (g)"),
         yaxis = list(type = "log", exponentformat = "power",
                      title_text = "biomass (g/m^3)",
                      range = c(-8, 0)))



# animation of flux
species <- "Target species"

flux_series <- lapply(seq_along(time), function(i) {
  n_at_t <- sim_reduced@n[i, species, ]
  data.frame(
    Time = time[i],
    Weight = w,
    Flux = n_at_t * grr
  )
}) %>% bind_rows()

flux_series$Species <- species

plot_ly(flux_series) %>% 
  add_lines(x ~Weight, 
            y ~Flux, 
            color = ~Species, 
            frame = ~Time, 
            line = list(simplify = FALSE)) %>%  
  layout(params, xaxis = list(type = "log", exponentformat = "power",
                         title_text = "body mass (g)"),
         yaxis = list(type = "log", exponentformat = "power",
                      title_text = "flux (g/year)", 
                      range = c(-8, 2)))

