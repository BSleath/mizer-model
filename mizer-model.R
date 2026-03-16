# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)
library(reticulate)
library(scales)
library(processx)

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


# Fishing gear ------------------------------------------------------------

# Set up fishing gear using sigmoid_weight() selectivity function

gear_params(params_double) <- data.frame(
  gear = "gear",
  species = "Target species",
  catchability = 1,
  sel_func = "sigmoid_weight",
  sigmoidal_weight = 25,
  sigmoidal_sigma = 10)

params_double <- setFishing(params_double, gear_params = gear_params)

gear_params(params_double)


# Simulation -------------------------------------------------

# Simulate biomass density when initial biomass is doubled and initial resources
# are reduced
sim_double <- project(params_double, t_max = 50, effort = 0.125)
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

plot_ly(initial_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))

# limited weight range

w[89:101]

initial_limited_flux <- data.frame(Weight = w[61:101], 
                                Flux = flux[61:101])

plot_ly(initial_limited_flux) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))



#initial_flux_plot <- ggplot(initial_flux_data, aes(x = Weight, y = Flux)) +
#  geom_line() +
#  scale_x_log10() +
#  scale_y_log10() +
#  labs(
#    x = paste0("Weight (g)"),
#    y = paste0("Flux (gm<sup>-2</sup>year<sup>-1</sup>)")) +
#  theme_test()

#plotly::ggplotly(initial_flux_plot) 

#ggsave("figures/flux_plot_initial.png",
#       plot = initial_flux_plot,
#       device = "png",
#       width = 8,
#       height = 6,
#       units = "in",
#       dpi = 300)



# Reduce rr
rr <- resource_rate(params_double)
w_full <- w_full(params_double)
w_full[182:242]
rr[182:242] <- rr[182:242] / 2


# Plot resource rate against weight to check reduction
rr_data <- data.frame(
  Weight = params_double@w_full,
  rr = rr)

rr_plot <- ggplot(rr_data, 
                  aes(x = Weight, y = rr)) +
  geom_line() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rr_plot


# Reduce resource capacity
rc <- resource_capacity(params_double)
w_full[182:242]
rc[182:242] <- rc[182:242] / 5

# Plot resource capacity over weights to see reduction
rc_data <- data.frame(
  Weight = params_double@w_full,
  rc = rc
)

rc_plot <- ggplot(rc_data, 
                  aes(x = Weight, y = rc)) +
  geom_line() +
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

large_pred_kernel <- ggplot(melt(pred_kernel_reduced)) +
                            geom_line(aes(x = w_prey, y = value)) +
                            scale_x_log10("Weight of prey (g)", limits = c(1e-5, 1e03)) +
                            scale_y_continuous("Proportion", limits = c(0,1)) +
                            theme_bw()



select(species_params(params_reduced), beta, sigma)

params <- params_reduced

given_species_params(params)$sigma <- 0.8
given_species_params(params)$beta <- 1000

getPredKernel(params)[, 89, , drop = FALSE] %>% 
  melt() %>% 
  ggplot() +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10("Weight of prey (g)", limits = c(1e-5, 1e03)) +
  scale_y_continuous("Proportion")+
  theme_bw()


## reducing feeding level
#feeding <- getFeedingLevel(params)
#w[81:88]
#feeding[81:88] <- feeding[81:88] / 5
#feeding


# Calculate new flux
sim_reduced <- project(params, t_max = 50, effort = 0.75)

N_reduced <- finalN(sim_reduced)["Target species", , drop = TRUE]
w <- w(params)

E_growth_reduced <- getEGrowth(params)["Target species", , drop = TRUE]
grr <- w * E_growth_reduced
flux_reduced <- grr * N_reduced


reduced_flux_data <- data.frame(Weight = w, 
                             Flux = flux_reduced)

plot_ly(reduced_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))


## Plot flux graph over a limited weight range
# Flux between weight of maturity and maximum weight

limited_flux_data <- data.frame(Weight = w[61:101], 
                                Flux = flux_reduced[61:101])

plot_ly(limited_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))

  


## growth rate plot
growth_data <- data.frame(growth = E_growth_reduced,
                          weight = w)

plot_ly(growth_data) |> 
  add_lines(x = ~weight, y = ~growth) |> 
  layout(yaxis = list(title_text = "Growth Rate (g/year)", 
                      range = c(0, 20),
                      showline = TRUE, showgrid = FALSE, zeroline = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))



## plot growth rate for larger fish
# plotted from the weight of maturity to maximum weight

limited_growth_data <- data.frame(growth = E_growth_reduced[89:101],
                          weight = w[89:101])

plot_ly(limited_growth_data) |> 
  add_lines(x = ~weight, y = ~growth) |> 
  layout(yaxis = list(title_text = "Growth Rate (g/year)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE, zeroline = FALSE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))


# animation 
nf <- melt(sim_reduced@n)
n_ppf <- melt(sim_reduced@n_pp)
n_ppf$sp <- "Resource"
nf <- rbind(nf, n_ppf)

plot_ly(nf) %>%
  # show only part of plankton spectrum
  filter(w > 10^-5) %>% 
  # start at time 20
  #filter(time >= 50) %>% 
  #filter(time <= 100) %>% 
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
time_flux <- getTimes(sim_reduced)

flux_series <- lapply(seq_along(time), function(i) {
  n_at_t <- sim_reduced@n[i, species, ]
  data.frame(
    Time = time[i],
    Weight = w,
    Flux = n_at_t * grr)
}) %>% bind_rows()

flux_series$Species <- species

plot_ly(flux_series) %>% 
  filter(time >= 50) %>% 
  #filter(time <= 100) %>%
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


## altered fishing effort
# flux with lower fishing effort than MSY

sim_lower <- project(params, t_max = 50, effort = 0.1)

N_lower <- finalN(sim_lower)["Target species", , drop = TRUE]
flux_lower <- grr * N_lower

lower_flux_data <- data.frame(Weight = w, 
                                Flux = flux_lower)

plot_ly(lower_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))


# flux with higher fishing effort

sim_higher <- project(params, t_max = 50, effort = 1.5)

N_higher <- finalN(sim_higher)["Target species", , drop = TRUE]
flux_higher <- grr * N_higher

higher_flux_data <- data.frame(Weight = w, 
                                Flux = flux_higher)
plot_ly(higher_flux_data) |> 
  add_lines(x = ~Weight, y = ~Flux) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))


# compare plots on same axis

reduced_flux_data$Effort = "MSY"
lower_flux_data$Effort = "Lower fishing effort"
higher_flux_data$Effort = "Higher fishing effort"


fishing_flux <- rbind(higher_flux_data, lower_flux_data, reduced_flux_data)

plot_ly(fishing_flux) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~Effort) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 50))


# biomass ratio between fishing effort

biomass_reduced <- N_reduced * w
biomass_medium <- N_lower * w
biomass_heavy <- N_higher * w

# biomass ratio - low and medium

biomass_ratio_lm <- biomass_medium / biomass_reduced

biomass_ratio_lm_data <- data.frame(Weight = w(params),
                                    Ratio = biomass_ratio_lm)

plot_ly(biomass_ratio_lm_data) |> 
  add_lines(x = ~Weight, y = ~Ratio, type = "scatter") |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
    title_text = "Biomass ratio"),
    xaxis = list(#type = "log", exponentformat = "power", 
      title_text = "Weight (g)"))

# biomass ratio - medium and high

biomass_ratio_mh <- biomass_heavy / biomass_medium

biomass_ratio_mh_data <- data.frame(Weight = w(params),
                                    Ratio = biomass_ratio_mh)

plot_ly(biomass_ratio_mh_data) |> 
  add_lines(x = ~Weight, y = ~Ratio, type = "scatter") |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
    title_text = "Biomass ratio"),
    xaxis = list(#type = "log", exponentformat = "power", 
      title_text = "Weight (g)"))

# biomass ratio - low and high
biomass_ratio_lh <- biomass_heavy / biomass_reduced

biomass_ratio_lh_data <- data.frame(Weight = w(params),
                                 Ratio = biomass_ratio_lh)

plot_ly(biomass_ratio_lh_data) |> 
  add_lines(x = ~Weight, y = ~Ratio, type = "scatter") |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
                      title_text = "Biomass ratio"),
         xaxis = list(#type = "log", exponentformat = "power", 
                      title_text = "Weight (g)"))


# Yield graphs ------------------------------------------------------------

## Yield graph for MSY
sim_MSY <- project(params, t_max = 150, effort = 0.75)

yield <- getYield(sim_MSY)
time <- getTimes(sim_MSY)

yield_data <- data.frame(Time = time, 
                         Yield = yield)

plot_ly(yield_data) |> 
  filter(time >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))


# yield with medium fishing effort
sim_low_yield <- project(params, t_max = 150, effort = 0.1)
yield_lower <- getYield(sim_low_yield)

yield_lower_data <- data.frame(Time = time, 
                               Yield = yield_lower)

plot_ly(yield_lower_data) |> 
  filter(time >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))


## Yield graph with high fishing effort
sim_high_yield <- project(params, t_max = 150, effort = 1.5)
yield_heavy <- getYield(sim_high_yield)

yield_heavy_data <- data.frame(Time = time, 
                         Yield = yield_heavy)

plot_ly(yield_heavy_data) |> 
  filter(time >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))



## time dependent fishing effort

t_total <- 500
low_effort <- 0.75
high_effort <- 1.5
threshold_value <- 0.000025
current_effort <- low_effort


sim_combined <- project(params, t_max = 1, effort = current_effort)

last_yield <- getYield(sim_combined)[idxFinalT(sim_combined), "Target species"]
if (last_yield >= threshold_value) {
  current_effort <- high_effort
} else {
  current_effort <- low_effort
}

for (t in 2:t_total) {
  sim_combined <- project(sim_combined, t_max = 1, effort = current_effort, append = TRUE)
  
  last_yield <- getYield(sim_combined)[idxFinalT(sim_combined), "Target species"]
  
  if (last_yield >= threshold_value) {
    current_effort <- high_effort
  } else {
    current_effort <- low_effort
  }
}


# changing effort over time
effort_df <- melt(getEffort(sim_combined))

plot_ly(effort_df) |> 
  add_lines(x = ~time, y = ~value) |> 
  layout(yaxis = list(title_text = "Effort",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Time (year)", range = c(300, 500),
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))


# Yield with changing effort

yield_fishing <- getYield(sim_combined)
time_fishing <- getTimes(sim_combined)

yield_fishing_data <- data.frame(Time = time_fishing, 
                                Yield = yield_fishing)

plot_ly(yield_fishing_data) |> 
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Year", range = c(300, 500),
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))






gf_original <- melt(getEGrowth(params_initial))
gf_original$Model <- "Original"
gf_starved <- melt(getEGrowth(params))
gf_starved$Model <- "Less prey"
gf_total <- rbind(gf_original, gf_starved)


plot_ly(gf_total) |> 
  add_lines(x = ~w, y = ~value, color = ~Model) |> 
  layout(yaxis = list(title_text = "Growth rate (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE, zeroline = FALSE,
                      range = c(0,40)),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 50))




plotSpectra(sim_reduced, power = 2, wlim = c(1e-8, NA), ylim = c(1e-8, NA),
            time_range = 50)+
  theme_test()



plotYieldVsF(
  params,
  'Target species',
  F_range,
  F_max = 1,
  F_min = 0,
  no_steps = 25,
  distance_func = distanceSSLogN,
  tol = 0.0001,
  t_max = 250) +
  theme_test()

getSelectivity(params)


plotYieldVsF(
  params_initial,
  'Target species',
  F_range,
  F_max = 1,
  F_min = 0,
  no_steps = 25,
  distance_func = distanceSSLogN,
  tol = 0.001,
  t_max = 250) +
  theme_test()
