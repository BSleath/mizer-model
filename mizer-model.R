# Library -----------------------------------------------------------------

library(mizer)
library(mizerExperimental)
library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)

# Setting parameters ------------------------------------------------------

# Set parameters as a single species model
params_single <- newSingleSpeciesParams() 
params_single <- setResource(params_single, 
                             resource_dynamics = "resource_semichemostat")

params_initial <- setResource(params_single,
                              resource_level = 0.05)

initialN(params_initial) <- 2*initialN(params_initial)
initialNResource(params_initial) <- initialNResource(params_initial)/2


# Fishing gear ------------------------------------------------------------

# Set up fishing gear using sigmoid_weight() selectivity function

gear_params(params_initial) <- data.frame(
  gear = "gear",
  species = "Target species",
  catchability = 1,
  sel_func = "sigmoid_weight",
  sigmoidal_weight = 25,
  sigmoidal_sigma = 10)

params_initial <- setFishing(params_initial, gear_params = gear_params)

gear_params(params_initial)


# Initial model -------------------------------------------------

# Calculate MSY for initial single species model

plotYieldVsF(
  params_single,
  'Target species',
  F_range,
  F_max = 1,
  F_min = 0,
  no_steps = 25,
  distance_func = distanceSSLogN,
  tol = 0.001,
  t_max = 250) +
  theme_test()

# Plot flux before reductions in resource dynamics
sim_initial <- project(params_initial, t_max = 34, t_save = 0.1, effort = 0.125)

N <- finalN(sim_initial)["Target species", , drop = TRUE]
w <- w(params_initial)

E_growth <- getEGrowth(params_initial)["Target species", , drop = TRUE]
gr <- w * E_growth
flux <- gr * N

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



# Changing parameters -----------------------------------------------------

# Make local reductions in resource replenishment rate and resource capacity

# Reduce rr
rr <- resource_rate(params_initial)
w_full <- w_full(params_initial)
w_full[182:242]
rr[182:242] <- rr[182:242] / 2


# Plot resource rate against weight to check reduction
rr_data <- data.frame(
  Weight = params_initial@w_full,
  rr = rr)

rr_plot <- ggplot(rr_data, 
                  aes(x = Weight, y = rr)) +
  geom_line() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rr_plot


# Reduce resource capacity
rc <- resource_capacity(params_initial)
w_full[182:242]
rc[182:242] <- rc[182:242] / 5

# Plot resource capacity over weights to see reduction
rc_data <- data.frame(Weight = params_initial@w_full,
                      rc = rc)

rc_plot <- ggplot(rc_data, 
                  aes(x = Weight, y = rc)) +
  geom_line() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()
rc_plot

params_reduced <- setResource(params_initial,
                              resource_capacity = rc,
                              resource_rate = rr,
                              balance = FALSE)


# Narrow predation kernel 
pred_kernel <- getPredKernel(params_reduced)

pred_kernel_reduced <- pred_kernel[, 89, , drop = FALSE]

ggplot(melt(pred_kernel_reduced)) +
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


## Calculate MSY after changed parameters

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



# Calculate new flux
sim_reduced <- project(params, t_max = 34, t_save = 0.1, effort = 4.6)

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


## comparative flux on same graph
reduced_flux_data$Sim = "Altered"
initial_flux_data$Sim = "Original"

models_flux <- rbind(reduced_flux_data, initial_flux_data)

plot_ly(models_flux) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~Sim, linetype = ~Sim) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE, range = c(-12,-2)),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 75)) |> 
  hide_legend()


## comparative limited flux
limited_flux_data$Sim = "Altered"
initial_limited_flux$Sim = "Original"

compare_limited <- rbind(limited_flux_data, initial_limited_flux)

plot_ly(compare_limited) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~Sim, linetype = ~Sim) |> 
  layout(yaxis = list(#type = "log", exponentformat = "power",
    title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
    showline = TRUE, showgrid = FALSE, mirror = TRUE, range = c(0,0.0025)),
    xaxis = list(type = "log", exponentformat = "power", 
                 title_text = "Weight (g)", 
                 showline = TRUE, showgrid = FALSE, mirror = TRUE), 
    margin = list(b = 65, l = 75))


## growth rate plot
growth_data <- data.frame(growth = E_growth_reduced,
                          weight = w)


plot_ly(growth_data) |> 
  add_lines(x = ~weight, y = ~growth, ) |> 
  layout(yaxis = list(title_text = "Growth Rate (g/year)", 
                      range = c(0, 20),
                      showline = TRUE, showgrid = FALSE, zeroline = FALSE, mirror = TRUE),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 65))


## Comparative growth rate

gf_original <- melt(getEGrowth(params_single))
gf_original$Model <- "Original"
gf_starved <- melt(getEGrowth(params))
gf_starved$Model <- "Less prey"
gf_total <- rbind(gf_original, gf_starved)


plot_ly(gf_total) |> 
  add_lines(x = ~w, y = ~value, color = ~Model, linetype = ~Model) |> 
  layout(yaxis = list(title_text = "Growth rate (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE, zeroline = FALSE,
                      range = c(0,40)),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 65))


# Biomass density spectra
plotSpectra(sim_reduced, power = 2, wlim = c(1e-8, NA), ylim = c(1e-8, NA),
            time_range = 36)+
  theme_test()


# Comparative biomass density spectra
plotSpectra2(sim_initial, name1 = "Original",
             sim_reduced, name2 = "Less prey",
             power = 2, time_range = 36) + 
  theme_test()


# Animations --------------------------------------------------------------


# animation of change in biomass density over time
nf <- melt(sim_reduced@n)
n_ppf <- melt(sim_reduced@n_pp)
n_ppf$sp <- "Resource"
nf <- rbind(nf, n_ppf)

plot_ly(nf) |> 
  filter(w > 10^-5) |>  
  mutate(b = value * w^2) |> 
  add_lines(
    x = ~w, y = ~b,
    color = ~sp,
    frame = ~time,
    line = list(simplify = FALSE)) |> 
  layout(params, xaxis = list(type = "log", exponentformat = "power",
                         title_text = "body mass (g)"),
         yaxis = list(type = "log", exponentformat = "power",
                      title_text = "biomass (g/m^3)",
                      range = c(-9, 0)))



# animation of flux over time
species <- "Target species"
time <- getTimes(sim_reduced)

flux_series <- lapply(seq_along(time), function(i) {
  n_at_t <- sim_reduced@n[i, species, ]
  data.frame(
    Time = time[i],
    Weight = w,
    Flux = n_at_t * grr)}) |>  bind_rows()

flux_series$Species <- species

plot_ly(flux_series) |> 
  add_lines(x ~Weight, 
            y ~Flux, 
            color = ~Species, 
            frame = ~Time, 
            line = list(simplify = FALSE)) |>  
  layout(params, xaxis = list(type = "log", exponentformat = "power",
                         title_text = "body mass (g)"),
         yaxis = list(type = "log", exponentformat = "power",
                      title_text = "flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      range = c(-9, 0)))




# Altered fishing effort --------------------------------------------------

## flux with lower fishing effort than MSY
sim_lower <- project(params, t_max = 34, t_save = 0.1, effort = 0.5)

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


## flux with higher fishing effort
sim_higher <- project(params, t_max = 34, t_save = 0.1, effort = 10)

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
         margin = list(b = 65, l = 75))


## compare plots on same axis
reduced_flux_data$Sim <- NULL
reduced_flux_data$Effort = "MSY"
lower_flux_data$Effort = "Lower fishing effort"
higher_flux_data$Effort = "Higher fishing effort"


fishing_flux <- rbind(higher_flux_data, lower_flux_data, reduced_flux_data)

plot_ly(fishing_flux) |> 
  add_lines(x = ~Weight, y = ~Flux, color = ~Effort, linetype = ~Effort) |> 
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "Flux (gm<sup>-2</sup>year<sup>-1</sup>)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE, range = c(-12,-2)),
         xaxis = list(type = "log", exponentformat = "power", 
                      title_text = "Weight (g)", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE), 
         margin = list(b = 65, l = 75)) |> 
  hide_legend()


# Yield graphs ------------------------------------------------------------

## Yield graph for medium effort
sim_medium <- project(params, t_max = 100, t_save = 0.1, effort = 4.6)

yield <- getYield(sim_medium)
times <- getTimes(sim_medium)

yield_data <- data.frame(Time = times, 
                         Yield = yield)

plot_ly(yield_data) |> 
  filter(times >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE,
                        range = c(-7, -2)),
         xaxis = list(title = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 70))


# yield with lower fishing effort
sim_low_yield <- project(params, t_max = 100, t_save = 0.1, effort = 0.5)
yield_lower <- getYield(sim_low_yield)

yield_lower_data <- data.frame(Time = times, 
                               Yield = yield_lower)

plot_ly(yield_lower_data) |> 
  filter(times >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE,
                      range = c(-7, -2)),
         xaxis = list(title_text = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 70))



## Yield graph with high fishing effort
sim_high_yield <- project(params, t_max = 100, t_save = 0.1,  effort = 10)
yield_heavy <- getYield(sim_high_yield)

yield_heavy_data <- data.frame(Time = times, 
                         Yield = yield_heavy)

plot_ly(yield_heavy_data) |> 
  filter(times >= 50) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE,
                      range = c(-7, -2)),
         xaxis = list(title_text = "Year", 
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 70))



## time dependent fishing effort

t_total <- 200
med_effort <- 4.6
low_effort <- 0.5
threshold_value <- 0.00022
current_effort <- med_effort


sim_combined <- project(params, t_max = 1, t_save = 0.1, effort = current_effort)

last_yield <- getYield(sim_combined)[idxFinalT(sim_combined), "Target species"]
if (last_yield >= threshold_value) {
  current_effort <- low_effort
} else {
  current_effort <- med_effort
}

for (t in 2:t_total) {
  sim_combined <- project(sim_combined, t_max = 1, t_save = 0.1, effort = current_effort, append = TRUE)
  
  last_yield <- getYield(sim_combined)[idxFinalT(sim_combined), "Target species"]
  
  if (last_yield >= threshold_value) {
    current_effort <- low_effort
  } else {
    current_effort <- med_effort
  }
}

# Yield with changing effort
yield_fishing <- getYield(sim_combined)

time = getTimes(sim_combined)

yield_fishing_data <- data.frame(Time = time, 
                                 Yield = yield_fishing)

plot_ly(yield_fishing_data) |> 
  filter(time >= 150) |>
  add_lines(x = ~Time, y = ~Target.species) |> 
  layout(yaxis = list(type = "log", exponentformat = 'power', 
                      title_text = "Yield (g/year)",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE,
                      range = c(-7, -2)),
         xaxis = list(title_text = "Year",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 70))

# changing effort over time
effort_df <- melt(getEffort(sim_combined))

plot_ly(effort_df) |> 
  filter(time >= 100) |> 
  add_lines(x = ~time, y = ~value) |> 
  layout(yaxis = list(title_text = "Effort",
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         xaxis = list(title_text = "Time (year)", range = c(100, 200),
                      showline = TRUE, showgrid = FALSE, mirror = TRUE),
         margin = list(b = 65, l = 70))

