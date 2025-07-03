# SEIP_functions.R

library(deSolve)

# ---- ODE Definitions ----

# Basic SEIP
SEIP_basic <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dS <- -beta * S * I
    dE <- beta * S * I - alpha * E
    dI <- alpha * E - gamma * I
    dP <- m * I - delta * P
    
    list(c(dS, dE, dI, dP))
  })
}

# SEIP with temperature-dependent beta
SEIP_temp_beta <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    temp_t <- mean_temp + amplitude * sin(2 * pi * t / 365 + phase_shift)
    beta_prob_t <- 1 / (1 + exp(-(b0 + b1 * temp_t)))
    
    effective_beta <- beta * beta_prob_t
    
    dS <- -effective_beta * S * I
    dE <- effective_beta * S * I - alpha * E
    dI <- alpha * E - gamma * I
    dP <- m * I - delta * P
    
    list(c(dS, dE, dI, dP))
  })
}

# SEIP with temperature-dependent delta
SEIP_temp_delta <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    temp_t <- mean_temp + amplitude * sin(2 * pi * t / 365 + phase_shift)
    delta_t <- exp(a_delta + b_delta * temp_t)
    
    dS <- -beta * S * I
    dE <- beta * S * I - alpha * E
    dI <- alpha * E - gamma * I
    dP <- m * I - delta_t * P
    
    list(c(dS, dE, dI, dP))
  })
}

# SEIP with both temperature-dependent parameters
SEIP_temp_both <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    temp_t <- mean_temp + amplitude * sin(2 * pi * t / 365 + phase_shift)
    delta_t <- exp(a_delta + b_delta * temp_t)
    beta_prob_t <- 1 / (1 + exp(-(b0 + b1 * temp_t)))
    
    effective_beta <- beta * beta_prob_t
    
    dS <- -effective_beta * S * I
    dE <- effective_beta * S * I - alpha * E
    dI <- alpha * E - gamma * I
    dP <- m * I - delta_t * P
    
    list(c(dS, dE, dI, dP))
  })
}

# ---- Solver Wrapper ----

run_SEIP <- function(ode_func, params, init, times) {
  out <- deSolve::ode(
    y = init,
    times = times,
    func = ode_func,
    parms = params
  )
  as.data.frame(out)
}

# ---- Optional Helpers ----

plot_SEIP_trajectories <- function(sim_data, title = "SEIP Trajectories") {
  library(ggplot2)
  df_long <- tidyr::pivot_longer(sim_data, -time, names_to = "Compartment", values_to = "Value")
  
  ggplot(df_long, aes(x = time, y = Value, color = Compartment)) +
    geom_line() +
    labs(title = title, x = "Time (days)", y = "Proportion or Count") +
    theme_minimal(base_family = "Arial") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title=element_text(size=16))
}
