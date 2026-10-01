# ---------------------------------------------------------------------------
# Hodgkin-Huxley model of the squid giant axon, implemented in R
#
# Four coupled first-order differential equations:
#   dV/dt = (I_ext - I_Na - I_K - I_L) / C_m      membrane potential
#   dm/dt = alpha_m (1 - m) - beta_m m            Na+ activation gate
#   dh/dt = alpha_h (1 - h) - beta_h h            Na+ inactivation gate
#   dn/dt = alpha_n (1 - n) - beta_n n            K+ activation gate
#
# with I_Na = g_Na m^3 h (V - E_Na), I_K = g_K n^4 (V - E_K), I_L = g_L (V - E_L)
#
# Units: mV, ms, uA/cm^2, mS/cm^2, uF/cm^2. Resting potential is about -65 mV
# (the modern convention; the 1952 paper measured voltage relative to rest).
#
# Reference: Hodgkin AL, Huxley AF. A quantitative description of membrane
# current and its application to conduction and excitation in nerve.
# J Physiol. 1952;117(4):500-544.
# ---------------------------------------------------------------------------

library(deSolve)

hh_params <- list(
  C_m  = 1,        # membrane capacitance (uF/cm^2)
  g_Na = 120,      # maximum sodium conductance (mS/cm^2)
  g_K  = 36,       # maximum potassium conductance (mS/cm^2)
  g_L  = 0.3,      # leak conductance (mS/cm^2)
  E_Na = 50,       # sodium reversal potential (mV)
  E_K  = -77,      # potassium reversal potential (mV)
  E_L  = -54.387   # leak reversal potential (mV)
)

# x / (exp(x / y) - 1), with its limit (y) handled at x = 0. The rate
# functions below contain this form and are undefined exactly at V = -55
# and V = -40 mV unless the limit is handled.
vtrap <- function(x, y) {
  ifelse(abs(x / y) < 1e-6, y * (1 - x / y / 2), x / (exp(x / y) - 1))
}

# Voltage-dependent opening and closing rates (1/ms) of each gate
alpha_n <- function(V) 0.01  * vtrap(-(V + 55), 10)
beta_n  <- function(V) 0.125 * exp(-(V + 65) / 80)
alpha_m <- function(V) 0.1   * vtrap(-(V + 40), 10)
beta_m  <- function(V) 4     * exp(-(V + 65) / 18)
alpha_h <- function(V) 0.07  * exp(-(V + 65) / 20)
beta_h  <- function(V) 1 / (1 + exp(-(V + 35) / 10))

# Steady-state value of a gate at a fixed voltage: alpha / (alpha + beta)
steady_state <- function(V) {
  c(
    m = alpha_m(V) / (alpha_m(V) + beta_m(V)),
    h = alpha_h(V) / (alpha_h(V) + beta_h(V)),
    n = alpha_n(V) / (alpha_n(V) + beta_n(V))
  )
}

# Stimulus current: a step of amplitude I0 (uA/cm^2) between t_on and t_off (ms)
stimulus_current <- function(t, I0, t_on, t_off) {
  ifelse(t >= t_on & t < t_off, I0, 0)
}

hh_derivatives <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    I_ext <- stimulus_current(t, I0, t_on, t_off)
    I_Na  <- g_Na * m^3 * h * (V - E_Na)
    I_K   <- g_K  * n^4     * (V - E_K)
    I_L   <- g_L            * (V - E_L)

    dV <- (I_ext - I_Na - I_K - I_L) / C_m
    dm <- alpha_m(V) * (1 - m) - beta_m(V) * m
    dh <- alpha_h(V) * (1 - h) - beta_h(V) * h
    dn <- alpha_n(V) * (1 - n) - beta_n(V) * n

    list(c(dV, dm, dh, dn), I_ext = I_ext, I_Na = I_Na, I_K = I_K, I_L = I_L)
  })
}

# Simulate the neuron's response to a current step.
#   I0     stimulus current amplitude (uA/cm^2)
#   t_on   time the stimulus starts (ms)
#   t_off  time the stimulus ends (ms)
#   t_end  end of the simulation (ms)
#   dt     output time step (ms)
simulate_hh <- function(I0, t_on = 10, t_off = 40, t_end = 60, dt = 0.01) {
  V0    <- -65
  start <- c(V = V0, steady_state(V0))
  times <- seq(0, t_end, by = dt)
  parms <- c(hh_params, I0 = I0, t_on = t_on, t_off = t_off)

  out <- ode(
    y = start, times = times, func = hh_derivatives, parms = parms,
    method = "lsoda", rtol = 1e-8, atol = 1e-8,
    tcrit = NULL, hmax = 0.05
  )
  as.data.frame(out)
}

# Count spikes as upward crossings of 0 mV
count_spikes <- function(sim) {
  sum(sim$V[-1] >= 0 & sim$V[-nrow(sim)] < 0)
}
