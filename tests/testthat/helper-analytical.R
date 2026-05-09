# Analytical reference solutions from Crank, "The Mathematics of Diffusion"
# (2nd ed., 1975). All time in minutes, lengths in um, diffusivity in um^2/min.
# Concentrations are in mg/um^3 internally; helpers below take and return mg
# per um^3 unless explicitly noted.
#
# These helpers exist only to give the test suite an absolute-truth target.
# They live in tests/testthat/helper-*.R so they are loaded by testthat for
# every test run but are not shipped in the package.

# ---- Crank section 4.3.3 (eq. 4.24a) ---------------------------------------
# Single homogeneous slab of thickness l, initial concentration zero.
# Donor side held at constant C_donor; receiver (perfect sink) held at zero.
# Returns Q(t), the cumulative mass per unit area that has crossed into the
# receiver, summed over n_terms of the Fourier series.
#
# Q(t) / (l * C_donor) =  D*t/l^2 - 1/6
#                       - (2/pi^2) * sum_{n=1..inf} ((-1)^n / n^2)
#                                    * exp(-D * n^2 * pi^2 * t / l^2)
crank_single_slab_Q <- function(t, l, D, C_donor, n_terms = 200) {
  vapply(t, function(ti) {
    if (ti <= 0) return(0)
    series <- 0
    for (n in seq_len(n_terms)) {
      exponent <- -D * n^2 * pi^2 * ti / l^2
      if (exponent < -700) break  # underflow guard
      term <- ((-1)^n) / n^2 * exp(exponent)
      series <- series + term
      if (abs(term) < 1e-18) break
    }
    l * C_donor * (D * ti / l^2 - 1 / 6 - (2 / pi^2) * series)
  }, numeric(1))
}

# Long-time linear asymptote of crank_single_slab_Q. Useful for checking the
# late-time slope and the lag time t_L = l^2 / (6*D) directly.
crank_single_slab_lag_time <- function(l, D) l^2 / (6 * D)
crank_single_slab_slope <- function(l, D, C_donor) C_donor * D / l

# Crank section 4.3.1 (eq. 4.16). Concentration profile c(x, t) inside a
# single homogeneous slab of thickness l with Dirichlet donor at x=0 (held
# at C_donor), perfect sink at x=l, and zero initial concentration:
#
#   c(x, t) = C_donor [1 - x/l - (2/π) Σ_{n=1}^∞ (1/n) sin(n π x/l) exp(-n²π²Dt/l²)]
#
# Returns a vector matching `x` at scalar `t`.
crank_single_slab_C <- function(x, t, l, D, C_donor, n_terms = 200) {
  if (t <= 0) {
    return(ifelse(x <= 0, C_donor, 0))
  }
  vapply(x, function(xi) {
    if (xi <= 0) return(C_donor)
    if (xi >= l) return(0)
    series <- 0
    for (n in seq_len(n_terms)) {
      exponent <- -n^2 * pi^2 * D * t / l^2
      if (exponent < -700) break
      term <- (1/n) * sin(n * pi * xi / l) * exp(exponent)
      series <- series + term
      if (abs(term) < 1e-20) break
    }
    C_donor * (1 - xi/l - (2/pi) * series)
  }, numeric(1))
}

# ---- Crank section 12 / Ash & Barrer ---------------------------------------
# Steady-state two-layer slab. Layer 1 of thickness l1 with diffusivity D1
# and partition K1, layer 2 of thickness l2 with D2, K2. Donor concentration
# C_donor at x=0 (reference frame: vehicle); receiver concentration zero at
# x = l1+l2. K_i is c_i / c_vehicle at equilibrium, so c at the top of layer
# 1 in steady state is K1 * C_donor.
#
# Returns the steady-state flux F (mass per unit area per unit time) and a
# function profile(x) giving the concentration at depth x in [0, l1+l2].
crank_two_layer_steady <- function(l1, D1, K1, l2, D2, K2, C_donor) {
  # Effective resistance R_i = l_i / (D_i * K_i / K_vehicle); with K_vehicle=1
  # the algebra below gives the steady-state flux F that satisfies
  #   K1 * C_donor = F * (l1 / D1 + (K1 / K2) * l2 / D2).
  R1 <- l1 / D1
  R2 <- (K1 / K2) * l2 / D2
  F  <- K1 * C_donor / (R1 + R2)

  c1_top <- K1 * C_donor                 # concentration at x = 0 (top of L1)
  c1_bot <- c1_top - F * l1 / D1         # concentration at x = l1 - eps
  c2_top <- c1_bot * K2 / K1             # partition jump at the interface
  c2_bot <- 0                            # perfect sink at x = l1 + l2

  profile <- function(x) {
    vapply(x, function(xi) {
      if (xi < 0 || xi > l1 + l2) return(NA_real_)
      if (xi <= l1) c1_top - (F / D1) * xi
      else          c2_top - (F / D2) * (xi - l1)
    }, numeric(1))
  }

  list(flux = F, profile = profile,
       c1_top = c1_top, c1_bot = c1_bot, c2_top = c2_top, c2_bot = c2_bot)
}

# ---- Kasting 2001 (J Pharm Sci 90:202) -------------------------------------
# Finite-dose / well-mixed-vehicle / homogeneous-membrane / perfect-sink BVP.
# Vehicle thickness h_v with concentration C_v0(t) (well-mixed: D_v -> infty),
# membrane thickness h with diffusivity D and partition K_mv across the
# vehicle/membrane interface, sink at x = h held at C = 0.
#
# β = K_mv · h / h_v is the dimensionless dose parameter.
#
# Eigenvalues α_n satisfy  α tan α = β  (eq. 9). The n-th root lies in
# ((n-1) π, (n-1/2) π) where the LHS sweeps continuously from 0 to +∞.

kasting_eigenvalues <- function(beta, n_terms = 100) {
  if (beta <= 0) stop("beta must be positive")
  alphas <- numeric(n_terms)
  eps <- 1e-10
  for (n in seq_len(n_terms)) {
    lo <- (n - 1) * pi + eps
    hi <- (n - 0.5) * pi - eps
    alphas[n] <- stats::uniroot(
      function(a) a * tan(a) - beta,
      lower = lo, upper = hi,
      tol = 1e-12
    )$root
  }
  alphas
}

# Cumulative permeated fraction M_t / M_∞ at the receptor (Kasting eq. 10):
#
#   M_t / M_∞ = 1 - 2 β Σ_n  exp(-α_n² D t / h²) / (cos α_n · (β + β² + α_n²))
#
# Returns a vector matching `t`. Hits 0 at t=0 and approaches 1 at t -> ∞.
kasting_M_over_Minf <- function(t, beta, D, h, n_terms = 100) {
  alphas <- kasting_eigenvalues(beta, n_terms)
  denom <- cos(alphas) * (beta + beta^2 + alphas^2)
  vapply(t, function(ti) {
    if (ti <= 0) return(0)
    exponents <- -alphas^2 * D * ti / h^2
    keep <- exponents > -700  # underflow guard
    if (!any(keep)) return(1)
    series <- sum(exp(exponents[keep]) / denom[keep])
    1 - 2 * beta * series
  }, numeric(1))
}

# Membrane concentration profile c(x, t) for the same Kasting BVP. Derived
# from the eigenfunction expansion implicit in eq. 6 by matching
# J(h, t) = -D ∂c/∂x|_{x=h} to the series in eq. 6:
#
#   c(x, t) = 2 K_mv C_v0 Σ_n  (α_n / (cos α_n · (β + β² + α_n²)))
#                              · sin(α_n (h - x) / h)
#                              · exp(-α_n² D t / h²)
#
# Sanity checks:
#   * at x = h:    sin(0) = 0          -> c(h, t) = 0   (sink BC)
#   * at x = 0:    sin(α_n) cos⁻¹(α_n) = tan(α_n) = β/α_n
#                  → c(0, t) = 2 β K_mv C_v0 Σ exp/(β+β²+α_n²)
#                            = K_mv · C_v(t)             (vehicle equilibrium)
#                  and C_v(0)/C_v0 = 2β Σ 1/(β+β²+α_n²) = 1.
#
# The series is conditionally convergent at t = 0 (the IC is discontinuous
# at x=0); the formula is intended for t > 0.
kasting_membrane_C <- function(x, t, beta, D, h, K_mv, C_v0, n_terms = 100) {
  if (t <= 0) {
    return(ifelse(x <= 0, K_mv * C_v0, 0))
  }
  alphas <- kasting_eigenvalues(beta, n_terms)
  coefs  <- alphas / (cos(alphas) * (beta + beta^2 + alphas^2))
  exps   <- exp(-alphas^2 * D * t / h^2)
  vapply(x, function(xi) {
    if (xi >= h) return(0)
    sins <- sin(alphas * (h - xi) / h)
    2 * K_mv * C_v0 * sum(coefs * sins * exps)
  }, numeric(1))
}


# ---- Crank section 3.3 (eq. 3.13) ------------------------------------------
# Semi-infinite medium, initial concentration zero, surface (x=0) held at C0.
# c(x, t) = C0 * erfc(x / (2 * sqrt(D*t)))
# Valid while the diffusion front has not yet reached the far boundary, i.e.
# while sqrt(D*t) << L (the actual finite domain length).
crank_semi_infinite <- function(x, t, D, C0) {
  if (t <= 0) return(rep(0, length(x)))
  C0 * erfc(x / (2 * sqrt(D * t)))
}

erfc <- function(z) 2 * pnorm(-z * sqrt(2))


# ---- helpers for unit conversion when comparing R-side mass output ---------

# Volume of a slab compartment in ml, given app_area in cm^2 and height in um.
slab_volume_ml <- function(app_area_cm2, height_um) {
  # 1 cm^2 = 1e8 um^2, 1 ml = 1e12 um^3
  app_area_cm2 * height_um * 1e-4
}

# Convert a per-unit-area mass quantity Q (in mg/um^2) into a total mass in
# the chosen scaling unit, given app_area in cm^2 and a scaling string.
mg_per_um2_to_total <- function(Q, app_area_cm2, scaling = c("mg", "ug", "ng")) {
  scaling <- match.arg(scaling)
  scale <- switch(scaling, mg = 1, ug = 1e3, ng = 1e6)
  Q * (app_area_cm2 * 1e8) * scale
}
