################################################################################
# SIMULATION
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)
library(splines)
library(runjags)
library(coda)
library(MASS)
library(ggplot2)
library(gridExtra)

# Set seed
set.seed(123)

# ==============================================================================
# 1. SIMULATION SETTINGS
# ==============================================================================

# Panel dimensions
N   <- 150
T   <- 8
p_x <- 3
Q   <- 1

# True scalar parameters 
true_beta0 <- -0.90
true_tau   <- 6.0
true_phi   <- 9.0
true_xi    <- 1.5
true_rho   <- 0.002
true_psi   <- 2.5

# True smooth frontier functions: s_j(x_j)
true_s1_fun <- function(x) 0.7 * (x + 0.20 * x^2)
true_s2_fun <- function(x) 0.6 * x - 1.0 * x^2 - 0.10 * x^3
true_s3_fun <- function(x) 0.5 * sin(pi * x)

# True smooth inefficiency function: s_1(z_1)
true_sz_fun <- function(z) -0.7 * z - 0.25 * z^2

# ==============================================================================
# 2. FIRM LOCATIONS AND DISTANCES IN KM
# ==============================================================================

# Longitude and latitude
lon <- runif(N, min = 5.0, max = 7.0)
lat <- runif(N, min = 50.0, max = 53.0)

# Distance matrix
coords <- cbind(lon, lat)
dist_matrix <- rdist.earth(coords, miles = FALSE)

# Spatial covariance implied by exponential kernel
cov_spatial <- exp(-true_rho * dist_matrix)
diag(cov_spatial) <- diag(cov_spatial) + 1e-6

# Spatial effect f(s_i)
Sigma_f <- (1 / true_psi) * cov_spatial
f_true <- as.numeric(mvrnorm(1, mu = rep(0, N), Sigma = Sigma_f))
f_true <- f_true - mean(f_true)

# ==============================================================================
# 3a. SIMULATE COVARIATES
# ==============================================================================

x1 <- matrix(runif(T * N, -1, 1), nrow = T, ncol = N)
x2 <- matrix(runif(T * N, -1, 1), nrow = T, ncol = N)
x3 <- matrix(runif(T * N, -1, 1), nrow = T, ncol = N)

z1 <- matrix(runif(T * N, -1.5, 1.5), nrow = T, ncol = N)

# ==============================================================================
# 3b. CENTER TRUE DGP COMPONENTS
# ==============================================================================

# True frontier smooths s_j(x_j)
s1_true <- true_s1_fun(x1)
s2_true <- true_s2_fun(x2)
s3_true <- true_s3_fun(x3)

s1_true <- s1_true - mean(s1_true)
s2_true <- s2_true - mean(s2_true)
s3_true <- s3_true - mean(s3_true)

# True inefficiency smooth s_1(z_1)
sz_true <- true_sz_fun(z1)
sz_true <- sz_true - mean(sz_true)

# ==============================================================================
# 4. SIMULATE RANDOM EFFECTS, INEFFICIENCY, OUTPUT
# ==============================================================================

# Firm-specific heterogeneity omega_i
omega_true <- rnorm(N, mean = 0, sd = sqrt(1 / true_xi))

y      <- matrix(NA, nrow = T, ncol = N)
u      <- matrix(NA, nrow = T, ncol = N)
alpha  <- matrix(NA, nrow = T, ncol = N)
u_star <- matrix(NA, nrow = T, ncol = N)

for (i in 1:N) {
  for (t in 1:T) {
    
    # Frontier mean: beta0 + sum_j s_j(x_jit) + omega_i
    mu_frontier <-
      true_beta0 +
      omega_true[i] +
      s1_true[t, i] +
      s2_true[t, i] +
      s3_true[t, i]
    
    # Scaling factor alpha_it = exp(s_1(z_1it) + f(s_i))
    alpha[t, i] <- exp(f_true[i] + sz_true[t, i])
    
    # Baseline inefficiency u_it^*
    u_star[t, i] <- abs(rnorm(1, mean = 0, sd = sqrt(1 / true_phi)))
    u[t, i] <- alpha[t, i] * u_star[t, i]
    
    # Noise term v_it
    v_it <- rnorm(1, mean = 0, sd = sqrt(1 / true_tau))
    
    # Output
    y[t, i] <- mu_frontier - u[t, i] + v_it
  }
}

# ==============================================================================
# 5. LONG DATA VECTORS FOR SPLINES
# ==============================================================================

x1_long <- as.vector(x1)
x2_long <- as.vector(x2)
x3_long <- as.vector(x3)
z1_long <- as.vector(z1)

# ==============================================================================
# 6. HELPER: CENTERED B-SPLINE BASIS
# ==============================================================================

make_B_array_specs <- function(x_long, T, N, nb, deg) {
  B0 <- bs(x_long, df = nb, degree = deg, intercept = TRUE)
  
  specs <- list(
    knots          = attr(B0, "knots"),
    Boundary.knots = attr(B0, "Boundary.knots"),
    degree         = attr(B0, "degree"),
    intercept      = TRUE,
    center         = colMeans(B0)
  )
  
  Bc <- sweep(B0, 2, specs$center, "-")
  list(B = array(Bc, dim = c(T, N, nb)), specs = specs)
}

# ==============================================================================
# 7. SPLINE BASIS FOR X-COVARIATES
# ==============================================================================

x_covars <- list(
  x1 = x1_long,
  x2 = x2_long,
  x3 = x3_long
)

nb_x  <- 6
deg_x <- 3

Bx_list <- vector("list", p_x)
spline_specs_x <- vector("list", p_x)

for (j in seq_len(p_x)) {
  tmp <- make_B_array_specs(x_covars[[j]], T, N, nb_x, deg_x)
  Bx_list[[j]] <- tmp$B
  spline_specs_x[[j]] <- tmp$specs
  spline_specs_x[[j]]$name <- names(x_covars)[j]
}

Bx <- array(0, dim = c(T, N, p_x * nb_x))
for (j in seq_len(p_x)) {
  idx <- ((j - 1) * nb_x + 1):(j * nb_x)
  Bx[, , idx] <- Bx_list[[j]]
}

Kx <- p_x * nb_x
idx_start_x <- (seq_len(p_x) - 1) * nb_x + 1
idx_end_x   <- seq_len(p_x) * nb_x

# ==============================================================================
# 8. SPLINE BASIS FOR Z-COVARIATE
# ==============================================================================

nb_z  <- 6
deg_z <- 3

tmpz <- make_B_array_specs(z1_long, T, N, nb_z, deg_z)
Bz <- tmpz$B
spline_specs_z <- tmpz$specs
spline_specs_z$name <- "z1"

# ==============================================================================
# 9. JAGS MODEL
# ==============================================================================

model_string <- "
model {

  # Firm-specific heterogeneity: omega_i
  for (i in 1:N) {
    omega[i] ~ dnorm(0, xi)
  }

  for (i in 1:N) {
    for (t in 1:T) {

      # Baseline inefficiency u_it^*
      u_star[t,i] ~ dnorm(0, phi) T(0, )

      # Scaling: alpha_it = exp(f(s_i) + s_1(z_1it))
      spl_z[t,i] <- inprod(beta_z[], Bz[t,i,])
      alpha[t,i] <- exp(f_spatial[i] + spl_z[t,i])

      u[t,i] <- alpha[t,i] * u_star[t,i]

      # Frontier smooths
      spl_x[t,i] <- inprod(beta_x[], Bx[t,i,])

      # Output equation
      meanrhsy[t,i] <- beta0 + omega[i] + spl_x[t,i] - u[t,i]
      y[t,i] ~ dnorm(meanrhsy[t,i], tau)
    }
  }

  # Spatial effect: f(s_i)
  f_spatial[1:N] ~ dmnorm(rep(0, N), precision_matrix)

  for (i in 1:N) {
    for (j in 1:N) {
      cov[i,j] <- exp(-rho * d[i,j]) + ifelse(i == j, 1e-6, 0)
    }
  }
  precision_matrix <- psi * inverse(cov)

  # Priors
  beta0 ~ dnorm(0, 1.0E-6)

  tau ~ dgamma(atau, btau)
  phi ~ dgamma(7.0, 0.5)

  xi ~ dgamma(a_xi, b_xi)

  psi ~ dgamma(2, 2)
  rho ~ dunif(0.001, 0.2)

  # RW2 priors for frontier splines
  for (j in 1:p_x) {
    beta_x[idx_start_x[j]]     ~ dnorm(0, 1.0E-6)
    beta_x[idx_start_x[j] + 1] ~ dnorm(0, 1.0E-6)

    for (k in (idx_start_x[j] + 2):idx_end_x[j]) {
      beta_x[k] ~ dnorm(2*beta_x[k-1] - beta_x[k-2], tau_x[j])
    }
    tau_x[j] ~ dgamma(0.001, 0.001)
  }

  # RW2 prior for inefficiency spline
  beta_z[1] ~ dnorm(0, 1.0E-6)
  beta_z[2] ~ dnorm(0, 1.0E-6)

  for (k in 3:nb_z) {
    beta_z[k] ~ dnorm(2*beta_z[k-1] - beta_z[k-2], tau_z)
  }

  tau_z ~ dgamma(0.001, 0.001)
}
"

# ==============================================================================
# 10. DATA FOR JAGS
# ==============================================================================

atau <- 0.001
btau <- 0.001
a_xi <- 0.001
b_xi <- 0.001

data_jags <- list(
  N = N,
  T = T,
  y = y,
  
  Bx = Bx,
  p_x = p_x,
  nb_x = nb_x,
  idx_start_x = idx_start_x,
  idx_end_x = idx_end_x,
  
  Bz = Bz,
  nb_z = nb_z,
  
  d = dist_matrix,
  
  atau = atau,
  btau = btau,
  a_xi = a_xi,
  b_xi = b_xi
)

# ==============================================================================
# 11. INITIAL VALUES
# ==============================================================================

allInits <- function() {
  list(
    beta0 = rnorm(1, true_beta0, 0.1),
    
    tau = rgamma(1, shape = 2, rate = 2 / true_tau),
    phi = rgamma(1, shape = 2, rate = 2 / true_phi),
    
    xi = rgamma(1, shape = 2, rate = 2 / true_xi),
    omega = rnorm(N, 0, 0.3),
    
    psi = rgamma(1, shape = 2, rate = 2 / true_psi),
    rho = runif(1, 0.001, 0.01),
    
    beta_x = rnorm(Kx, 0, 0.1),
    tau_x  = rgamma(p_x, 1, 1),
    
    beta_z = rnorm(nb_z, 0, 0.1),
    tau_z  = rgamma(1, 1, 1)
  )
}

# ==============================================================================
# 12. PARAMETERS TO MONITOR
# ==============================================================================

parameters <- c(
  "beta0",
  "tau",
  "phi",
  "xi",
  "psi",
  "rho",
  "beta_x",
  "beta_z",
  "tau_x",
  "tau_z",
  "deviance"
)

# ==============================================================================
# 13. RUN JAGS
# ==============================================================================

burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin

fit_sim <- run.jags(
  model    = model_string,
  data     = data_jags,
  inits    = allInits,
  monitor  = parameters,
  n.chains = chains,
  burnin   = burnin,
  sample   = sample,
  thin     = thin,
  method   = "parallel"
)

# ==============================================================================
# 14. SUMMARY TABLE
# ==============================================================================

summ_jags <- function(fit, pars, digits = 3) {
  m_list <- as.mcmc.list(fit$mcmc)
  keep   <- intersect(pars, varnames(m_list))
  
  m <- m_list[, keep, drop = FALSE]
  s <- summary(m)
  rh <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  ess <- effectiveSize(m)
  
  out <- data.frame(
    mean = s$statistics[, "Mean"],
    sd   = s$statistics[, "SD"],
    lo95 = s$quantiles[, "2.5%"],
    hi95 = s$quantiles[, "97.5%"],
    Rhat = rh,
    ESS  = as.numeric(ess),
    MCSE = s$statistics[, "SD"] / sqrt(as.numeric(ess)),
    check.names = FALSE
  )
  
  f <- function(x) formatC(x, format = "f", digits = digits)
  out[] <- lapply(out, function(x) if (is.numeric(x)) f(x) else x)
  out
}

pars_report <- c(
  "beta0",
  "tau",
  "phi",
  "xi",
  "psi",
  "rho",
  paste0("tau_x[", 1:p_x, "]"),
  "tau_z",
  "deviance"
)

tab_sim <- summ_jags(fit_sim, pars_report, digits = 3)
print(tab_sim, row.names = TRUE)

# ------------------------------------------------------------------------------
# 15. MCMC matrix
# ------------------------------------------------------------------------------
m_list <- as.mcmc.list(fit_sim$mcmc)
M <- as.matrix(m_list)

# helper to extract indexed parameters
get_param_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

# ------------------------------------------------------------------------------
# 16. Plot function: estimated smooth with 95% band + true function
# ------------------------------------------------------------------------------
plot_true_vs_estimated_smooth <- function(b_draws, x_obs, specs, true_fun,
                                          main = "", ngrid = 200,
                                          legend_pos = "topleft") {
  
  xg <- seq(min(x_obs, na.rm = TRUE), max(x_obs, na.rm = TRUE), length.out = ngrid)
  
  Bg0 <- bs(
    xg,
    knots = specs$knots,
    Boundary.knots = specs$Boundary.knots,
    degree = specs$degree,
    intercept = specs$intercept
  )
  
  Bg <- sweep(Bg0, 2, specs$center, "-")
  
  f_draws <- b_draws %*% t(Bg)
  f_draws <- sweep(f_draws, 1, rowMeans(f_draws), "-")
  
  f_mean <- colMeans(f_draws)
  f_lo   <- apply(f_draws, 2, quantile, 0.025)
  f_hi   <- apply(f_draws, 2, quantile, 0.975)
  
  f_true <- true_fun(xg)
  f_true <- f_true - mean(f_true)
  
  plot(xg, f_mean,
       type = "l", lwd = 2,
       xlab = "", ylab = "",
       main = main,
       ylim = range(c(f_lo, f_hi, f_true)))
  
  polygon(c(xg, rev(xg)),
          c(f_lo, rev(f_hi)),
          col = rgb(0, 0, 1, 0.15),
          border = NA)
  
  lines(xg, f_mean, lwd = 2, col = "blue")
  lines(xg, f_true, lwd = 2, lty = 2, col = "red")
  abline(h = 0, lty = 3)
  
  legend(legend_pos,
         legend = c("Estimated", "95% band", "True"),
         lwd = c(2, NA, 2),
         lty = c(1, NA, 2),
         pch = c(NA, 15, NA),
         pt.cex = c(NA, 2, NA),
         col = c("blue", rgb(0, 0, 1, 0.15), "red"),
         bty = "n")
}

# ------------------------------------------------------------------------------
# 17. Extract beta_x and beta_z posterior draws
# ------------------------------------------------------------------------------
draw_beta_x <- get_param_block(M, "beta_x", 1:Kx)
draw_beta_z <- get_param_block(M, "beta_z", 1:nb_z)

beta_x_blocks <- vector("list", p_x)
for (j in 1:p_x) {
  idxj <- idx_start_x[j]:idx_end_x[j]
  beta_x_blocks[[j]] <- draw_beta_x[, idxj, drop = FALSE]
}

beta_z_block <- draw_beta_z

# ------------------------------------------------------------------------------
# 18. Plot true vs estimated smooths for x1, x2, x3, z1
# ------------------------------------------------------------------------------
pdf("smooths_all.pdf", width = 8, height = 8)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot_true_vs_estimated_smooth(
  b_draws  = beta_x_blocks[[1]],
  x_obs    = x1_long,
  specs    = spline_specs_x[[1]],
  true_fun = true_s1_fun,
  main     = expression(s[1](x[1]))
)

plot_true_vs_estimated_smooth(
  b_draws  = beta_x_blocks[[2]],
  x_obs    = x2_long,
  specs    = spline_specs_x[[2]],
  true_fun = true_s2_fun,
  main     = expression(s[2](x[2]))
)

plot_true_vs_estimated_smooth(
  b_draws  = beta_x_blocks[[3]],
  x_obs    = x3_long,
  specs    = spline_specs_x[[3]],
  true_fun = true_s3_fun,
  main     = expression(s[3](x[3]))
)

plot_true_vs_estimated_smooth(
  b_draws    = beta_z_block,
  x_obs      = z1_long,
  specs      = spline_specs_z,
  true_fun   = true_sz_fun,
  main       = expression(s[1](z[1])),
  legend_pos = "bottomleft"
)

dev.off()
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 19. Posterior density plots for scalar parameters
# ------------------------------------------------------------------------------
posterior_df <- as.data.frame(as.matrix(m_list))

true_values <- c(
  beta0 = true_beta0,
  tau   = true_tau,
  phi   = true_phi,
  xi    = true_xi,
  psi   = true_psi,
  rho   = true_rho
)

param_names <- names(true_values)

display_names <- c(
  beta0 = expression(beta[0]),
  tau   = expression(tau),
  phi   = expression(phi),
  xi    = expression(xi),
  psi   = expression(psi),
  rho   = expression(rho)
)

plot_posterior <- function(param_name) {
  
  posterior_samples <- posterior_df[[param_name]]
  true_value <- true_values[param_name]
  
  lower_limit <- min(posterior_samples, na.rm = TRUE)
  upper_limit <- max(posterior_samples, na.rm = TRUE)
  range_width <- upper_limit - lower_limit
  
  lower_limit <- lower_limit - 0.15 * range_width
  upper_limit <- upper_limit + 0.15 * range_width
  
  if (abs(true_value) > 0) {
    lower_limit <- min(lower_limit, true_value - 0.10 * abs(true_value))
    upper_limit <- max(upper_limit, true_value + 0.10 * abs(true_value))
  } else {
    lower_limit <- min(lower_limit, true_value - 0.1)
    upper_limit <- max(upper_limit, true_value + 0.1)
  }
  
  positive_params <- c("tau", "phi", "xi", "psi", "rho")
  if (param_name %in% positive_params) {
    lower_limit <- max(0, lower_limit)
  }
  
  ggplot(posterior_df, aes(x = .data[[param_name]])) +
    geom_density(fill = "blue", alpha = 0.3) +
    geom_vline(xintercept = true_value,
               color = "red",
               linetype = "dashed",
               linewidth = 1) +
    scale_x_continuous(limits = c(lower_limit, upper_limit)) +
    labs(
      x = display_names[param_name],
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_line(linewidth = 0.3, color = "gray80"),
      panel.grid.minor = element_blank()
    )
}

plots <- lapply(param_names, plot_posterior)

pdf("posterior_densities.pdf", width = 10, height = 10)
grid.arrange(grobs = plots, ncol = 2, nrow = 3)
dev.off()