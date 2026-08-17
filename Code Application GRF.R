################################################################################
# APPLICATION (MODEL 1)
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)
library(splines)
library(runjags)
library(coda)
library(ggplot2)
library(ggtext)

# Import data
farm <- read.table("spatial_GRF.txt", header = TRUE)

# Define data dimensions
NT <- nrow(farm)
T  <- 8
N  <- NT / T
stopifnot(N == as.integer(N))

# Response variable
y <- as.numeric(farm$log_y)
dim(y) <- c(T, N)

# Create distance matrix
coords <- cbind(farm$longi, farm$lat)
dist_matrix <- rdist.earth(coords, miles = FALSE)

# ------------------------------------------------------------------------------
# Helper: centered B-spline basis + specs
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Frontier covariates: x_it
# ------------------------------------------------------------------------------
x_covars <- list(
  logK  = farm$log_K,
  logL  = farm$log_L,
  logA  = farm$log_A,
  logS  = farm$log_S,
  logI  = farm$log_I,
  logF  = farm$log_F,
  trend = farm$trend
)

p_x   <- length(x_covars)
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

# ------------------------------------------------------------------------------
# Inefficiency covariate: z_it
# ------------------------------------------------------------------------------
subsidies <- matrix(farm$subsidiesallha, nrow = T, ncol = N, byrow = TRUE)
subsidies_std <- (subsidies - mean(subsidies)) / sd(subsidies)
z_long <- as.vector(subsidies_std)

nb_z  <- 6
deg_z <- 3

tmpz <- make_B_array_specs(z_long, T, N, nb_z, deg_z)
Bz <- tmpz$B
spline_specs_z <- tmpz$specs
spline_specs_z$name <- "subsidies_std"

# ------------------------------------------------------------------------------
# JAGS model string
# ------------------------------------------------------------------------------
model_string <- "
model {

  # Firm-specific heterogeneity
  for (i in 1:N) {
    omega[i] ~ dnorm(0, xi)
  }

  for (i in 1:N) {
    for (t in 1:T) {

      # Half-normal inefficiency: u_it^*
      u_star[t,i] ~ dnorm(0, phi) T(0, )

      # Scaling: alpha_it = exp(f(s_i) + s_z(z_it))
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

  # GRF prior 
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

  # rho
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

# ------------------------------------------------------------------------------
# Data for JAGS
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Initial values
# ------------------------------------------------------------------------------
allInits <- function() {
  list(
    beta0 = rnorm(1, 0, 0.1),
    
    tau = 100,
    phi = 10,
    
    xi = 1,
    omega = rnorm(N, 0, 0.1),
    
    psi = rgamma(1, 2, 2),
    rho = runif(1, 0.001, 0.2),
    
    beta_x = rnorm(Kx, 0, 0.05),
    tau_x  = rgamma(p_x, 1, 1),
    
    beta_z = rnorm(nb_z, 0, 0.05),
    tau_z  = rgamma(1, 1, 1)
  )
}

# ------------------------------------------------------------------------------
# Parameters to monitor
# ------------------------------------------------------------------------------
parameters <- c(
  "beta0", "tau", "phi", "xi", "psi", "rho",
  "beta_x", "beta_z", "tau_x", "tau_z",
  "deviance", "f_spatial", "u"
)

# ------------------------------------------------------------------------------
# Run JAGS
# ------------------------------------------------------------------------------
burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin
nc <- min(chains, parallel::detectCores())
runjags.options(ncores = nc)

t0 <- Sys.time()

fit <- run.jags(
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

print(Sys.time() - t0)

# ------------------------------------------------------------------------------
# Summary report
# ------------------------------------------------------------------------------
summ_jags <- function(fit, pars, digits = 3) {
  m_list <- as.mcmc.list(fit$mcmc)

  keep <- intersect(pars, varnames(m_list))
  if (length(keep) == 0) stop("None of the requested parameters were found in fit$mcmc.")

  m <- m_list[, keep, drop = FALSE]

  s   <- summary(m)
  rh  <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
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
  "xi",
  "tau",
  "phi",
  "psi",
  "rho",
  paste0("tau_x[", 1:p_x, "]"),
  "tau_z",
  "deviance"
)

tab <- summ_jags(fit, pars_report, digits = 3)
print(tab, row.names = TRUE)

# ------------------------------------------------------------------------------
# Traceplots
# ------------------------------------------------------------------------------
m_list <- as.mcmc.list(fit$mcmc)

chain_cols <- c(
  rgb(0.2, 0.4, 0.8, 0.5),
  rgb(0.8, 0.3, 0.3, 0.5),
  rgb(0.3, 0.6, 0.4, 0.5),
  rgb(0.8, 0.6, 0.3, 0.5)
)

nice_names <- list(
  beta0 = expression(beta[0]),
  xi    = expression(xi),
  tau   = expression(tau),
  phi   = expression(phi),
  psi   = expression(psi),
  rho   = expression(rho),
  tau_z = expression(tau[z])
)

for (j in 1:p_x) {
  nice_names[[paste0("tau_x[", j, "]")]] <- bquote(tau[x*","*.(j)])
}

plot_trace_group_fixed <- function(m_list, params, nice_names, nrow, ncol) {
  
  keep <- intersect(params, varnames(m_list))
  M_all <- as.matrix(m_list)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(nrow, ncol),
      mar = c(3, 3, 2, 1),
      cex.main = 1.3,
      cex.axis = 1.1,
      cex.lab = 1.1)
  
  for (p in keep) {
    yvals <- M_all[, p]
    
    plot(NULL,
         xlim = c(1, niter(m_list[[1]])),
         ylim = range(yvals, na.rm = TRUE),
         xlab = "Iteration",
         ylab = "",
         main = nice_names[[p]])
    
    for (ch in seq_along(m_list)) {
      lines(m_list[[ch]][, p], col = chain_cols[ch])
    }
  }
  
  n_empty <- nrow * ncol - length(keep)
  if (n_empty > 0) {
    for (i in seq_len(n_empty)) plot.new()
  }
}

pars_structural <- c("beta0", "xi", "tau", "phi", "psi", "rho")

pdf("traceplots_structural.pdf", width = 10, height = 6)
plot_trace_group_fixed(m_list, pars_structural, nice_names, nrow = 2, ncol = 3)
dev.off()

pars_smoothing <- c(
  paste0("tau_x[", 1:p_x, "]"),
  "tau_z"
)

pdf("traceplots_smoothing.pdf", width = 10, height = 10)
plot_trace_group_fixed(m_list, pars_smoothing, nice_names, nrow = 3, ncol = 3)
dev.off()

# ------------------------------------------------------------------------------
# Smooth plots
# ------------------------------------------------------------------------------
pdf("all_splines.pdf", width = 10, height = 10)

m_list <- as.mcmc.list(fit$mcmc)
M <- as.matrix(m_list)

get_param_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

plot_smooth <- function(b_draws, x_raw, specs, main, ngrid = 200) {
  
  xg <- seq(min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE), length.out = ngrid)
  
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
  
  line_col <- rgb(0.2, 0.4, 0.8, 0.9)
  band_col <- rgb(0.2, 0.4, 0.8, 0.2)
  
  plot(xg, f_mean,
       type = "l", lwd = 2,
       col = line_col,
       xlab = main,
       ylab = "",
       main = "")
  
  polygon(c(xg, rev(xg)),
          c(f_lo, rev(f_hi)),
          border = NA,
          col = band_col)
  
  lines(xg, f_mean, lwd = 2, col = line_col)
  abline(h = 0, lty = 2)
}

x_labels <- c("log K", "log L", "log A", "log S", "log I", "log F", "t")

oldpar <- par(no.readonly = TRUE)
par(cex.lab = 1.6, cex.axis = 1.5, cex.main = 1.2)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (j in seq_len(p_x)) {
  idx <- idx_start_x[j]:idx_end_x[j]
  b_j <- get_param_block(M, "beta_x", idx)
  
  plot_smooth(
    b_draws = b_j,
    x_raw   = x_covars[[j]],
    specs   = spline_specs_x[[j]],
    main    = x_labels[j]
  )
}

b_z <- get_param_block(M, "beta_z", 1:nb_z)

plot_smooth(
  b_draws = b_z,
  x_raw   = z_long,
  specs   = spline_specs_z,
  main    = "subsidies"
)

n_empty <- 9 - (p_x + 1)
if (n_empty > 0) {
  for (i in seq_len(n_empty)) plot.new()
}

par(oldpar)
dev.off()

# ------------------------------------------------------------------------------
# Marginal effects
# ------------------------------------------------------------------------------
eval_centered_basis <- function(x, specs) {
  B0 <- bs(
    x,
    knots = specs$knots,
    Boundary.knots = specs$Boundary.knots,
    degree = specs$degree,
    intercept = specs$intercept
  )
  sweep(B0, 2, specs$center, "-")
}

get_derivative_draws <- function(b_draws, x0, specs, h = 1e-4) {
  
  lower <- specs$Boundary.knots[1]
  upper <- specs$Boundary.knots[2]
  
  x_plus  <- min(x0 + h, upper - 1e-8)
  x_minus <- max(x0 - h, lower + 1e-8)
  
  B_plus  <- eval_centered_basis(x_plus, specs)
  B_minus <- eval_centered_basis(x_minus, specs)
  
  f_plus  <- as.vector(b_draws %*% t(B_plus))
  f_minus <- as.vector(b_draws %*% t(B_minus))
  
  (f_plus - f_minus) / (x_plus - x_minus)
}

x_labels <- c("log K", "log L", "log A", "log S", "log I", "log F", "t")

effects_list <- list()

for (j in seq_len(p_x)) {
  
  idx <- idx_start_x[j]:idx_end_x[j]
  b_j <- get_param_block(M, "beta_x", idx)
  
  x0 <- mean(x_covars[[j]], na.rm = TRUE)
  d_draws <- get_derivative_draws(b_j, x0, spline_specs_x[[j]])
  
  effects_list[[length(effects_list) + 1]] <- data.frame(
    Variable = x_labels[j],
    Mean     = mean(d_draws),
    Q2.5     = quantile(d_draws, 0.025),
    Q97.5    = quantile(d_draws, 0.975),
    ProbPos  = mean(d_draws > 0),
    ProbNeg  = mean(d_draws < 0),
    row.names = NULL
  )
}

b_z <- get_param_block(M, "beta_z", 1:nb_z)

z0 <- mean(z_long, na.rm = TRUE)
d_draws_z <- get_derivative_draws(b_z, z0, spline_specs_z)

effects_list[[length(effects_list) + 1]] <- data.frame(
  Variable = "subsidies",
  Mean     = mean(d_draws_z),
  Q2.5     = quantile(d_draws_z, 0.025),
  Q97.5    = quantile(d_draws_z, 0.975),
  ProbPos  = mean(d_draws_z > 0),
  ProbNeg  = mean(d_draws_z < 0),
  row.names = NULL
)

marginal_effects_tab <- do.call(rbind, effects_list)

marginal_effects_tab$Mean    <- round(marginal_effects_tab$Mean, 4)
marginal_effects_tab$Q2.5    <- round(marginal_effects_tab$Q2.5, 4)
marginal_effects_tab$Q97.5   <- round(marginal_effects_tab$Q97.5, 4)
marginal_effects_tab$ProbPos <- round(marginal_effects_tab$ProbPos, 3)
marginal_effects_tab$ProbNeg <- round(marginal_effects_tab$ProbNeg, 3)

print(marginal_effects_tab, row.names = FALSE)

# ------------------------------------------------------------------------------
# Inefficiency
# ------------------------------------------------------------------------------
dim_y <- dim(y)
T <- dim_y[1]
N <- dim_y[2]

u_names <- paste0("u[", rep(1:T, times = N), ",", rep(1:N, each = T), "]")
u_names <- u_names[u_names %in% colnames(M)]

length(u_names)

u_means_vec <- colMeans(M[, u_names, drop = FALSE])
u_means <- matrix(u_means_vec, nrow = T, ncol = N, byrow = FALSE)
u_all <- as.vector(u_means)

print(dim(u_means))

u_summary <- data.frame(
  mean  = mean(u_all),
  sd    = sd(u_all),
  q2.5  = quantile(u_all, 0.025),
  q97.5 = quantile(u_all, 0.975)
)

print(u_summary)

TE <- exp(-u_all)
summary(TE)

u_M1 <- as.matrix(as.mcmc.list(fit$mcmc))[ , grep("^u\\[", colnames(as.matrix(as.mcmc.list(fit$mcmc)))) ]
u_mean_M1 <- colMeans(u_M1)
summary(u_mean_M1)

uM1_summary <- c(
  mean  = mean(u_mean_M1),
  sd    = sd(u_mean_M1),
  q2.5  = quantile(u_mean_M1, 0.025),
  q97.5 = quantile(u_mean_M1, 0.975)
)

print(uM1_summary)

u_means_vector <- as.vector(u_means)

pdf("inefficiency_density.pdf", width = 6, height = 4)

ggplot(data = data.frame(Inefficiency = u_means_vector), aes(x = Inefficiency)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.7, size = 1) +
  labs(
    x = expression("Inefficiency " * u[it]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

dev.off()

summary(u_means_vector)
quantile(u_means_vector, probs = c(0.025, 0.975))
sd(u_means_vector)

# ------------------------------------------------------------------------------
# Spatial decay plot
# ------------------------------------------------------------------------------
rho_samples <- as.numeric(as.matrix(m_list)[,"rho"])

spatial_correlation_with_ci <- function(distance, rho_samples) {
  
  corrs <- sapply(rho_samples, function(rho) exp(-rho * distance))
  
  mean_corr <- rowMeans(corrs)
  lower_ci  <- apply(corrs, 1, quantile, 0.025)
  upper_ci  <- apply(corrs, 1, quantile, 0.975)
  
  list(mean = mean_corr, lower = lower_ci, upper = upper_ci)
}

distances <- seq(0, max(dist_matrix), length.out = 200)
correlation_results <- spatial_correlation_with_ci(distances, rho_samples)

data <- data.frame(
  Distance = distances,
  Correlation = correlation_results$mean,
  LowerCI = correlation_results$lower,
  UpperCI = correlation_results$upper
)

pdf("spatial_decay.pdf", width = 6, height = 4)

ggplot(data, aes(x = Distance, y = Correlation)) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI),
              fill = "blue", alpha = 0.25) +
  geom_line(color = "red", linewidth = 1.2) +
  labs(
    x = "Distance (km)",
    y = expression(exp(-rho %.% d(s[i], s[j])))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.1)
  )

dev.off()

# ------------------------------------------------------------------------------
# Histogram of spatial effects
# ------------------------------------------------------------------------------
f_names <- paste0("f_spatial[", 1:N, "]")
f_spatial <- summary(m_list[, f_names, drop = FALSE])$statistics[, "Mean"]

par(mar = c(5, 5, 2, 2))

hist(
  f_spatial,
  breaks = 20,
  col = "#69b3a2",
  border = "white",
  main = "",
  xlab = expression(bold(f(s[i]))),
  ylab = expression(bold("Frequency")),
  xlim = range(f_spatial),
  las = 1,
  cex.lab = 1.2,
  cex.axis = 1.1
)
box(col = "gray", lwd = 2)

df <- data.frame(SpatialEffect = f_spatial)

pdf("spatial_effect_histogram.pdf", width = 6, height = 4)

ggplot(df, aes(x = SpatialEffect)) +
  geom_histogram(aes(y = ..density..),
                 fill = "#69b3a2", color = "white",
                 bins = 30, alpha = 0.8) +
  labs(
    title = "",
    x = expression(f(s[i])),
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    panel.grid.major = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.minor = element_blank()
  )

dev.off()

# ------------------------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------------------------
pdf("spatial_effect_map_normalized.pdf", width = 6, height = 5)

coords_df <- data.frame(
  latitude  = farm$lat,
  longitude = farm$longi,
  transient = f_spatial
)

coords_df$longitude_norm <- (coords_df$longitude - min(coords_df$longitude)) /
  (max(coords_df$longitude) - min(coords_df$longitude))

coords_df$latitude_norm <- (coords_df$latitude - min(coords_df$latitude)) /
  (max(coords_df$latitude) - min(coords_df$latitude))

ggplot(coords_df, aes(x = longitude_norm, y = latitude_norm)) +
  geom_point(aes(color = transient), size = 5, shape = 16) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "",
    x = "Normalized Longitude",
    y = "Normalized Latitude",
    color = expression(f(s[i]))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

dev.off()

# ------------------------------------------------------------------------------
# Marginal log likelihood
# ------------------------------------------------------------------------------
get_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

draw_beta0 <- M[, "beta0"]
draw_xi    <- M[, "xi"]
draw_tau   <- M[, "tau"]
draw_phi   <- M[, "phi"]
draw_psi   <- M[, "psi"]
draw_rho   <- M[, "rho"]

draw_tau_x <- get_block(M, "tau_x", 1:p_x)
draw_tau_z <- M[, "tau_z"]

draw_beta_x <- get_block(M, "beta_x", 1:Kx)
draw_beta_z <- get_block(M, "beta_z", 1:nb_z)

draw_dev <- M[, "deviance"]

H <- cbind(
  beta0 = draw_beta0,
  xi    = draw_xi,
  tau   = draw_tau,
  phi   = draw_phi,
  psi   = draw_psi,
  rho   = draw_rho,
  draw_tau_x,
  tau_z = draw_tau_z,
  draw_beta_x,
  draw_beta_z
)

k <- ncol(H)

row_star <- which.min(draw_dev)
theta_star <- H[row_star, ]

i <- 0
take <- function(n = 1){
  out <- theta_star[(i+1):(i+n)]
  i <<- i + n
  out
}

beta0_star <- take(1)
xi_star    <- take(1)
tau_star   <- take(1)
phi_star   <- take(1)
psi_star   <- take(1)
rho_star   <- take(1)

tau_x_star <- take(p_x)
tau_z_star <- take(1)

beta_x_star <- take(Kx)
beta_z_star <- take(nb_z)

logp_beta0 <- dnorm(beta0_star, 0, sqrt(1/1e-6), log = TRUE)

logp_tau <- dgamma(tau_star, atau, btau, log = TRUE)
logp_xi  <- dgamma(xi_star, a_xi, b_xi, log = TRUE)
logp_phi <- dgamma(phi_star, 7, 0.5, log = TRUE)
logp_psi <- dgamma(psi_star, 2, 2, log = TRUE)

logp_rho <- if (rho_star >= 0.001 && rho_star <= 0.2) {
  log(1 / (0.2 - 0.001))
} else {
  -Inf
}

logp_beta_x <- 0
logp_tau_x <- 0

for (j in 1:p_x) {
  
  idxj <- idx_start_x[j]:idx_end_x[j]
  bj <- beta_x_star[idxj]
  
  logp_beta_x <- logp_beta_x +
    dnorm(bj[1], 0, sqrt(1/1e-6), log = TRUE) +
    dnorm(bj[2], 0, sqrt(1/1e-6), log = TRUE)
  
  for (k2 in 3:length(bj)) {
    mu <- 2 * bj[k2 - 1] - bj[k2 - 2]
    logp_beta_x <- logp_beta_x +
      dnorm(bj[k2], mu, 1 / sqrt(tau_x_star[j]), log = TRUE)
  }
  
  logp_tau_x <- logp_tau_x +
    dgamma(tau_x_star[j], 0.001, 0.001, log = TRUE)
}

logp_beta_z <- 0

logp_beta_z <- logp_beta_z +
  dnorm(beta_z_star[1], 0, sqrt(1/1e-6), log = TRUE) +
  dnorm(beta_z_star[2], 0, sqrt(1/1e-6), log = TRUE)

for (k2 in 3:nb_z) {
  mu <- 2 * beta_z_star[k2 - 1] - beta_z_star[k2 - 2]
  logp_beta_z <- logp_beta_z +
    dnorm(beta_z_star[k2], mu, 1 / sqrt(tau_z_star), log = TRUE)
}

logp_tau_z <- dgamma(tau_z_star, 0.001, 0.001, log = TRUE)

log_prior <- logp_beta0 + logp_xi + logp_tau + logp_phi +
  logp_psi + logp_rho + logp_beta_x + logp_tau_x +
  logp_beta_z + logp_tau_z

cov_H <- cov(H) + diag(1e-6, k)
logDetH <- -(determinant(cov_H, logarithm = TRUE)$modulus[1])

logL <- -0.5 * draw_dev[row_star]

logML <- 0.5 * k * log(2 * pi) + 0.5 * logDetH + log_prior + logL
print(unname(logML))

################################################################################
# MODEL 1 With Mundlak
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)
library(splines)
library(runjags)
library(coda)

# =========================
# Data
# =========================
farm <- read.table("spatial_GRF.txt", header = TRUE)

NT <- nrow(farm)
T  <- 8
N  <- NT / T
stopifnot(N == as.integer(N))

y <- as.numeric(farm$log_y)
dim(y) <- c(T, N)

# =========================
# Distance matrix 
# =========================
coords <- cbind(farm$longi, farm$lat)
dist_matrix <- rdist.earth(coords, miles = FALSE)

# =========================
# Helper: centered B-spline basis
# =========================
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

# =========================
# Frontier covariates: x_it
# =========================
x_covars <- list(
  logK  = farm$log_K,
  logL  = farm$log_L,
  logA  = farm$log_A,
  logS  = farm$log_S,
  logI  = farm$log_I,
  logF  = farm$log_F,
  trend = farm$trend
)

p_x   <- length(x_covars)
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

# =========================
# Correlated random effects (Mundlak)
# =========================
x_cre_names <- c("logK", "logL", "logA", "logS", "logI", "logF")
p_cre <- length(x_cre_names)

Xbar <- matrix(NA, nrow = N, ncol = p_cre)

for (j in seq_len(p_cre)) {
  x_mat <- matrix(x_covars[[x_cre_names[j]]], nrow = T, ncol = N, byrow = TRUE)
  Xbar[, j] <- colMeans(x_mat)
}

colnames(Xbar) <- paste0("mean_", x_cre_names)

# =========================
# Inefficiency covariate
# =========================
subsidies <- matrix(farm$subsidiesallha, nrow = T, ncol = N, byrow = TRUE)
subsidies_std <- (subsidies - mean(subsidies)) / sd(subsidies)
z_long <- as.vector(subsidies_std)

nb_z  <- 6
deg_z <- 3

tmpz <- make_B_array_specs(z_long, T, N, nb_z, deg_z)
Bz <- tmpz$B
spline_specs_z <- tmpz$specs
spline_specs_z$name <- "subsidies_std"

# =========================
# JAGS model
# =========================
model_string <- "
model {

  for (i in 1:N) {

    # Firm-specific heterogeneity
    omega[i] ~ dnorm(0, xi)

    # Correlated random effects (Mundlak) in the frontier
    cre_frontier[i] <- inprod(gamma[], Xbar[i,])

    for (t in 1:T) {

      # Half-normal transient inefficiency: u_it^*
      u_star[t,i] ~ dnorm(0, phi) T(0, )

      # Scaling function: alpha_it = exp(f(s_i) + s_z(z_it))
      spl_z[t,i] <- inprod(beta_z[], Bz[t,i,])
      alpha[t,i] <- exp(f_spatial[i] + spl_z[t,i])

      u[t,i] <- alpha[t,i] * u_star[t,i]

      # Frontier spline
      spl_x[t,i] <- inprod(beta_x[], Bx[t,i,])

      # Output equation
      meanrhsy[t,i] <- beta0 + omega[i] + cre_frontier[i] + spl_x[t,i] - u[t,i]
      y[t,i] ~ dnorm(meanrhsy[t,i], tau)
    }
  }

  # -----------------------------
  # GRF prior
  # -----------------------------
  f_spatial[1:N] ~ dmnorm(rep(0, N), precision_matrix)

  for (i in 1:N) {
    for (j in 1:N) {
      cov[i,j] <- exp(-rho * d[i,j]) + ifelse(i == j, 1e-6, 0)
    }
  }

  precision_matrix <- psi * inverse(cov)

  # -----------------------------
  # Priors
  # -----------------------------
  beta0 ~ dnorm(0, 1.0E-6)

  tau ~ dgamma(atau, btau)
  phi ~ dgamma(7.0, 0.5)

  xi ~ dgamma(a_xi, b_xi)

  psi ~ dgamma(2, 2)
  rho ~ dunif(0.001, 0.2)

  # Mundlak coefficients
  for (j in 1:p_cre) {
    gamma[j] ~ dnorm(0, 1.0E-6)
  }

  # -----------------------------
  # RW2 priors for frontier splines
  # -----------------------------
  for (j in 1:p_x) {

    beta_x[idx_start_x[j]]     ~ dnorm(0, 1.0E-6)
    beta_x[idx_start_x[j] + 1] ~ dnorm(0, 1.0E-6)

    for (k in (idx_start_x[j] + 2):idx_end_x[j]) {
      beta_x[k] ~ dnorm(2 * beta_x[k-1] - beta_x[k-2], tau_x[j])
    }

    tau_x[j] ~ dgamma(0.001, 0.001)
  }

  # -----------------------------
  # RW2 prior for inefficiency spline
  # -----------------------------
  beta_z[1] ~ dnorm(0, 1.0E-6)
  beta_z[2] ~ dnorm(0, 1.0E-6)

  for (k in 3:nb_z) {
    beta_z[k] ~ dnorm(2 * beta_z[k-1] - beta_z[k-2], tau_z)
  }

  tau_z ~ dgamma(0.001, 0.001)
}
"

# =========================
# Data for JAGS
# =========================
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
  
  Xbar = Xbar,
  p_cre = p_cre,
  
  d = dist_matrix,
  
  atau = atau,
  btau = btau,
  a_xi = a_xi,
  b_xi = b_xi
)

# =========================
# Initial values
# =========================
allInits <- function() {
  list(
    beta0 = rnorm(1, 0, 0.1),
    
    tau = 100,
    phi = 10,
    
    xi = 1,
    omega = rnorm(N, 0, 0.1),
    
    psi = rgamma(1, 2, 2),
    rho = runif(1, 0.001, 0.2),
    
    gamma = rnorm(p_cre, 0, 0.05),
    
    beta_x = rnorm(Kx, 0, 0.05),
    tau_x  = rgamma(p_x, 1, 1),
    
    beta_z = rnorm(nb_z, 0, 0.05),
    tau_z  = rgamma(1, 1, 1)
  )
}

# ------------------------------------------------------------------------------
# Parameters to monitor
# ------------------------------------------------------------------------------
parameters <- c(
  "beta0", "tau", "phi", "xi", "psi", "rho",
  "beta_x", "beta_z", "tau_x", "tau_z", "gamma",
  "deviance", "f_spatial", "u"
)

# =========================
# Run JAGS
# =========================
burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin

nc <- min(chains, parallel::detectCores())
runjags.options(ncores = nc)

t0 <- Sys.time()

fit <- run.jags(
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

print(Sys.time() - t0)

# =========================
# Summary table
# =========================
summ_jags <- function(fit, pars, digits = 3) {
  m_list <- as.mcmc.list(fit$mcmc)
  keep   <- intersect(pars, varnames(m_list))
  
  if (length(keep) == 0) stop("None of the requested parameters were found in fit$mcmc.")
  
  m <- m_list[, keep, drop = FALSE]
  
  s   <- summary(m)
  rh  <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
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
  paste0("gamma[", 1:p_cre, "]"),
  "xi",
  "tau",
  "phi",
  "psi",
  "rho",
  paste0("tau_x[", 1:p_x, "]"),
  "tau_z",
  "deviance"
)

tab <- summ_jags(fit, pars_report, digits = 3)
print(tab, row.names = TRUE)

# ------------------------------------------------------------------------------
# Marginal log likelihood
# ------------------------------------------------------------------------------

# -------------------------------
# Stack chains
# -------------------------------
m_list <- as.mcmc.list(fit$mcmc)
M <- as.matrix(m_list)

# -------------------------------
# Helper to pull indexed params
# -------------------------------
get_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

# -------------------------------
# Extract posterior draws
# -------------------------------
draw_beta0 <- M[, "beta0"]
draw_gamma <- get_block(M, "gamma", 1:p_cre)
draw_xi    <- M[, "xi"]
draw_tau   <- M[, "tau"]
draw_phi   <- M[, "phi"]
draw_psi   <- M[, "psi"]
draw_rho   <- M[, "rho"]

draw_tau_x <- get_block(M, "tau_x", 1:p_x)
draw_tau_z <- M[, "tau_z"]

draw_beta_x <- get_block(M, "beta_x", 1:Kx)
draw_beta_z <- get_block(M, "beta_z", 1:nb_z)

draw_dev <- M[, "deviance"]

# -------------------------------
# Parameter block
# -------------------------------
H <- cbind(
  beta0 = draw_beta0,
  draw_gamma,
  xi    = draw_xi,
  tau   = draw_tau,
  phi   = draw_phi,
  psi   = draw_psi,
  rho   = draw_rho,
  draw_tau_x,
  tau_z = draw_tau_z,
  draw_beta_x,
  draw_beta_z
)

k <- ncol(H)

# -------------------------------
# Minimum-deviance draw
# -------------------------------
row_star <- which.min(draw_dev)
theta_star <- H[row_star, ]

# helper to unpack theta_star
i <- 0
take <- function(n = 1) {
  out <- theta_star[(i + 1):(i + n)]
  i <<- i + n
  out
}

beta0_star <- take(1)
gamma_star <- take(p_cre)
xi_star    <- take(1)
tau_star   <- take(1)
phi_star   <- take(1)
psi_star   <- take(1)
rho_star   <- take(1)

tau_x_star <- take(p_x)
tau_z_star <- take(1)

beta_x_star <- take(Kx)
beta_z_star <- take(nb_z)

# -------------------------------
# Log-priors at star point
# -------------------------------

# beta0 ~ dnorm(0, 1e-6) [precision]
logp_beta0 <- dnorm(beta0_star, mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)

# gamma[j] ~ dnorm(0, 1e-6)
logp_gamma <- sum(dnorm(gamma_star, mean = 0, sd = sqrt(1 / 1e-6), log = TRUE))

# tau ~ dgamma(atau, btau)
logp_tau <- dgamma(tau_star, shape = atau, rate = btau, log = TRUE)

# xi ~ dgamma(a_xi, b_xi)
logp_xi <- dgamma(xi_star, shape = a_xi, rate = b_xi, log = TRUE)

# phi ~ dgamma(7, 0.5)
logp_phi <- dgamma(phi_star, shape = 7, rate = 0.5, log = TRUE)

# psi ~ dgamma(2, 2)
logp_psi <- dgamma(psi_star, shape = 2, rate = 2, log = TRUE)

# rho ~ dunif(0.001, 0.2)
logp_rho <- if (rho_star >= 0.001 && rho_star <= 0.2) {
  log(1 / (0.2 - 0.001))
} else {
  -Inf
}

# -------------------------------
# RW2 priors for beta_x and tau_x
# -------------------------------
logp_beta_x <- 0
logp_tau_x  <- 0

for (j in 1:p_x) {
  idxj <- idx_start_x[j]:idx_end_x[j]
  bj   <- beta_x_star[idxj]
  
  # first two coefficients
  logp_beta_x <- logp_beta_x + dnorm(bj[1], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
  logp_beta_x <- logp_beta_x + dnorm(bj[2], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
  
  # RW2 prior for the rest
  for (kk in 3:length(bj)) {
    mu_k <- 2 * bj[kk - 1] - bj[kk - 2]
    logp_beta_x <- logp_beta_x + dnorm(bj[kk], mean = mu_k, sd = 1 / sqrt(tau_x_star[j]), log = TRUE)
  }
  
  # tau_x[j] ~ dgamma(0.001, 0.001)
  logp_tau_x <- logp_tau_x + dgamma(tau_x_star[j], shape = 0.001, rate = 0.001, log = TRUE)
}

# -------------------------------
# RW2 priors for beta_z and tau_z
# -------------------------------
logp_beta_z <- 0

# first two coefficients
logp_beta_z <- logp_beta_z + dnorm(beta_z_star[1], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
logp_beta_z <- logp_beta_z + dnorm(beta_z_star[2], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)

# RW2 prior for the rest
for (kk in 3:nb_z) {
  mu_k <- 2 * beta_z_star[kk - 1] - beta_z_star[kk - 2]
  logp_beta_z <- logp_beta_z + dnorm(beta_z_star[kk], mean = mu_k, sd = 1 / sqrt(tau_z_star), log = TRUE)
}

# tau_z ~ dgamma(0.001, 0.001)
logp_tau_z <- dgamma(tau_z_star, shape = 0.001, rate = 0.001, log = TRUE)

# -------------------------------
# Total log prior
# -------------------------------
log_prior <- logp_beta0 + logp_gamma + logp_xi + logp_tau +
  logp_phi + logp_psi + logp_rho + logp_beta_x + logp_tau_x +
  logp_beta_z + logp_tau_z

# -------------------------------
# Laplace covariance term
# -------------------------------
cov_H <- cov(H) + diag(1e-6, k)   # small jitter for stability
logDetH <- -(determinant(cov_H, logarithm = TRUE)$modulus[1])

# -------------------------------
# Log-likelihood using minimum deviance
# -------------------------------
logL <- -0.5 * draw_dev[row_star]

# -------------------------------
# Marginal log likelihood
# -------------------------------
logML <- 0.5 * k * log(2 * pi) + 0.5 * logDetH + log_prior + logL
print(unname(logML))

################################################################################
# MODEL 2
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)
library(splines)
library(runjags)
library(coda)

# =========================
# Data
# =========================
farm <- read.table("spatial_GRF.txt", header = TRUE)

NT <- nrow(farm)
T  <- 8
N  <- NT / T
stopifnot(N == as.integer(N))

y <- as.numeric(farm$log_y)
dim(y) <- c(T, N)

# =========================
# Helper: centered B-spline basis + specs
# =========================
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

# =========================
# X splines (frontier covariates)
# =========================
x_covars <- list(
  logK  = farm$log_K,
  logL  = farm$log_L,
  logA  = farm$log_A,
  logS  = farm$log_S,
  logI  = farm$log_I,
  logF  = farm$log_F,
  trend = farm$trend
)

p_x   <- length(x_covars)
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

# =========================
# Z spline
# =========================
subsidies <- matrix(farm$subsidiesallha, nrow = T, ncol = N, byrow = TRUE)
subsidies_std <- (subsidies - mean(subsidies)) / sd(subsidies)
z_long <- as.vector(subsidies_std)

nb_z  <- 6
deg_z <- 3

tmpz <- make_B_array_specs(z_long, T, N, nb_z, deg_z)
Bz <- tmpz$B
spline_specs_z <- tmpz$specs
spline_specs_z$name <- "subsidies_std"

# =========================
# JAGS model string
# =========================
model_string <- "
model {

  # Firm-specific random effects
  for (i in 1:N) {
    omega[i] ~ dnorm(0, xi)
  }

  for (i in 1:N) {
    for (t in 1:T) {

      # Half-normal inefficiency: u_it^*
      u_star[t,i] ~ dnorm(0, phi) T(0, )

      # Scaling without spatial field
      spl_z[t,i] <- inprod(beta_z[], Bz[t,i,])
      alpha[t,i] <- exp(spl_z[t,i])

      u[t,i] <- alpha[t,i] * u_star[t,i]

      # Frontier spline for x's
      spl_x[t,i] <- inprod(beta_x[], Bx[t,i,])

      # Output equation
      meanrhsy[t,i] <- beta0 + omega[i] + spl_x[t,i] - u[t,i]
      y[t,i] ~ dnorm(meanrhsy[t,i], tau)
    }
  }

  # -----------------------------
  # Priors
  # -----------------------------
  beta0 ~ dnorm(0, 1.0E-6)

  tau ~ dgamma(atau, btau)
  phi ~ dgamma(7.0, 0.5)

  xi ~ dgamma(a_xi, b_xi)

  # -----------------------------
  # RW2 priors for x-splines
  # -----------------------------
  for (j in 1:p_x) {
    beta_x[idx_start_x[j]]     ~ dnorm(0, 1.0E-6)
    beta_x[idx_start_x[j] + 1] ~ dnorm(0, 1.0E-6)

    for (k in (idx_start_x[j] + 2):idx_end_x[j]) {
      beta_x[k] ~ dnorm(2*beta_x[k-1] - beta_x[k-2], tau_x[j])
    }
    tau_x[j] ~ dgamma(0.001, 0.001)
  }

  # -----------------------------
  # RW2 prior for z-spline
  # -----------------------------
  beta_z[1] ~ dnorm(0, 1.0E-6)
  beta_z[2] ~ dnorm(0, 1.0E-6)

  for (k in 3:nb_z) {
    beta_z[k] ~ dnorm(2*beta_z[k-1] - beta_z[k-2], tau_z)
  }

  tau_z ~ dgamma(0.001, 0.001)
}
"

# =========================
# Data for JAGS
# =========================
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
  
  atau = atau,
  btau = btau,
  a_xi = a_xi,
  b_xi = b_xi
)

# =========================
# Initial values
# =========================
allInits <- function() {
  list(
    beta0 = rnorm(1, 0, 0.1),
    
    tau = 100,
    phi = 10,
    
    xi = 1,
    omega = rnorm(N, 0, 0.1),
    
    beta_x = rnorm(Kx, 0, 0.05),
    tau_x  = rgamma(p_x, 1, 1),
    
    beta_z = rnorm(nb_z, 0, 0.05),
    tau_z  = rgamma(1, 1, 1)
  )
}

# ------------------------------------------------------------------------------
# Parameters to monitor
# ------------------------------------------------------------------------------
parameters <- c(
  "beta0", "tau", "phi", "xi", "psi", "rho",
  "beta_x", "beta_z", "tau_x", "tau_z",
  "deviance", "u"

# =========================
# Run JAGS 
# =========================
burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin
nc <- min(chains, parallel::detectCores())
runjags.options(ncores = nc)

t0 <- Sys.time()

fitnospat <- run.jags(
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

print(Sys.time() - t0)

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------
summ_jags <- function(fitnospat, pars, digits = 3) {
  m_list <- as.mcmc.list(fitnospat$mcmc)
  keep   <- intersect(pars, varnames(m_list))
  
  if (length(keep) == 0) stop("None of the requested parameters were found in fit$mcmc.")
  
  m <- m_list[, keep, drop = FALSE]
  
  s   <- summary(m)
  rh  <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
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

pars_report_nospat <- c(
  "beta0",
  "xi",
  "tau",
  "phi",
  paste0("tau_x[", 1:p_x, "]"),
  "tau_z",
  "deviance"
)

tab_nospat <- summ_jags(fitnospat, pars_report_nospat, digits = 3)
print(tab_nospat, row.names = TRUE)

# ------------------------------------------------------------------------------
# Marginal log likelihood
# ------------------------------------------------------------------------------

# -------------------------------
# Stack chains
# -------------------------------
m_list <- as.mcmc.list(fitnospat$mcmc)
M <- as.matrix(m_list)

# -------------------------------
# Helper to pull indexed params
# -------------------------------
get_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

# -------------------------------
# Extract posterior draws
# -------------------------------
draw_beta0 <- M[, "beta0"]
draw_xi    <- M[, "xi"]
draw_tau   <- M[, "tau"]
draw_phi   <- M[, "phi"]

draw_tau_x <- get_block(M, "tau_x", 1:p_x)
draw_tau_z <- M[, "tau_z"]

draw_beta_x <- get_block(M, "beta_x", 1:Kx)
draw_beta_z <- get_block(M, "beta_z", 1:nb_z)

draw_dev <- M[, "deviance"]

# -------------------------------
# Star point
# -------------------------------
i_star <- which.min(draw_dev)

beta0_star <- draw_beta0[i_star]
xi_star    <- draw_xi[i_star]
tau_star   <- draw_tau[i_star]
phi_star   <- draw_phi[i_star]

tau_x_star <- draw_tau_x[i_star, ]
tau_z_star <- draw_tau_z[i_star]

beta_x_star <- draw_beta_x[i_star, ]
beta_z_star <- draw_beta_z[i_star, ]

# -------------------------------
# Log priors at star point
# -------------------------------

logp_beta0 <- dnorm(beta0_star, 0, sqrt(1/1e-6), log = TRUE)

logp_tau <- dgamma(tau_star, shape = atau, rate = btau, log = TRUE)

logp_xi <- dgamma(xi_star, shape = a_xi, rate = b_xi, log = TRUE)

logp_phi <- dgamma(phi_star, shape = 7, rate = 0.5, log = TRUE)

# -------------------------------
# RW2 priors beta_x
# -------------------------------
logp_beta_x <- 0
logp_tau_x  <- 0

for (j in 1:p_x) {
  
  idxj <- idx_start_x[j]:idx_end_x[j]
  bj   <- beta_x_star[idxj]
  
  logp_beta_x <- logp_beta_x +
    dnorm(bj[1], 0, sqrt(1/1e-6), log = TRUE) +
    dnorm(bj[2], 0, sqrt(1/1e-6), log = TRUE)
  
  for (k in 3:length(bj)) {
    
    mu_k <- 2 * bj[k - 1] - bj[k - 2]
    
    logp_beta_x <- logp_beta_x +
      dnorm(bj[k], mu_k, 1 / sqrt(tau_x_star[j]), log = TRUE)
  }
  
  logp_tau_x <- logp_tau_x +
    dgamma(tau_x_star[j], 0.001, 0.001, log = TRUE)
}

# -------------------------------
# RW2 priors beta_z
# -------------------------------
logp_beta_z <- 0

logp_beta_z <- logp_beta_z +
  dnorm(beta_z_star[1], 0, sqrt(1/1e-6), log = TRUE) +
  dnorm(beta_z_star[2], 0, sqrt(1/1e-6), log = TRUE)

for (k in 3:nb_z) {
  
  mu_k <- 2 * beta_z_star[k - 1] - beta_z_star[k - 2]
  
  logp_beta_z <- logp_beta_z +
    dnorm(beta_z_star[k], mu_k, 1 / sqrt(tau_z_star), log = TRUE)
}

logp_tau_z <- dgamma(tau_z_star, 0.001, 0.001, log = TRUE)

# -------------------------------
# Total log prior
# -------------------------------
log_prior <- logp_beta0 + logp_xi + logp_tau + logp_phi +
  logp_beta_x + logp_tau_x + logp_beta_z + logp_tau_z

# -------------------------------
# Laplace covariance term
# -------------------------------
H <- cbind(
  beta0 = draw_beta0,
  xi    = draw_xi,
  tau   = draw_tau,
  phi   = draw_phi,
  draw_tau_x,
  tau_z = draw_tau_z,
  draw_beta_x,
  draw_beta_z
)

k <- ncol(H)

cov_H <- cov(H) + diag(1e-6, k)

logDetH <- -(determinant(cov_H, logarithm = TRUE)$modulus[1])

# -------------------------------
# Log likelihood
# -------------------------------
dev_star <- draw_dev[i_star]
logL <- -0.5 * dev_star

# -------------------------------
# Marginal log likelihood
# -------------------------------
logML <- 0.5 * k * log(2 * pi) + 0.5 * logDetH + log_prior + logL
print(unname(logML))

################################################################################
# MODEL 3
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)   
library(splines) 
library(runjags)
library(coda)

# =========================
# Data
# =========================
farm <- read.table("spatial_GRF.txt", header = TRUE)

NT <- nrow(farm)
T  <- 8
N  <- NT / T
stopifnot(N == as.integer(N))

y <- as.numeric(farm$log_y)
dim(y) <- c(T, N)

# =========================
# Distance matrix
# =========================
coords <- cbind(farm$longi, farm$lat)
dist_matrix <- rdist.earth(coords, miles = FALSE)

# =========================
# Helper: centered B-spline basis + specs
# =========================
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

# =========================
# X splines (frontier covariates)
# =========================
x_covars <- list(
  logK  = farm$log_K,
  logL  = farm$log_L,
  logA  = farm$log_A,
  logS  = farm$log_S,
  logI  = farm$log_I,
  logF  = farm$log_F,
  trend = farm$trend
)

p_x   <- length(x_covars)
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

# =========================
# JAGS model string
# =========================
model_string <- "
model {

  # Firm-specific random effects
  for (i in 1:N) {
    omega[i] ~ dnorm(0, xi)
  }

  for (i in 1:N) {
    for (t in 1:T) {

      # Half-normal inefficiency
      u_star[t,i] ~ dnorm(0, phi) T(0, )

      # Scaling depends only on spatial effect
      alpha[t,i] <- exp(f_spatial[i])

      u[t,i] <- alpha[t,i] * u_star[t,i]

      # Frontier spline for x's
      spl_x[t,i] <- inprod(beta_x[], Bx[t,i,])

      # Output equation
      meanrhsy[t,i] <- beta0 + omega[i] + spl_x[t,i] - u[t,i]
      y[t,i] ~ dnorm(meanrhsy[t,i], tau)
    }
  }

  # -----------------------------
  # GRF prior
  # -----------------------------
  f_spatial[1:N] ~ dmnorm(rep(0, N), precision_matrix)

  for (i in 1:N) {
    for (j in 1:N) {
      cov[i,j] <- exp(-rho * d[i,j]) + ifelse(i == j, 1e-6, 0)
    }
  }
  precision_matrix <- psi * inverse(cov)

  # -----------------------------
  # Priors
  # -----------------------------
  beta0 ~ dnorm(0, 1.0E-6)

  tau ~ dgamma(atau, btau)
  phi ~ dgamma(7.0, 0.5)

  xi ~ dgamma(a_xi, b_xi)

  psi ~ dgamma(2, 2)

  # rho in 1/km
  rho ~ dunif(0.001, 0.2)

  # -----------------------------
  # RW2 priors for x-splines
  # -----------------------------
  for (j in 1:p_x) {
    beta_x[idx_start_x[j]]     ~ dnorm(0, 1.0E-6)
    beta_x[idx_start_x[j] + 1] ~ dnorm(0, 1.0E-6)

    for (k in (idx_start_x[j] + 2):idx_end_x[j]) {
      beta_x[k] ~ dnorm(2*beta_x[k-1] - beta_x[k-2], tau_x[j])
    }
    tau_x[j] ~ dgamma(0.001, 0.001)
  }
}
"

# =========================
# Data for JAGS
# =========================
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
  
  d = dist_matrix,
  
  atau = atau,
  btau = btau,
  a_xi = a_xi,
  b_xi = b_xi
)

# =========================
# Initial values
# =========================
allInits <- function() {
  list(
    beta0 = rnorm(1, 0, 0.1),
    
    tau = 100,
    phi = 10,
    
    xi = 1,
    omega = rnorm(N, 0, 0.1),
    
    psi = rgamma(1, 2, 2),
    rho = runif(1, 0.001, 0.2),
    
    beta_x = rnorm(Kx, 0, 0.05),
    tau_x  = rgamma(p_x, 1, 1)
  )
}

# =========================
# Parameters to monitor
# =========================

parameters <- c(
  "beta0",
  "xi",
  "tau",
  "phi",
  "psi",
  "rho",
  "beta_x",
  "tau_x",
  "deviance",
  "f_spatial",
  "u"	
)

# =========================
# Run JAGS 
# =========================
burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin

t0 <- Sys.time()

fitonlysi <- run.jags(
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

print(Sys.time() - t0)

------------------------------------------------------------------------------
# SUMMARY TABLE
------------------------------------------------------------------------------

summ_jags <- function(fit, pars, digits = 3) {
  m_list <- as.mcmc.list(fit$mcmc)
  keep   <- intersect(pars, varnames(m_list))
  
  if (length(keep) == 0) stop("None of the requested parameters were found in fit$mcmc.")
  
  m <- m_list[, keep, drop = FALSE]
  
  s   <- summary(m)
  rh  <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
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

pars_report_onlysi <- c(
  "beta0",
  "xi",
  "tau",
  "phi",
  "psi",
  "rho",
  paste0("tau_x[", 1:p_x, "]"),
  "deviance"
)

tab_onlysi <- summ_jags(fitonlysi, pars_report_onlysi, digits = 3)
print(tab_onlysi, row.names = TRUE)

# ------------------------------------------------------------------------------
# Marginal log likelihood
# ------------------------------------------------------------------------------

m_list <- as.mcmc.list(fitonlysi$mcmc)
M <- as.matrix(m_list)

get_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

# Extract posterior draws
draw_beta0 <- M[, "beta0"]
draw_xi    <- M[, "xi"]
draw_tau   <- M[, "tau"]
draw_phi   <- M[, "phi"]
draw_psi   <- M[, "psi"]
draw_rho   <- M[, "rho"]

draw_tau_x  <- get_block(M, "tau_x", 1:p_x)
draw_beta_x <- get_block(M, "beta_x", 1:Kx)

draw_dev <- M[, "deviance"]

# Parameter block
H <- cbind(
  beta0 = draw_beta0,
  xi    = draw_xi,
  tau   = draw_tau,
  phi   = draw_phi,
  psi   = draw_psi,
  rho   = draw_rho,
  draw_tau_x,
  draw_beta_x
)

k <- ncol(H)

# Minimum-deviance draw
row_star <- which.min(draw_dev)
theta_star <- H[row_star, ]

# unpack parameters
i <- 0
take <- function(n = 1) {
  out <- theta_star[(i + 1):(i + n)]
  i <<- i + n
  out
}

beta0_star <- take(1)
xi_star    <- take(1)
tau_star   <- take(1)
phi_star   <- take(1)
psi_star   <- take(1)
rho_star   <- take(1)

tau_x_star  <- take(p_x)
beta_x_star <- take(Kx)

# Log priors
logp_beta0 <- dnorm(beta0_star, mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
logp_tau   <- dgamma(tau_star, shape = atau, rate = btau, log = TRUE)
logp_xi    <- dgamma(xi_star, shape = a_xi, rate = b_xi, log = TRUE)
logp_phi   <- dgamma(phi_star, shape = 7, rate = 0.5, log = TRUE)
logp_psi   <- dgamma(psi_star, shape = 2, rate = 2, log = TRUE)

logp_rho <- if (rho_star >= 0.001 && rho_star <= 0.2) {
  log(1 / (0.2 - 0.001))
} else {
  -Inf
}

# RW2 priors for beta_x and tau_x
logp_beta_x <- 0
logp_tau_x  <- 0

for (j in 1:p_x) {
  idxj <- idx_start_x[j]:idx_end_x[j]
  bj   <- beta_x_star[idxj]
  
  logp_beta_x <- logp_beta_x +
    dnorm(bj[1], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE) +
    dnorm(bj[2], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
  
  for (kk in 3:length(bj)) {
    mu_k <- 2 * bj[kk - 1] - bj[kk - 2]
    logp_beta_x <- logp_beta_x +
      dnorm(bj[kk], mean = mu_k, sd = 1 / sqrt(tau_x_star[j]), log = TRUE)
  }
  
  logp_tau_x <- logp_tau_x +
    dgamma(tau_x_star[j], shape = 0.001, rate = 0.001, log = TRUE)
}

# Total log prior
log_prior <- logp_beta0 + logp_xi + logp_tau + logp_phi +
  logp_psi + logp_rho + logp_beta_x + logp_tau_x

# Laplace covariance term
cov_H <- cov(H) + diag(1e-6, k)
logDetH <- -(determinant(cov_H, logarithm = TRUE)$modulus[1])

# Log likelihood 
logL <- -0.5 * draw_dev[row_star]

# Marginal log likelihood
logML_onlysi <- 0.5 * k * log(2 * pi) + 0.5 * logDetH + log_prior + logL
print(unname(logML_onlysi))

################################################################################
# MODEL 4
################################################################################

# Clear memory
rm(list = ls(all = TRUE))

# Load packages
library(fields)
library(splines)
library(runjags)
library(coda)

# =========================
# Data
# =========================
farm <- read.table("spatial_GRF.txt", header = TRUE)

NT <- nrow(farm)
T  <- 8
N  <- NT / T
stopifnot(N == as.integer(N))

y <- as.numeric(farm$log_y)
dim(y) <- c(T, N)

# =========================
# Helper: centered B-spline basis + specs
# =========================
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

# =========================
# X splines (frontier covariates)
# =========================
x_covars <- list(
  logK  = farm$log_K,
  logL  = farm$log_L,
  logA  = farm$log_A,
  logS  = farm$log_S,
  logI  = farm$log_I,
  logF  = farm$log_F,
  trend = farm$trend
)

p_x   <- length(x_covars)
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

# =========================
# JAGS model string
# =========================
model_string <- "
model {

  # Firm-specific random effects
  for (i in 1:N) {
    omega[i] ~ dnorm(0, xi)
  }

  for (i in 1:N) {
    for (t in 1:T) {

      # Half-normal inefficiency without scaling
      u[t,i] ~ dnorm(0, phi) T(0, )

      # Frontier spline for x's
      spl_x[t,i] <- inprod(beta_x[], Bx[t,i,])

      # Output equation
      meanrhsy[t,i] <- beta0 + omega[i] + spl_x[t,i] - u[t,i]
      y[t,i] ~ dnorm(meanrhsy[t,i], tau)
    }
  }

  # -----------------------------
  # Priors
  # -----------------------------
  beta0 ~ dnorm(0, 1.0E-6)

  tau ~ dgamma(atau, btau)
  phi ~ dgamma(7.0, 0.5)

  xi ~ dgamma(a_xi, b_xi)

  # -----------------------------
  # RW2 priors for x-splines
  # -----------------------------
  for (j in 1:p_x) {
    beta_x[idx_start_x[j]]     ~ dnorm(0, 1.0E-6)
    beta_x[idx_start_x[j] + 1] ~ dnorm(0, 1.0E-6)

    for (k in (idx_start_x[j] + 2):idx_end_x[j]) {
      beta_x[k] ~ dnorm(2 * beta_x[k-1] - beta_x[k-2], tau_x[j])
    }

    tau_x[j] ~ dgamma(0.001, 0.001)
  }
}
"

# =========================
# Data for JAGS
# =========================
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
  
  atau = atau,
  btau = btau,
  a_xi = a_xi,
  b_xi = b_xi
)

# =========================
# Initial values
# =========================
allInits <- function() {
  list(
    beta0 = rnorm(1, 0, 0.1),
    
    tau = 100,
    phi = 10,
    
    xi = 1,
    omega = rnorm(N, 0, 0.1),
    
    beta_x = rnorm(Kx, 0, 0.05),
    tau_x  = rgamma(p_x, 1, 1)
  )
}

# =========================
# Parameters to monitor
# =========================

parameters <- c(
  "beta0",
  "xi",
  "tau",
  "phi",
  "beta_x",
  "tau_x",
  "deviance",
  "u"
)


# =========================
# Run JAGS
# =========================
burnin <- 30000
iters  <- 100000
thin   <- 1
chains <- 4

sample <- iters - burnin
nc <- min(chains, parallel::detectCores())
runjags.options(ncores = nc)

t0 <- Sys.time()

fitnospatscal <- run.jags(
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

print(Sys.time() - t0)

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------

summ_jags <- function(fit, pars, digits = 3) {
  m_list <- as.mcmc.list(fit$mcmc)
  keep   <- intersect(pars, varnames(m_list))
  
  if (length(keep) == 0) stop("None of the requested parameters were found in fit$mcmc.")
  
  m <- m_list[, keep, drop = FALSE]
  
  s   <- summary(m)
  rh  <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
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

pars_report_nospatscal <- c(
  "beta0",
  "xi",
  "tau",
  "phi",
  paste0("tau_x[", 1:p_x, "]"),
  "deviance"
)

tab_nospatscal <- summ_jags(fitnospatscal, pars_report_nospatscal, digits = 3)
print(tab_nospatscal, row.names = TRUE)

# ------------------------------------------------------------------------------
# Marginal log likelihood
# ------------------------------------------------------------------------------

# -------------------------------
# Stack chains
# -------------------------------
m_list <- as.mcmc.list(fitnospatscal$mcmc)
M <- as.matrix(m_list)

# -------------------------------
# Helper to pull indexed JAGS params
# -------------------------------
get_block <- function(M, base, idx) {
  cols <- paste0(base, "[", idx, "]")
  cols <- cols[cols %in% colnames(M)]
  if (length(cols) == 0) stop(paste("No columns found for", base))
  M[, cols, drop = FALSE]
}

# -------------------------------
# Extract posterior draws
# -------------------------------
draw_beta0 <- M[, "beta0"]
draw_xi    <- M[, "xi"]
draw_tau   <- M[, "tau"]
draw_phi   <- M[, "phi"]

draw_tau_x  <- get_block(M, "tau_x", 1:p_x)
draw_beta_x <- get_block(M, "beta_x", 1:Kx)

draw_dev <- M[, "deviance"]

# -------------------------------
# Build parameter block
# -------------------------------
H <- cbind(
  beta0 = draw_beta0,
  xi    = draw_xi,
  tau   = draw_tau,
  phi   = draw_phi,
  draw_tau_x,
  draw_beta_x
)

k <- ncol(H)

# -------------------------------
# Minimum-deviance draw as star point
# -------------------------------
row_star <- which.min(draw_dev)
theta_star <- H[row_star, ]

# helper to unpack theta_star
i <- 0
take <- function(n = 1) {
  out <- theta_star[(i + 1):(i + n)]
  i <<- i + n
  out
}

beta0_star <- take(1)
xi_star    <- take(1)
tau_star   <- take(1)
phi_star   <- take(1)
tau_x_star <- take(p_x)
beta_x_star <- take(Kx)

# -------------------------------
# Log-priors at star point
# -------------------------------

# beta0 ~ dnorm(0, 1e-6) [precision]
logp_beta0 <- dnorm(beta0_star, mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)

# tau ~ dgamma(atau, btau)
logp_tau <- dgamma(tau_star, shape = atau, rate = btau, log = TRUE)

# xi ~ dgamma(a_xi, b_xi)
logp_xi <- dgamma(xi_star, shape = a_xi, rate = b_xi, log = TRUE)

# phi ~ dgamma(7, 0.5)
logp_phi <- dgamma(phi_star, shape = 7, rate = 0.5, log = TRUE)

# -------------------------------
# RW2 priors for beta_x and tau_x
# -------------------------------
logp_beta_x <- 0
logp_tau_x  <- 0

for (j in 1:p_x) {
  idxj <- idx_start_x[j]:idx_end_x[j]
  bj   <- beta_x_star[idxj]
  
  # first two coefficients
  logp_beta_x <- logp_beta_x + dnorm(bj[1], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
  logp_beta_x <- logp_beta_x + dnorm(bj[2], mean = 0, sd = sqrt(1 / 1e-6), log = TRUE)
  
  # RW2 prior for remaining coefficients
  for (kk in 3:length(bj)) {
    mu_k <- 2 * bj[kk - 1] - bj[kk - 2]
    logp_beta_x <- logp_beta_x + dnorm(bj[kk], mean = mu_k, sd = 1 / sqrt(tau_x_star[j]), log = TRUE)
  }
  
  # tau_x[j] ~ dgamma(0.001, 0.001)
  logp_tau_x <- logp_tau_x + dgamma(tau_x_star[j], shape = 0.001, rate = 0.001, log = TRUE)
}

# -------------------------------
# Total log prior
# -------------------------------
log_prior <- logp_beta0 + logp_xi + logp_tau + logp_phi +
  logp_beta_x + logp_tau_x

# -------------------------------
# Laplace covariance term
# -------------------------------
cov_H <- cov(H) + diag(1e-6, k)   # small jitter for numerical stability
logDetH <- -(determinant(cov_H, logarithm = TRUE)$modulus[1])

# -------------------------------
# Log-likelihood using MINIMUM deviance
# -------------------------------
logL <- -0.5 * draw_dev[row_star]

# -------------------------------
# Marginal log likelihood
# -------------------------------
logML <- 0.5 * k * log(2 * pi) + 0.5 * logDetH + log_prior + logL
print(unname(logML))

#######################################################################################
# INEFFICIENCY COMPARISON ACROSS MODELS (PROVIDED ALL RESULTS ARE SAVED POSTESTIMATION)
#######################################################################################

### COMPARISONS M1 VS M2rm(list=ls(all=TRUE))setwd("~/Documents")env1 <- new.env()env2 <- new.env()load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M1u.RData", envir = env1)load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M2u.RData", envir = env2)library(coda)get_fit <- function(env) {  nms <- ls(env)  for (nm in nms) {    x <- get(nm, envir = env)    if (is.list(x) && "mcmc" %in% names(x)) return(x)  }  stop("No fitted object found.")}# memory-safe posterior meansextract_u_mean <- function(fit_obj) {  mcmc_list <- fit_obj$mcmc  u_names <- grep("^u\\[", varnames(mcmc_list), value = TRUE)    sapply(u_names, function(v) {    mean(as.matrix(mcmc_list[, v]))  })}fit1 <- get_fit(env1)fit2 <- get_fit(env2)u_M1_mean <- extract_u_mean(fit1)u_M2_mean <- extract_u_mean(fit2)# plot directlyplot(u_M1_mean, u_M2_mean,     pch = 16, col = rgb(0,0,0,0.4),     xlab = "Inefficiency (M1)",     ylab = "Inefficiency (M2)")abline(0, 1, col = "red", lwd = 2)diff_u <- u_M2_mean - u_M1_meancol_vec <- colorRampPalette(c("blue", "white", "red"))(100)pdf("diff_M1_M2.pdf", width = 5, height = 5)idx <- cut(diff_u, breaks = 100, labels = FALSE)plot(u_M1_mean, diff_u,     pch = 16,     col = col_vec[idx],     ylim = c(-0.006, 0.006),     xlab = expression(paste("Inefficiency ", u[it], " (Model 1)")),     ylab = expression(paste("Difference ", u[it], ": Model 2 - Model 1")),     cex.lab = 1.2,   # axis titles bigger     cex.axis = 1.2   # tick labels bigger)abline(h = 0, col = "black", lwd = 2)dev.off()### COMPARISONS M1 VS M3rm(list=ls(all=TRUE))setwd("~/Documents")env1 <- new.env()env3 <- new.env()load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M1u.RData", envir = env1)load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M3u.RData", envir = env3)library(coda)get_fit <- function(env) {  nms <- ls(env)  for (nm in nms) {    x <- get(nm, envir = env)    if (is.list(x) && "mcmc" %in% names(x)) return(x)  }  stop("No fitted object found.")}# memory-safe posterior meansextract_u_mean <- function(fit_obj) {  mcmc_list <- fit_obj$mcmc  u_names <- grep("^u\\[", varnames(mcmc_list), value = TRUE)    sapply(u_names, function(v) {    mean(as.matrix(mcmc_list[, v]))  })}fit1 <- get_fit(env1)fit3 <- get_fit(env3)u_M1_mean <- extract_u_mean(fit1)u_M3_mean <- extract_u_mean(fit3)# plot directlyplot(u_M1_mean, u_M3_mean,     pch = 16, col = rgb(0,0,0,0.4),     xlab = "Inefficiency (M1)",     ylab = "Inefficiency (M3)")abline(0, 1, col = "red", lwd = 2)diff_u <- u_M3_mean - u_M1_meancol_vec <- colorRampPalette(c("blue", "white", "red"))(100)pdf("diff_M1_M3.pdf", width = 5, height = 5)# map colors to differencesidx <- cut(diff_u, breaks = 100, labels = FALSE)plot(u_M1_mean, diff_u,     pch = 16,     col = col_vec[idx],     xlab = expression(paste("Inefficiency ", u[it], " (Model 1)")),     ylab = expression(paste("Difference ", u[it], ": Model 3 - Model 1")),     cex.lab = 1.2,   # axis titles bigger     cex.axis = 1.2   # tick labels bigger)abline(h = 0, col = "black", lwd = 2)dev.off()### COMPARISONS M1 VS M4rm(list=ls(all=TRUE))setwd("~/Documents")env1 <- new.env()env4 <- new.env()load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M1u.RData", envir = env1)load("/Users/giovanni/Desktop/Research/Papers/GRF/revision/results/M4u.RData", envir = env4)library(coda)get_fit <- function(env) {  nms <- ls(env)  for (nm in nms) {    x <- get(nm, envir = env)    if (is.list(x) && "mcmc" %in% names(x)) return(x)  }  stop("No fitted object found.")}# memory-safe posterior meansextract_u_mean <- function(fit_obj) {  mcmc_list <- fit_obj$mcmc  u_names <- grep("^u\\[", varnames(mcmc_list), value = TRUE)    sapply(u_names, function(v) {    mean(as.matrix(mcmc_list[, v]))  })}fit1 <- get_fit(env1)fit4 <- get_fit(env4)u_M1_mean <- extract_u_mean(fit1)u_M4_mean <- extract_u_mean(fit4)# plot directlyplot(u_M1_mean, u_M4_mean,     pch = 16, col = rgb(0,0,0,0.4),     xlab = "Inefficiency (M1)",     ylab = "Inefficiency (M4)")abline(0, 1, col = "red", lwd = 2)diff_u <- u_M4_mean - u_M1_meancol_vec <- colorRampPalette(c("blue", "white", "red"))(100)pdf("diff_M1_M4.pdf", width = 5, height = 5)# map colors to differencesidx <- cut(diff_u, breaks = 100, labels = FALSE)plot(u_M1_mean, diff_u,     pch = 16,     col = col_vec[idx],     ylim = c(-0.006, 0.006),     xlab = expression(paste("Inefficiency ", u[it], " (Model 1)")),     ylab = expression(paste("Difference ", u[it], ": Model 4 - Model 1")),     cex.lab = 1.2,   # axis titles bigger     cex.axis = 1.2   # tick labels bigger)abline(h = 0, col = "black", lwd = 2)dev.off()