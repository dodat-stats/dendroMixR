## Figure 4: Model selection comparison (AIC, BIC, ICL, DSC) under
## skew-normal misspecification (Section 5.2, Misspecified Regime 2).
## True: 0.4*SN(5, 1, 20) + 0.6*SN(8, 1.5, 20)  [location, scale, skewness]
## Fitting: location-scale Gaussian (mclust "V") with k_over = 10 components.
## DSC uses rescaled formula (Remark 2): omega_n = log(n).

library(dendroMixR)
library(mclust)
library(ggplot2)
library(patchwork)
library(ggdendro)

# --- Parameters --------------------------------------------------------------

ps_true  <- c(0.4, 0.6)
mus_true <- c(5,   8)
sds_true <- c(1,   1.5)
alp_true <- c(20,  20)

k_max       <- 10
log10_n_seq <- seq(2, 4, by = 0.25)
n_rep       <- 100
n_demo      <- 2000

# --- Skew-normal helpers (no extra package needed) ---------------------------
# SN(mu, sigma, alpha): f(x) = (2/sigma) * phi((x-mu)/sigma) * Phi(alpha*(x-mu)/sigma)
# Generation: Z = delta*|W1| + sqrt(1-delta^2)*W2, X = mu + sigma*Z, delta = alpha/sqrt(1+alpha^2)

rskewnorm <- function(n, mu, sigma, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  z     <- delta * abs(rnorm(n)) + sqrt(1 - delta^2) * rnorm(n)
  mu + sigma * z
}

dskewnorm <- function(x, mu, sigma, alpha) {
  z <- (x - mu) / sigma
  (2 / sigma) * dnorm(z) * pnorm(alpha * z)
}

gen_data <- function(n) {
  comp <- sample(2L, n, replace = TRUE, prob = ps_true)
  x    <- numeric(n)
  for (g in 1:2) {
    idx <- comp == g
    if (any(idx))
      x[idx] <- rskewnorm(sum(idx), mus_true[g], sds_true[g], alp_true[g])
  }
  x
}

true_density <- function(x) {
  ps_true[1] * dskewnorm(x, mus_true[1], sds_true[1], alp_true[1]) +
  ps_true[2] * dskewnorm(x, mus_true[2], sds_true[2], alp_true[2])
}

# --- Helpers -----------------------------------------------------------------

log_lik_1d_mix <- function(x, G) {
  k   <- length(G$ps)
  mus <- as.vector(G$thetas)
  sds <- sqrt(pmax(as.vector(G$sigmas), 1e-10))
  log_mat <- matrix(0, length(x), k)
  for (j in seq_len(k))
    log_mat[, j] <- log(G$ps[j]) + dnorm(x, mus[j], sds[j], log = TRUE)
  row_max <- apply(log_mat, 1, max)
  mean(row_max + log(rowSums(exp(log_mat - row_max))))
}

# DSC (rescaled, Remark 2):  DSC(kappa) = -(omega_n * ll(kappa)/|ll(k)| + d(kappa)/max_d)
# d_n^(kappa) = merge height FROM kappa TO kappa-1 components = hc$height[k_over - kappa + 1]
# G_n^(kappa) = Gs[[k_over - kappa + 1]]
compute_dsc <- function(x, dmm, k_over) {
  omega_n  <- log(length(x))
  heights  <- dmm$hc$height
  max_h    <- max(heights)
  ll_kover <- log_lik_1d_mix(x, dmm$Gs[[1]])

  dsc_vals <- setNames(numeric(k_over - 1), as.character(2:k_over))
  for (kappa in 2:k_over) {
    G_kappa  <- dmm$Gs[[k_over - kappa + 1]]
    ll_kappa <- log_lik_1d_mix(x, G_kappa)
    dn_kappa <- heights[k_over - kappa + 1]
    dsc_vals[as.character(kappa)] <-
      -(omega_n * ll_kappa / abs(ll_kover) + dn_kappa / max_h)
  }
  as.integer(names(which.min(dsc_vals)))
}

select_aic_bic_icl <- function(x, k_max) {
  n <- length(x)
  aic_v <- bic_v <- icl_v <- rep(NA_real_, k_max)
  for (k in seq_len(k_max)) {
    fit <- tryCatch(Mclust(x, G = k, modelNames = "V", verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(fit) || length(fit$parameters$pro) != k) next
    bic_v[k] <- fit$bic
    aic_v[k] <- fit$bic + fit$df * (log(n) - 2)
    z   <- fit$z
    ent <- -sum(z * ifelse(z > 0, log(z), 0))
    icl_v[k] <- fit$bic - 2 * ent
  }
  c(aic = which.max(aic_v), bic = which.max(bic_v), icl = which.max(icl_v))
}

# --- Demo fit for panels (a) and (b) ----------------------------------------

set.seed(42)
x_demo   <- gen_data(n_demo)
fit_demo <- Mclust(x_demo, G = k_max, modelNames = "V", verbose = FALSE)
dmm_demo <- dendrogram_mixing(fit_demo$parameters$pro,
                               fit_demo$parameters$mean,
                               fit_demo$parameters$variance$sigmasq)

# --- Panel (a): Histogram + true density -------------------------------------

xgrid   <- seq(min(x_demo) - 0.5, max(x_demo) + 0.5, length.out = 500)
df_dens <- data.frame(x = xgrid, y = true_density(xgrid))

p_hist <- ggplot(data.frame(x = x_demo), aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "grey80", color = "white", linewidth = 0.2) +
  geom_line(data = df_dens, aes(x = x, y = y), linewidth = 0.9) +
  labs(x = "x", y = "Density", title = "(a) Histogram with true density") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# --- Panel (b): Dendrogram ---------------------------------------------------

ddata <- dendro_data(as.dendrogram(dmm_demo$hc), type = "rectangle")

p_dendro <- ggplot() +
  geom_segment(data = ggdendro::segment(ddata),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.5) +
  geom_text(data = ggdendro::label(ddata),
            aes(x = x, y = y, label = label),
            vjust = 1.5, size = 3) +
  labs(x = NULL, y = "Dissimilarity",
       title = "(b) Dendrogram of mixing measure") +
  theme_classic(base_size = 12) +
  theme(plot.title   = element_text(hjust = 0.5, size = 11),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

# --- Simulation for panel (c) ------------------------------------------------

results <- data.frame()
set.seed(2024)

for (log10_n in log10_n_seq) {
  n <- round(10^log10_n)
  cat(sprintf("n = %d (log10 = %.2f)\n", n, log10_n))

  for (rep in seq_len(n_rep)) {
    x   <- gen_data(n)
    abc <- tryCatch(select_aic_bic_icl(x, k_max),
                   error = function(e) c(aic = NA, bic = NA, icl = NA))

    fit_over <- tryCatch(Mclust(x, G = k_max, modelNames = "V", verbose = FALSE),
                         error = function(e) NULL)
    k_dsc <- NA_integer_
    if (!is.null(fit_over) && length(fit_over$parameters$pro) == k_max) {
      dmm   <- dendrogram_mixing(fit_over$parameters$pro,
                                 fit_over$parameters$mean,
                                 fit_over$parameters$variance$sigmasq)
      k_dsc <- tryCatch(compute_dsc(x, dmm, k_max), error = function(e) NA_integer_)
    }

    results <- rbind(results, data.frame(
      log10_n = log10_n, n = n, rep = rep,
      k_aic = abc["aic"], k_bic = abc["bic"],
      k_icl = abc["icl"], k_dsc = k_dsc
    ))
  }
}

# --- Panel (c): Average chosen k vs log10(n) ---------------------------------

agg <- aggregate(cbind(k_aic, k_bic, k_icl, k_dsc) ~ log10_n,
                 data = results,
                 FUN  = function(v) c(mean = mean(v, na.rm = TRUE),
                                      sd   = sd(v,   na.rm = TRUE)))

plot_df <- rbind(
  data.frame(log10_n = agg$log10_n,
             mean = agg$k_aic[, "mean"], sd = agg$k_aic[, "sd"], method = "AIC"),
  data.frame(log10_n = agg$log10_n,
             mean = agg$k_bic[, "mean"], sd = agg$k_bic[, "sd"], method = "BIC"),
  data.frame(log10_n = agg$log10_n,
             mean = agg$k_icl[, "mean"], sd = agg$k_icl[, "sd"], method = "ICL"),
  data.frame(log10_n = agg$log10_n,
             mean = agg$k_dsc[, "mean"], sd = agg$k_dsc[, "sd"], method = "DSC")
)
plot_df$method <- factor(plot_df$method, levels = c("AIC", "BIC", "ICL", "DSC"))

cols_m  <- c(AIC = "#1f77b4", BIC = "#ff7f0e", ICL = "#2ca02c", DSC = "#d62728")
shape_m <- c(AIC = 16L, BIC = 17L, ICL = 15L, DSC = 18L)

p_avg <- ggplot(plot_df, aes(x = log10_n, y = mean,
                              color = method, shape = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.06, linewidth = 0.4) +
  geom_hline(yintercept = 2, linetype = "dotted", color = "grey50", linewidth = 0.6) +
  scale_color_manual(values = cols_m) +
  scale_shape_manual(values = shape_m) +
  scale_y_continuous(limits = c(1, k_max), breaks = seq(2, k_max, by = 2)) +
  labs(x    = expression(log[10](n)),
       y    = "Average no. of components",
       title = "(c) Average choices no. of components",
       color = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title      = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(color = "grey50", fill = "white",
                                         linewidth = 0.3),
        legend.margin   = margin(3, 5, 3, 5))

# --- Combine and save --------------------------------------------------------

fig4 <- p_hist | p_dendro | p_avg

quartz(width = 15, height = 5)
print(fig4)

ggsave("code_reproduce_paper/R_fig4_skewnormal.pdf",
       plot = fig4, width = 15, height = 5)
cat("\nSaved: code_reproduce_paper/R_fig4_skewnormal.pdf\n")
