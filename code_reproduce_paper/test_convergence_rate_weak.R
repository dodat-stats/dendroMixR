## Verify Theorem 4 (weakly identifiable setting, location-scale Gaussian):
## Overfitted MLE has slow W6 rate ~n^{-1/12}, but the merged measure
## and exact-fitted MLE both achieve the fast W1 rate ~n^{-1/2}.
## Reproduces Fig. 2(b) of the paper.

library(dendroMixR)
library(mclust)

# --- Wasserstein distance for location-scale Gaussian mixtures ------------
# Ground metric: ||(mu_i, Sigma_i) - (mu_j, Sigma_j)|| =
#   sqrt(||mu_i - mu_j||^2 + ||Sigma_i - Sigma_j||_F)

wasserstein_r_gauss <- function(ps, mus, sigmas, ps0, mus0, sigmas0, r = 1) {
     k <- length(ps)
     k0 <- length(ps0)
     cost <- matrix(0, k, k0)
     for (i in seq_len(k)) {
          for (j in seq_len(k0)) {
               d_mu <- sqrt(sum((mus[i, ] - mus0[j, ])^2))
               d_sig <- sqrt(sum((sigmas[[i]] - sigmas0[[j]])^2))
               cost[i, j] <- (d_mu + d_sig)^r
          }
     }
     f.obj <- as.vector(t(cost))
     A_row <- matrix(0, k, k * k0)
     for (i in seq_len(k)) {
          cols <- (i - 1) * k0 + seq_len(k0)
          A_row[i, cols] <- 1
     }
     A_col <- matrix(0, k0, k * k0)
     for (j in seq_len(k0)) {
          cols <- seq(j, k * k0, by = k0)
          A_col[j, cols] <- 1
     }
     A <- rbind(A_row, A_col)
     b <- c(ps, ps0)
     dir <- rep("=", k + k0)
     sol <- lpSolve::lp("min", f.obj, A, dir, b)
     return(sol$objval^(1 / r))
}

# --- Helper: extract mclust parameters into our format -------------------

extract_mclust <- function(fit) {
     k <- length(fit$parameters$pro)
     d <- nrow(fit$parameters$mean)
     ps <- fit$parameters$pro
     mus <- t(fit$parameters$mean)  # k x d
     sigmas <- lapply(seq_len(k), function(j) fit$parameters$variance$sigma[, , j])
     list(ps = ps, mus = mus, sigmas = sigmas)
}

# --- Simulation setup (matches Section 5.1, weak identifiability) --------

ps0 <- c(1/3, 1/3, 1/3)
mus0 <- matrix(c(2, 1,
                  0, 6,
                  -2, 1), nrow = 3, byrow = TRUE)
sigmas0 <- list(
     matrix(c(0.5, 0.5, 0.5, 1.0), 2, 2),
     matrix(c(0.5, -0.1, -0.1, 0.1), 2, 2),
     matrix(c(0.25, 0.5, 0.5, 2.0), 2, 2)
)

k0 <- 3
k_over <- 5
# r(k - k0 + 1) = r(3) = 6, so use W6 for overfitted
r_over <- 6

log10_n_seq <- seq(2, 4, by = 0.25)
n_rep <- 64

results <- data.frame()

set.seed(2024)

for (log10_n in log10_n_seq) {
     n <- round(10^log10_n)
     cat(sprintf("n = %d (log10 = %.2f)\n", n, log10_n))

     for (rep in seq_len(n_rep)) {
          # generate data from true mixture of location-scale Gaussians
          z <- sample(k0, n, replace = TRUE, prob = ps0)
          x <- matrix(0, n, 2)
          for (i in seq_len(n)) {
               x[i, ] <- MASS::mvrnorm(1, mu = mus0[z[i], ], Sigma = sigmas0[[z[i]]])
          }

          # exact-fitted MLE (k = k0 = 3) using mclust
          fit_exact_raw <- Mclust(x, G = k0, modelNames = "VVV", verbose = FALSE)
          if (is.null(fit_exact_raw)) next
          fit_exact <- extract_mclust(fit_exact_raw)

          # overfitted MLE (k = 5) using mclust
          fit_over_raw <- Mclust(x, G = k_over, modelNames = "VVV", verbose = FALSE)
          if (is.null(fit_over_raw)) next
          fit_over <- extract_mclust(fit_over_raw)

          # merge overfitted down to k0 = 3 atoms using weak identifiability
          dmm <- dendrogram_mixing(fit_over$ps, fit_over$mus, fit_over$sigmas)
          G_merged <- dmm$Gs[[k_over - k0 + 1]]
          merged_mus <- matrix(G_merged$thetas, ncol = 2)
          merged_sigmas <- G_merged$sigmas

          # compute distances
          w6_over <- wasserstein_r_gauss(fit_over$ps, fit_over$mus, fit_over$sigmas,
                                         ps0, mus0, sigmas0, r = r_over)
          w1_over <- wasserstein_r_gauss(fit_over$ps, fit_over$mus, fit_over$sigmas,
                                         ps0, mus0, sigmas0, r = 1)
          w1_exact <- wasserstein_r_gauss(fit_exact$ps, fit_exact$mus, fit_exact$sigmas,
                                          ps0, mus0, sigmas0, r = 1)
          w1_merged <- wasserstein_r_gauss(G_merged$ps, merged_mus, merged_sigmas,
                                           ps0, mus0, sigmas0, r = 1)

          results <- rbind(results, data.frame(
               log10_n = log10_n,
               n = n,
               rep = rep,
               w6_over = w6_over,
               w1_over = w1_over,
               w1_exact = w1_exact,
               w1_merged = w1_merged
          ))
     }
}

# --- Aggregate and plot (reproduces Fig. 2(b)) ----------------------------

agg <- aggregate(
     cbind(w6_over, w1_over, w1_exact, w1_merged) ~ log10_n,
     data = results,
     FUN = function(x) c(mean = mean(log10(x)),
                         q25 = quantile(log10(x), 0.25),
                         q75 = quantile(log10(x), 0.75))
)

log10_n <- agg$log10_n
mean_w6_over <- agg$w6_over[, "mean"]
mean_w1_over <- agg$w1_over[, "mean"]
mean_w1_exact <- agg$w1_exact[, "mean"]
mean_w1_merged <- agg$w1_merged[, "mean"]

q25_w6_over <- agg$w6_over[, "q25.25%"]
q75_w6_over <- agg$w6_over[, "q75.75%"]
q25_w1_over <- agg$w1_over[, "q25.25%"]
q75_w1_over <- agg$w1_over[, "q75.75%"]
q25_w1_exact <- agg$w1_exact[, "q25.25%"]
q75_w1_exact <- agg$w1_exact[, "q75.75%"]
q25_w1_merged <- agg$w1_merged[, "q25.25%"]
q75_w1_merged <- agg$w1_merged[, "q75.75%"]

# fit linear regression for slope (expected: -1/12 ~ -0.08 for W6, -0.5 for W1 exact/merged)
fit_w6_over <- lm(mean_w6_over ~ log10_n)
fit_w1_over <- lm(mean_w1_over ~ log10_n)
fit_w1_exact <- lm(mean_w1_exact ~ log10_n)
fit_w1_merged <- lm(mean_w1_merged ~ log10_n)

cat("\n=== Fitted slopes (theory: W6 overfitted ~ -0.08, W1 exact/merged ~ -0.50) ===\n")
cat(sprintf("W6(Go_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w6_over)[2], coef(fit_w6_over)[1]))
cat(sprintf("W1(Go_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_over)[2], coef(fit_w1_over)[1]))
cat(sprintf("W1(Ge_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_exact)[2], coef(fit_w1_exact)[1]))
cat(sprintf("W1(Gm_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_merged)[2], coef(fit_w1_merged)[1]))

# --- Plot with ggplot2 + latex2exp ----------------------------------------

library(ggplot2)
library(latex2exp)

col_exact <- "#d62728"    # red
col_w6_over <- "#1f77b4"  # blue
col_w1_over <- "#9467bd"  # purple
col_merged <- "#8c564b"   # brown

plot_df <- data.frame(
     log10_n = rep(log10_n, 4),
     mean = c(mean_w1_exact, mean_w6_over, mean_w1_over, mean_w1_merged),
     q25 = c(q25_w1_exact, q25_w6_over, q25_w1_over, q25_w1_merged),
     q75 = c(q75_w1_exact, q75_w6_over, q75_w1_over, q75_w1_merged),
     group = rep(c("w1_exact", "w6_over", "w1_over", "w1_merged"), each = length(log10_n))
)

s_e <- coef(fit_w1_exact); s_w6 <- coef(fit_w6_over)
s_w1 <- coef(fit_w1_over); s_m <- coef(fit_w1_merged)

labels <- c(
     w1_exact = sprintf("$\\log_{10}(W_1(\\hat{G}_n^e,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_e[2], s_e[1]),
     w6_over  = sprintf("$\\log_{10}(W_6(\\hat{G}_n^o,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_w6[2], s_w6[1]),
     w1_over  = sprintf("$\\log_{10}(W_1(\\hat{G}_n^o,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_w1[2], s_w1[1]),
     w1_merged = sprintf("$\\log_{10}(W_1(\\hat{G}_n^m,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_m[2], s_m[1])
)

plot_df$group <- factor(plot_df$group, levels = c("w1_exact", "w6_over", "w1_over", "w1_merged"))

p <- ggplot(plot_df, aes(x = log10_n, y = mean, color = group, shape = group)) +
     geom_point(size = 2.5) +
     geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.05, linewidth = 0.4) +
     geom_abline(intercept = s_e[1], slope = s_e[2], color = col_exact, linetype = "dashed", linewidth = 0.6) +
     geom_abline(intercept = s_w6[1], slope = s_w6[2], color = col_w6_over, linetype = "dotdash", linewidth = 0.6) +
     geom_abline(intercept = s_w1[1], slope = s_w1[2], color = col_w1_over, linetype = "dotdash", linewidth = 0.6) +
     geom_abline(intercept = s_m[1], slope = s_m[2], color = col_merged, linetype = "solid", linewidth = 0.6) +
     scale_color_manual(
          values = c(w1_exact = col_exact, w6_over = col_w6_over,
                     w1_over = col_w1_over, w1_merged = col_merged),
          labels = lapply(labels, TeX)
     ) +
     scale_shape_manual(
          values = c(w1_exact = 16, w6_over = 17, w1_over = 18, w1_merged = 15),
          labels = lapply(labels, TeX)
     ) +
     labs(x = TeX("$\\log_{10}(n)$"),
          y = TeX("$\\log_{10}(error)$"),
          # title = "Convergence rates (weakly identifiable)",
          color = NULL, shape = NULL) +
     theme_classic(base_size = 13) +
     theme(
       legend.position = c(0.565, 0.01),
          legend.justification = c(1, 0),
          legend.direction = "vertical",
          legend.text = element_text(size = 11),
          legend.background = element_rect(color = "grey50", fill = "white", linewidth = 0.3),
          legend.margin = margin(4, 6, 4, 6),
          plot.title = element_text(face = "bold", hjust = 0.5)
     ) +
     guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

quartz(width = 7, height = 6)
print(p)
p

ggsave("test/fig_convergence_weak.pdf", plot = p, width = 7, height = 6)
