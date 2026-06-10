## Verify Theorem 1 (strongly identifiable setting):
## Overfitted MLE has slow W2 rate ~n^{-1/4}, but the merged measure
## and exact-fitted MLE both achieve the fast W1 rate ~n^{-1/2}.
## Reproduces Fig. 2(a) of the paper.

library(dendroMixR)

# --- Wasserstein distance via optimal transport (linear program) -----------
# For discrete measures G = sum p_i delta_{theta_i}, G0 = sum p0_j delta_{theta0_j}
# W_r^r = min_{q in Pi(p, p0)} sum_{i,j} q_{ij} ||theta_i - theta0_j||^r

wasserstein_r <- function(ps, thetas, ps0, thetas0, r = 1) {
     thetas <- as.matrix(thetas)
     thetas0 <- as.matrix(thetas0)
     k <- length(ps)
     k0 <- length(ps0)
     # cost matrix: ||theta_i - theta0_j||^r
     cost <- matrix(0, k, k0)
     for (i in seq_len(k)) {
          for (j in seq_len(k0)) {
               cost[i, j] <- sqrt(sum((thetas[i, ] - thetas0[j, ])^2))^r
          }
     }
     # solve the transportation problem via LP
     # min sum c_{ij} q_{ij}
     # s.t. sum_j q_{ij} = p_i, sum_i q_{ij} = p0_j, q_{ij} >= 0
     # flatten q into a vector of length k * k0
     f.obj <- as.vector(t(cost))
     # row-sum constraints: sum_j q_{ij} = p_i for each i
     A_row <- matrix(0, k, k * k0)
     for (i in seq_len(k)) {
          cols <- (i - 1) * k0 + seq_len(k0)
          A_row[i, cols] <- 1
     }
     # col-sum constraints: sum_i q_{ij} = p0_j for each j
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

# --- EM for location Gaussian with KNOWN covariance (= Identity) ----------

em_location_gaussian <- function(x, k, max_iter = 100, tol = 1e-8, n_init = 10) {
     n <- nrow(x)
     d <- ncol(x)
     best_ll <- -Inf
     best_result <- NULL

     for (init in seq_len(n_init)) {
          # random initialization: pick k data points as centers
          idx <- sample(n, k)
          mu <- x[idx, , drop = FALSE]
          p <- rep(1 / k, k)

          ll_old <- -Inf
          for (iter in seq_len(max_iter)) {
               # E-step: compute responsibilities
               log_resp <- matrix(0, n, k)
               for (j in seq_len(k)) {
                    diff <- sweep(x, 2, mu[j, ])
                    log_resp[, j] <- log(p[j]) - 0.5 * rowSums(diff^2)
               }
               # log-sum-exp for numerical stability
               log_max <- apply(log_resp, 1, max)
               log_sum <- log_max + log(rowSums(exp(log_resp - log_max)))
               resp <- exp(log_resp - log_sum)

               ll_new <- sum(log_sum) - n * d / 2 * log(2 * pi)

               # M-step
               nk <- colSums(resp)
               p <- nk / n
               for (j in seq_len(k)) {
                    mu[j, ] <- colSums(resp[, j] * x) / nk[j]
               }

               if (abs(ll_new - ll_old) < tol) break
               ll_old <- ll_new
          }

          if (ll_new > best_ll) {
               best_ll <- ll_new
               best_result <- list(ps = p, thetas = mu)
          }
     }
     best_result
}

# --- Simulation setup (matches Section 5.1, strong identifiability) -------

# True mixing measure: 3 components, 2D, known covariance = I
ps0 <- c(1/3, 1/3, 1/3)
thetas0 <- matrix(c(2, 1,
                     0, 6,
                     -2, 1), nrow = 3, byrow = TRUE)
k0 <- 3
k_over <- 5

log10_n_seq <- seq(2, 4, by = 0.25)
n_rep <- 64

results <- data.frame()

set.seed(2024)

for (log10_n in log10_n_seq) {
     n <- round(10^log10_n)
     cat(sprintf("n = %d (log10 = %.2f)\n", n, log10_n))

     for (rep in seq_len(n_rep)) {
          # generate data from true mixture
          z <- sample(k0, n, replace = TRUE, prob = ps0)
          x <- matrix(rnorm(n * 2), n, 2)  # standard normal
          for (i in seq_len(n)) {
               x[i, ] <- x[i, ] + thetas0[z[i], ]
          }

          # exact-fitted MLE (k = k0 = 3)
          fit_exact <- em_location_gaussian(x, k = k0)
          # overfitted MLE (k = 5)
          fit_over <- em_location_gaussian(x, k = k_over)

          # merge overfitted down to k0 = 3 atoms
          dmm <- dendrogram_mixing(fit_over$ps, fit_over$thetas)
          # Gs[[1]] has k_over atoms, Gs[[k_over - k0 + 1]] has k0 atoms
          G_merged <- dmm$Gs[[k_over - k0 + 1]]

          # compute distances
          w2_over <- wasserstein_r(fit_over$ps, fit_over$thetas, ps0, thetas0, r = 2)
          w1_over <- wasserstein_r(fit_over$ps, fit_over$thetas, ps0, thetas0, r = 1)
          w1_exact <- wasserstein_r(fit_exact$ps, fit_exact$thetas, ps0, thetas0, r = 1)
          w1_merged <- wasserstein_r(G_merged$ps, as.matrix(G_merged$thetas), ps0, thetas0, r = 1)

          results <- rbind(results, data.frame(
               log10_n = log10_n,
               n = n,
               rep = rep,
               w2_over = w2_over,
               w1_over = w1_over,
               w1_exact = w1_exact,
               w1_merged = w1_merged
          ))
     }
}

# --- Aggregate and plot (reproduces Fig. 2(a)) ----------------------------

agg <- aggregate(
     cbind(w2_over, w1_over, w1_exact, w1_merged) ~ log10_n,
     data = results,
     FUN = function(x) c(mean = mean(log10(x)),
                         q25 = quantile(log10(x), 0.25),
                         q75 = quantile(log10(x), 0.75))
)

log10_n <- agg$log10_n
mean_w2_over <- agg$w2_over[, "mean"]
mean_w1_over <- agg$w1_over[, "mean"]
mean_w1_exact <- agg$w1_exact[, "mean"]
mean_w1_merged <- agg$w1_merged[, "mean"]

q25_w2_over <- agg$w2_over[, "q25.25%"]
q75_w2_over <- agg$w2_over[, "q75.75%"]
q25_w1_over <- agg$w1_over[, "q25.25%"]
q75_w1_over <- agg$w1_over[, "q75.75%"]
q25_w1_exact <- agg$w1_exact[, "q25.25%"]
q75_w1_exact <- agg$w1_exact[, "q75.75%"]
q25_w1_merged <- agg$w1_merged[, "q25.25%"]
q75_w1_merged <- agg$w1_merged[, "q75.75%"]

# fit linear regression for slope (expected: -0.25 for W2 overfitted, -0.5 for W1 exact/merged)
fit_w2_over <- lm(mean_w2_over ~ log10_n)
fit_w1_over <- lm(mean_w1_over ~ log10_n)
fit_w1_exact <- lm(mean_w1_exact ~ log10_n)
fit_w1_merged <- lm(mean_w1_merged ~ log10_n)

cat("\n=== Fitted slopes (theory: W2 overfitted ~ -0.25, W1 exact/merged ~ -0.50) ===\n")
cat(sprintf("W2(Go_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w2_over)[2], coef(fit_w2_over)[1]))
cat(sprintf("W1(Go_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_over)[2], coef(fit_w1_over)[1]))
cat(sprintf("W1(Ge_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_exact)[2], coef(fit_w1_exact)[1]))
cat(sprintf("W1(Gm_n, G0):  slope = %.2f, intercept = %.2f\n",
            coef(fit_w1_merged)[2], coef(fit_w1_merged)[1]))

# --- Plot with ggplot2 + latex2exp (reproduces Fig. 2(a)) -----------------

library(ggplot2)
library(latex2exp)

col_exact <- "#d62728"    # red
col_w2_over <- "#1f77b4"  # blue
col_w1_over <- "#9467bd"  # purple
col_merged <- "#8c564b"   # brown

# build data frame for points + error bars
plot_df <- data.frame(
     log10_n = rep(log10_n, 4),
     mean = c(mean_w1_exact, mean_w2_over, mean_w1_over, mean_w1_merged),
     q25 = c(q25_w1_exact, q25_w2_over, q25_w1_over, q25_w1_merged),
     q75 = c(q75_w1_exact, q75_w2_over, q75_w1_over, q75_w1_merged),
     group = rep(c("w1_exact", "w2_over", "w1_over", "w1_merged"), each = length(log10_n))
)

# fitted slopes and intercepts
s_e <- coef(fit_w1_exact); s_w2 <- coef(fit_w2_over)
s_w1 <- coef(fit_w1_over); s_m <- coef(fit_w1_merged)

# legend labels with LaTeX
labels <- c(
     w1_exact = sprintf("$\\log_{10}(W_1(\\hat{G}_n^e,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_e[2], s_e[1]),
     w2_over  = sprintf("$\\log_{10}(W_2(\\hat{G}_n^o,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_w2[2], s_w2[1]),
     w1_over  = sprintf("$\\log_{10}(W_1(\\hat{G}_n^o,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_w1[2], s_w1[1]),
     w1_merged = sprintf("$\\log_{10}(W_1(\\hat{G}_n^m,\\, G_0)) = %.2f \\, \\log_{10}(n) + %.2f$", s_m[2], s_m[1])
)

plot_df$group <- factor(plot_df$group, levels = c("w1_exact", "w2_over", "w1_over", "w1_merged"))

p <- ggplot(plot_df, aes(x = log10_n, y = mean, color = group, shape = group)) +
     geom_point(size = 2.5) +
     geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.05, linewidth = 0.4) +
     geom_abline(intercept = s_e[1], slope = s_e[2], color = col_exact, linetype = "dashed", linewidth = 0.6) +
     geom_abline(intercept = s_w2[1], slope = s_w2[2], color = col_w2_over, linetype = "dotdash", linewidth = 0.6) +
     geom_abline(intercept = s_w1[1], slope = s_w1[2], color = col_w1_over, linetype = "dotdash", linewidth = 0.6) +
     geom_abline(intercept = s_m[1], slope = s_m[2], color = col_merged, linetype = "solid", linewidth = 0.6) +
     scale_color_manual(
          values = c(w1_exact = col_exact, w2_over = col_w2_over,
                     w1_over = col_w1_over, w1_merged = col_merged),
          labels = lapply(labels, TeX)
     ) +
     scale_shape_manual(
          values = c(w1_exact = 16, w2_over = 17, w1_over = 18, w1_merged = 15),
          labels = lapply(labels, TeX)
     ) +
     labs(x = TeX("$\\log_{10}(n)$"),
          y = TeX("$\\log_{10}(error)$"),
          # title = "Convergence rates (strongly identifiable)",
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
print(p) # preview

ggsave("test/fig_convergence_strong.pdf", plot = p, width = 7, height = 6)
