## Figure 6: 2D demonstration of overfitting + dendrogram merging.
## Simulate n = 2000 observations from a 3-component 2D Gaussian mixture,
## overfit with k = 10, merge via dendrogram_mixing() down to k = 3, and
## show that the merged components closely recover the true Gaussians.

library(dendroMixR)
library(mclust)
library(ggplot2)
library(MASS)
library(patchwork)

# --- Helpers -----------------------------------------------------------------

# 95% confidence ellipse for a 2D Gaussian (mean mu, covariance Sigma)
gauss_ellipse <- function(mu, Sigma, level = 0.95, n_pts = 200) {
  eig  <- eigen(Sigma)
  r    <- sqrt(qchisq(level, df = 2))
  angs <- seq(0, 2 * pi, length.out = n_pts)
  pts  <- t(mu + r * eig$vectors %*% diag(sqrt(eig$values)) %*% rbind(cos(angs), sin(angs)))
  data.frame(x1 = pts[, 1], x2 = pts[, 2])
}

build_ellipses <- function(mus, sigmas, ps) {
  do.call(rbind, lapply(seq_along(ps), function(j) {
    ell      <- gauss_ellipse(mus[j, ], sigmas[[j]])
    ell$comp <- factor(j)
    ell$wt   <- ps[j]
    ell
  }))
}

# --- True parameters (Section 5.1, weakly identifiable setting) -------------

set.seed(42)
n <- 2000

ps0     <- c(1/3, 1/3, 1/3)
mus0    <- matrix(c( 2,  1,
                      0,  6,
                     -2,  1), nrow = 3, byrow = TRUE)
sigmas0 <- list(
  matrix(c( 0.50,  0.50,  0.50,  1.00), 2, 2),
  matrix(c( 0.50, -0.10, -0.10,  0.10), 2, 2),
  matrix(c( 0.25,  0.50,  0.50,  2.00), 2, 2)
)
k0 <- 3

# --- Simulate data -----------------------------------------------------------

z <- sample(k0, n, replace = TRUE, prob = ps0)
x <- matrix(0, n, 2)
for (i in seq_len(n)) x[i, ] <- MASS::mvrnorm(1, mus0[z[i], ], sigmas0[[z[i]]])
df_data <- data.frame(x1 = x[, 1], x2 = x[, 2])

# --- Overfitted MLE: k = 10 --------------------------------------------------

k_over <- 10
cat(sprintf("Fitting k = %d GMM via mclust (VVV)...\n", k_over))
fit_raw     <- Mclust(x, G = k_over, modelNames = "VVV", verbose = FALSE)
ps_over     <- fit_raw$parameters$pro
mus_over    <- t(fit_raw$parameters$mean)           # k_over x 2
sigmas_over <- lapply(seq_len(k_over), function(j)
                 fit_raw$parameters$variance$sigma[, , j])

# --- Dendrogram merging: k = 10 → k = 3 ------------------------------------
# Gs[[1]] = initial 10-component state; Gs[[k_over - k0 + 1]] = 3-component state

cat(sprintf("Merging from k = %d down to k = %d via dendrogram_mixing()...\n",
            k_over, k0))
dmm           <- dendrogram_mixing(ps_over, mus_over, sigmas_over)
G_merged      <- dmm$Gs[[k_over - k0 + 1]]          # after 7 merge steps
ps_merged     <- G_merged$ps
mus_merged    <- matrix(G_merged$thetas, ncol = 2)   # 3 x 2
sigmas_merged <- G_merged$sigmas

# --- Match merged components to true (by nearest-mean, greedy) --------------

dist_mat   <- as.matrix(dist(rbind(mus_merged, mus0)))[seq_len(k0), (k0 + 1):(2 * k0)]
final_match <- integer(k0)   # final_match[j] = true component index for merged j
used        <- logical(k0)
for (j in seq_len(k0)) {
  for (cand in order(dist_mat[j, ])) {
    if (!used[cand]) { final_match[j] <- cand; used[cand] <- TRUE; break }
  }
}
# perm: true component j corresponds to merged component order(final_match)[j]
perm <- order(final_match)

# --- Print parameter comparison ----------------------------------------------

cat("\n=== Merged vs. True component parameters ===\n")
for (j in seq_len(k0)) {
  jm <- perm[j]
  cat(sprintf(
    "Component %d | weight: true=%.3f  merged=%.3f\n",
    j, ps0[j], ps_merged[jm]
  ))
  cat(sprintf(
    "             | mean:   true=(%.3f, %.3f)  merged=(%.3f, %.3f)\n",
    mus0[j, 1], mus0[j, 2], mus_merged[jm, 1], mus_merged[jm, 2]
  ))
}

# --- Build ellipse data frames -----------------------------------------------

ell_true   <- build_ellipses(mus0, sigmas0, ps0)
ell_over   <- build_ellipses(mus_over, sigmas_over, ps_over)
# reorder merged so component j matches true component j
ell_merged <- build_ellipses(mus_merged[perm, , drop = FALSE],
                              sigmas_merged[perm],
                              ps_merged[perm])

ctr_true   <- data.frame(x1 = mus0[, 1],                    x2 = mus0[, 2],                    comp = factor(1:k0))
ctr_over   <- data.frame(x1 = mus_over[, 1],                x2 = mus_over[, 2],                comp = factor(seq_len(k_over)))
ctr_merged <- data.frame(x1 = mus_merged[perm, 1],          x2 = mus_merged[perm, 2],          comp = factor(1:k0))

# --- Colors ------------------------------------------------------------------

cols3 <- c("#e41a1c", "#377eb8", "#4daf4a")   # red / blue / green for k = 3
cols10 <- colorRampPalette(c(
  "#e41a1c", "#ff7f00", "#d4b400",
  "#4daf4a", "#00bcd4", "#377eb8",
  "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"
))(k_over)

xlim_r <- range(x[, 1]) + c(-0.5, 0.5)
ylim_r <- range(x[, 2]) + c(-0.5, 0.5)

# --- Base plot (shared scatter layer) ----------------------------------------

base_plot <- ggplot(df_data, aes(x = x1, y = x2)) +
  geom_point(color = "grey78", alpha = 0.25, size = 0.4) +
  coord_cartesian(xlim = xlim_r, ylim = ylim_r) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title         = element_text(size = 11)
  ) +
  labs(x = expression(X[1]), y = expression(X[2]))

# --- Panel 1: True mixture (k = 3) ------------------------------------------

p1 <- base_plot +
  geom_path(data  = ell_true,
            aes(x = x1, y = x2, group = comp, color = comp),
            linewidth = 1.3) +
  geom_point(data  = ctr_true,
             aes(x = x1, y = x2, color = comp),
             shape = 3, size = 5, stroke = 2.2) +
  scale_color_manual(values = setNames(cols3, 1:k0)) +
  ggtitle("True mixture  (k = 3)")

# --- Panel 2: Overfitted MLE (k = 10) ----------------------------------------

p2 <- base_plot +
  geom_path(data  = ell_over,
            aes(x = x1, y = x2, group = comp, color = comp),
            linewidth = 0.8, linetype = "dashed") +
  geom_point(data  = ctr_over,
             aes(x = x1, y = x2, color = comp),
             shape = 3, size = 3, stroke = 1.5) +
  scale_color_manual(values = setNames(cols10, seq_len(k_over))) +
  ggtitle("Overfitted MLE  (k = 10)")

# --- Panel 3: Merged to k = 3 ------------------------------------------------

p3 <- base_plot +
  geom_path(data  = ell_merged,
            aes(x = x1, y = x2, group = comp, color = comp),
            linewidth = 1.3) +
  geom_point(data  = ctr_merged,
             aes(x = x1, y = x2, color = comp),
             shape = 3, size = 5, stroke = 2.2) +
  scale_color_manual(values = setNames(cols3, 1:k0)) +
  ggtitle("Merged  (k = 3)")

# --- Combine and save --------------------------------------------------------

fig <- p1 | p2 | p3

quartz(width = 13, height = 4.5)
print(fig)

ggsave("code_reproduce_paper/2d-demonstration.pdf", plot = fig, width = 13, height = 4.5)
cat("\nSaved: code_reproduce_paper/2d-demonstration.pdf\n")
