#' Plot dendrogram of mixing measures
#'
#' @param dmm Output from \code{dendrogram_mixing}.
#' @param dim Dimension to plot.
#' @param point_size Point size.
#' @param line_width Line width.
#' @param palette Vector of colors.
#' @param show_legend Show colors associated with each step in the dendrogram.
#' @param main Title name
#' @return A base R plot.
#' @export
plot_dendrogram_mixing <- function(dmm, dim = 1,
                                      point_size = 3,
                                      line_width = 0.5,
                                      palette = NULL,
                                      show_legend = FALSE,
                                      main = NULL
) {
     stopifnot(dim >= 1)
     if (!requireNamespace("ggplot2", quietly = TRUE)) {
          stop("Package 'ggplot2' is required for this function. Install it with install.packages('ggplot2').")
     }

     Gs         <- dmm$Gs
     heights    <- dmm$hc$height
     merge_pair <- dmm$merge_pair
     K_bar      <- length(Gs)

     # cumulative x positions
     heights_cumsum <- cumsum(c(0, heights)[-(K_bar + 1)])

     # collect theta matrices
     Theta <- vector("list", K_bar)
     for (s in seq_len(K_bar)) {
          th <- as.matrix(Gs[[s]]$thetas)
          if (s == K_bar) th <- t(th)
          Theta[[s]] <- th
     }

     # y range
     y_all <- unlist(lapply(Theta, function(M) M[, dim]))
     y_rng <- range(y_all, finite = TRUE)

     # ---- Build segments (gray branches) ----
     segs <- list()
     for (step in seq_len(K_bar - 1)) {
          th_now  <- Theta[[step]]
          th_next <- Theta[[step + 1]]

          i_merge <- merge_pair[step, 1]
          j_merge <- merge_pair[step, 2]

          x1 <- heights_cumsum[step]
          x2 <- heights_cumsum[step + 1]

          # left block
          if (j_merge - 1 >= 1) {
               idx <- seq_len(j_merge - 1)
               y1  <- th_now[idx,  dim]
               y2  <- th_next[idx, dim]
               segs[[length(segs) + 1]] <- data.frame(x = x1, xend = x2, y = y1, yend = y1)
               segs[[length(segs) + 1]] <- data.frame(x = x2, xend = x2, y = y1, yend = y2)
          }

          # mid block (merge pair)
          y1 <- th_now[j_merge, dim]
          y2 <- th_next[i_merge, dim]
          segs[[length(segs) + 1]] <- data.frame(x = x1, xend = x2, y = y1, yend = y1)
          segs[[length(segs) + 1]] <- data.frame(x = x2, xend = x2, y = y1, yend = y2)

          # right block
          right_max <- K_bar - step + 1
          if (j_merge != right_max && step < K_bar - 1) {
               idx <- seq.int(j_merge + 1, right_max)
               y1  <- th_now[idx, dim]
               y2  <- th_next[idx - 1, dim]
               segs[[length(segs) + 1]] <- data.frame(x = x1, xend = x2, y = y1, yend = y1)
               segs[[length(segs) + 1]] <- data.frame(x = x2, xend = x2, y = y1, yend = y2)
          }
     }
     seg_df <- if (length(segs)) do.call(rbind, segs) else
          data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric())

     # ---- Points (colored) ----
     pts <- do.call(
          rbind,
          lapply(seq_len(K_bar), function(s) {
               data.frame(
                    x = rep(heights_cumsum[s], nrow(Theta[[s]])),
                    y = Theta[[s]][, dim],
                    step = s,
                    idx = seq_len(nrow(Theta[[s]]))
               )
          })
     )

     # categorical palette
     if (is.null(palette)) {
          # falls back to base HCL if randomcoloR unavailable
          if (requireNamespace("randomcoloR", quietly = TRUE)) {
               palette <- randomcoloR::distinctColorPalette(K_bar)
          } else {
               palette <- hcl.colors(K_bar, "Dark3")  # nice categorical palette
          }
     }

     # labels for the first step (like your text() call)
     lab_df <- data.frame(
          x = rep(heights_cumsum[1], nrow(Theta[[1]])),
          y = Theta[[1]][, dim],
          label = seq_len(nrow(Theta[[1]]))
     )


     p <- ggplot2::ggplot() +
          # gray dendrogram branches
          ggplot2::geom_segment(
               data = seg_df,
               ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = line_width,
               color = "grey60"
          ) +
          ggplot2::geom_point(
               data = pts,
               ggplot2::aes(x = x, y = y, color = factor(step)),
               size = point_size
          ) +
          ggplot2::geom_text(
               data = lab_df,
               ggplot2::aes(x = x, y = y, label = label),
               color = "white", size = max(2, point_size * 0.9)
          ) +
          ggplot2::scale_color_manual(values = palette) +
          ggplot2::labs(x = "Heights", y = paste0("theta[", dim, "]"),
               title = main,
               color = NULL) +
          ggplot2::coord_cartesian(ylim = y_rng, expand = TRUE) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
               panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
          ) +
          ggplot2::theme(
               plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
               legend.position = if (show_legend) "right" else "none"
          )

     p
}
