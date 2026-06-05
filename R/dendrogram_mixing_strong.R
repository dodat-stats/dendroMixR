#' Dendrogram Mixing (Strong Identifiability, 1D)
#' @noRd
dendrogram_mixing_1d_strong <- function(ps, thetas){
     ## univariate case
     n <- length(thetas)
     if (n < 2) stop("need n >= 2")
     harmonic_means = outer(ps, ps, "*") / outer(ps, ps, "+")
     origD <- harmonic_means * (outer(thetas, thetas, "-"))^2
     # dynamic cluster distance matrix (between current clusters)
     D <- origD
     diag(D) <- Inf
     # clusters hold the member indices for each current cluster
     clusters = lapply(c(1:n), function(i) list(members=i, p=ps[i], theta=thetas[i]))
     # cluster_labels hold numbers to be placed in merge: negative for leaves, positive for merge-row indices
     cluster_labels <- - (1:n)
     merge <- matrix(0L, n - 1, 2L)
     merge_pair <- matrix(0L, n - 1, 2L)
     height <- numeric(n - 1)
     G = list(ps = ps, thetas = thetas)
     Gs = list(G)
     for (step in seq_len(n - 2)) {
          # pick closest pair (first occurrence if ties)
          idx <- which(D == min(D), arr.ind = TRUE)[1, ]
          i <- idx[1]; j <- idx[2]
          if (i > j) { tmp <- i; i <- j; j <- tmp }
          # put the labels (negative leaves or positive merge-row numbers) into merge[step, ]
          merge[step, ] <- c(cluster_labels[i], cluster_labels[j])
          height[step] <- D[i, j]
          merge_pair[step, ] <- c(i, j)
          # create merged cluster
          new_members <- c(clusters[[i]]$members, clusters[[j]]$members)
          new_p = clusters[[i]]$p + clusters[[j]]$p
          new_theta = (clusters[[i]]$p * clusters[[i]]$theta + clusters[[j]]$p * clusters[[j]]$theta) / new_p
          new_cluster <- list(members=new_members, p=new_p, theta=new_theta)
          clusters[[i]] <- new_cluster
          clusters[[j]] <- NULL
          # this merged cluster is formed at row 'step' of merge
          cluster_labels[i] <- step
          cluster_labels <- cluster_labels[-j]
          ps <- ps[-j]
          ps[i] = new_p
          thetas <- (thetas[-j])
          thetas[i] = new_theta
          new_G = list(ps = ps, thetas = thetas)
          Gs = append(Gs, list(new_G))
          # remove j-th row/col from D
          # and recompute distances between new cluster i and all others
          D <- D[-j, -j, drop = FALSE]
          hm = (ps[i] * ps[-i]) / (ps[i] + ps[-i])
          dm = ((thetas[i] - thetas[-i])^2)
          D[i, -i] = hm * dm
          D[-i, i] = hm * dm
          D[i, i] = Inf
     }
     ## last step: only merge and don't need compute dissimilarity matrix
     step = n-1
     idx <- which(D == min(D), arr.ind = TRUE)[1, ]
     i <- idx[1]; j <- idx[2]
     if (i > j) { tmp <- i; i <- j; j <- tmp }
     # put the labels (negative leaves or positive merge-row numbers) into merge[step, ]
     merge[step, ] <- c(cluster_labels[i], cluster_labels[j])
     height[step] <- D[i, j]
     merge_pair[step, ] <- c(i, j)
     # create merged cluster
     new_members <- c(clusters[[i]]$members, clusters[[j]]$members)
     new_p = clusters[[i]]$p + clusters[[j]]$p
     new_theta = (clusters[[i]]$p * clusters[[i]]$theta + clusters[[j]]$p * clusters[[j]]$theta) / new_p
     new_cluster <- list(members=new_members, p=new_p, theta=new_theta)
     clusters[[i]] <- new_cluster
     clusters[[j]] <- NULL
     # this merged cluster is formed at row 'step' of merge
     cluster_labels[i] <- step
     cluster_labels <- cluster_labels[-j]
     ps <- ps[-j]
     ps[i] = new_p
     thetas <- thetas[-j]
     thetas[i] = new_theta
     new_G = list(ps = ps, thetas = thetas)
     Gs = append(Gs, list(new_G))
     # recursive order builder: positive integers are merge-row indices (1..n-1), negative are leaves
     get_order <- function(node) {
          if (node < 0) return(-node)
          left <- merge[node, 1]
          right <- merge[node, 2]
          c(get_order(left), get_order(right))
     }
     order <- get_order(n - 1)  # root is the last merge row
     labels <- if (!is.null(rownames(Gs[[1]]$thetas))) rownames(Gs[[1]]$thetas) else as.character(seq_len(n))
     hc = structure(
          list(merge = merge,
               height = height,
               order = order,
               labels = labels,
               call = match.call()),
          class = "hclust"
     )
     return(list(hc=hc, Gs=Gs, merge_pair=merge_pair))
}



#' Dendrogram Mixing (Strong Identifiability, Multi-D)
#' @noRd
dendrogram_mixing_multid_strong <- function(ps, thetas){
     ## multivariate case
     n <- nrow(thetas)
     if (n < 2) stop("need n >= 2")
     harmonic_means = outer(ps, ps, "*") / outer(ps, ps, "+")
     origD <- harmonic_means * as.matrix(dist(thetas)^2)
     # dynamic cluster distance matrix (between current clusters)
     D <- origD
     diag(D) <- Inf
     # clusters hold the member indices for each current cluster
     clusters = lapply(c(1:n), function(i) list(members=i, p=ps[i], theta=thetas[i, ]))
     # cluster_labels hold numbers to be placed in merge: negative for leaves, positive for merge-row indices
     cluster_labels <- - (1:n)
     merge <- matrix(0L, n - 1, 2L)
     merge_pair <- matrix(0L, n - 1, 2L)
     height <- numeric(n - 1)
     G = list(ps = ps, thetas = thetas)
     Gs = list(G)
     for (step in seq_len(n - 2)) {
          # pick closest pair (first occurrence if ties)
          idx <- which(D == min(D), arr.ind = TRUE)[1, ]
          i <- idx[1]; j <- idx[2]
          if (i > j) { tmp <- i; i <- j; j <- tmp }
          # put the labels (negative leaves or positive merge-row numbers) into merge[step, ]
          merge[step, ] <- c(cluster_labels[i], cluster_labels[j])
          height[step] <- D[i, j]
          merge_pair[step, ] <- c(i, j)
          # create merged cluster
          new_members <- c(clusters[[i]]$members, clusters[[j]]$members)
          new_p = clusters[[i]]$p + clusters[[j]]$p
          new_theta = (clusters[[i]]$p * clusters[[i]]$theta + clusters[[j]]$p * clusters[[j]]$theta) / new_p
          new_cluster <- list(members=new_members, p=new_p, theta=new_theta)
          clusters[[i]] <- new_cluster
          clusters[[j]] <- NULL
          # this merged cluster is formed at row 'step' of merge
          cluster_labels[i] <- step
          cluster_labels <- cluster_labels[-j]
          ps <- ps[-j]
          ps[i] = new_p
          thetas <- thetas[-j, ]
          thetas[i, ] = new_theta
          new_G = list(ps = ps, thetas = thetas)
          Gs = append(Gs, list(new_G))
          # remove j-th row/col from D
          # and recompute distances between new cluster i and all others
          D <- D[-j, -j, drop = FALSE]
          hm = ps[i] * ps[-i] / (ps[i] + ps[-i])
          if (step < n-2) dm = rowSums((thetas[rep(i, n - step - 1), ] - thetas[-i, ])^2) else dm = sum((thetas[i, ] - thetas[-i, ])^2)
          D[i, -i] = hm * dm
          D[-i, i] = hm * dm
          D[i, i] = Inf
     }
     ## last step: only merge and don't need compute dissimilarity matrix
     step = n-1
     idx <- which(D == min(D), arr.ind = TRUE)[1, ]
     i <- idx[1]; j <- idx[2]
     if (i > j) { tmp <- i; i <- j; j <- tmp }
     # put the labels (negative leaves or positive merge-row numbers) into merge[step, ]
     merge[step, ] <- c(cluster_labels[i], cluster_labels[j])
     height[step] <- D[i, j]
     merge_pair[step, ] <- c(i, j)
     # create merged cluster
     new_members <- c(clusters[[i]]$members, clusters[[j]]$members)
     new_p = clusters[[i]]$p + clusters[[j]]$p
     new_theta = (clusters[[i]]$p * clusters[[i]]$theta + clusters[[j]]$p * clusters[[j]]$theta) / new_p
     new_cluster <- list(members=new_members, p=new_p, theta=new_theta)
     clusters[[i]] <- new_cluster
     clusters[[j]] <- NULL
     # this merged cluster is formed at row 'step' of merge
     cluster_labels[i] <- step
     cluster_labels <- cluster_labels[-j]
     ps <- ps[-j]
     ps[i] = new_p
     thetas <- thetas[-j, ]
     thetas[i] = new_theta
     new_G = list(ps = ps, thetas = thetas)
     Gs = append(Gs, list(new_G))
     # recursive order builder: positive integers are merge-row indices (1..n-1), negative are leaves
     get_order <- function(node) {
          if (node < 0) return(-node)
          left <- merge[node, 1]
          right <- merge[node, 2]
          c(get_order(left), get_order(right))
     }
     order <- get_order(n - 1)  # root is the last merge row
     labels <- if (!is.null(rownames(Gs[[1]]$thetas))) rownames(Gs[[1]]$thetas) else as.character(seq_len(n))
     hc = structure(
          list(merge = merge,
               height = height,
               order = order,
               labels = labels,
               call = match.call()),
          class = "hclust"
     )
     return(list(hc=hc, Gs=Gs, merge_pair=merge_pair))
}
