#' Main function for dendrogram-based mixing
#'
#' This function performs dendrogram-based hierarchical clustering
#' for mixture model parameters.
#'
#' @param ps A numeric vector of mixture weights.
#' @param thetas A numeric matrix (or list) of component means.
#' @param sigmas Optional: a list or vector of variances/covariances.
#'
#' @return A list containing hierarchical clustering results (`hc`),
#'         intermediate mixing measures (`Gs`), and merge history.
#' @export
dendrogram_mixing <- function(ps, thetas, sigmas=NULL) {
     thetas <- if (is.list(thetas)) {
          do.call(rbind, thetas)
     } else {
          as.matrix(thetas)
     }

     if (nrow(thetas)!=length(ps)){
          print("Number of atoms theta must equal number of probability ps")
     }

     if(is.null(sigmas)==TRUE){
          if (ncol(thetas)==1){
               ## 1d strong identifiability
               return(dendrogram_mixing_1d_strong(ps, thetas[, 1]))
          } else {
               ## multi-dimensional strong identifiability
               return(dendrogram_mixing_multid_strong(ps, thetas))
          }
     }

     if(is.null(sigmas)==FALSE){
          if (ncol(thetas)==1){
               ## 1d location-scale Gaussian (weak identifiable)
               return(dendrogram_mixing_1d_weak(ps, thetas[, 1], sigmas))
          } else {
               ## multi-dimensional location-scale Gaussian (weak identifiable)
               return(dendrogram_mixing_multid_weak(ps, thetas, sigmas))
          }
     }
}
