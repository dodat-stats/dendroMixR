# Main function for dendrogram-based mixing

This function performs dendrogram-based hierarchical clustering for
mixture model parameters.

## Usage

``` r
dendrogram_mixing(ps, thetas, sigmas = NULL)
```

## Arguments

- ps:

  A numeric vector of mixture weights.

- thetas:

  A numeric matrix (or list) of component means.

- sigmas:

  Optional: a list or vector of variances/covariances.

## Value

A list containing hierarchical clustering results (`hc`), intermediate
mixing measures (`Gs`), and merge history.
