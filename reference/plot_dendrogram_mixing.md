# Plot dendrogram of mixing measures

Plot dendrogram of mixing measures

## Usage

``` r
plot_dendrogram_mixing(
  dmm,
  dim = 1,
  point_size = 3,
  line_width = 0.5,
  palette = NULL,
  show_legend = FALSE,
  main = NULL
)
```

## Arguments

- dmm:

  Output from `dendrogram_mixing`.

- dim:

  Dimension to plot.

- point_size:

  Point size.

- line_width:

  Line width.

- palette:

  Vector of colors.

- show_legend:

  Show colors associated with each step in the dendrogram.

- main:

  Title name

## Value

A base R plot.
