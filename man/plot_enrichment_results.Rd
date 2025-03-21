% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_enrichment_results.R
\name{plot_enrichment_results}
\alias{plot_enrichment_results}
\title{A plot function for TEA.}
\usage{
plot_enrichment_results(df, labels_inside_columns = NULL, max_bars = 15L)
}
\arguments{
\item{df}{dataframe as output by `enrichment_analysis()`.}

\item{labels_inside_columns}{write the labels on top of the bars, instead of on the side. Leave NULL for automatic selection.}

\item{max_bars}{maximal number of bars to be shown.}
}
\value{
a ggplot object
}
\description{
A plot function for TEA.
}
