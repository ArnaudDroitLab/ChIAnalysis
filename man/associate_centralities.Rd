% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotateChIA.R
\name{associate_centralities}
\alias{associate_centralities}
\title{Associates centrality scores and boolean centrality markers to regions.}
\usage{
associate_centralities(chia.obj, which.measures = c("Degree", "Betweenness",
  "Eigenvector", "Closeness"), weight.attr = NULL)
}
\arguments{
\item{chia.obj}{ChIA-PET data, as returned by \code{\link{annotate_chia}}.}

\item{which.measures}{A vector containing the names of the measures to be used
to assess centrality. Those can be "Degree", "Betweenness" and "Eigenvector".}

\item{weight.attr}{The anme of the edge attribute to be sued as a weight in 
centrality calculations.}
}
\value{
The annotated chia.obj.
}
\description{
Associates centrality scores and boolean centrality markers to regions.
}

