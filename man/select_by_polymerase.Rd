% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OperationsChIA.R
\name{select_by_polymerase}
\alias{select_by_polymerase}
\title{Returns a subset of a ChIA object containing only nodes with the specified 
transcription factors.}
\usage{
select_by_polymerase(chia.obj, polymerase.names)
}
\arguments{
\item{chia.obj}{The ChIA object whose regions must be assessed.}

\item{polymerase.names}{The polymerases whose regions must be returned.}
}
\value{
A subset of the ChIA object containing only nodes with the specified 
  polymerases.
}
\description{
Returns a subset of a ChIA object containing only nodes with the specified 
transcription factors.
}

