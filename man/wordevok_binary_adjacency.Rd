% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{wordevok_binary_adjacency}
\alias{wordevok_binary_adjacency}
\title{Construction of the binary adjacency matrix.}
\usage{
wordevok_binary_adjacency(wordevok)
}
\arguments{
\item{wordevok}{A wordevok class object.}
}
\value{
\code{wordevok_intersection_adjacency} returns the adjacency matrix
in which
each line and each column are one indiviual
and each cell is the sum of all weigths assigned to that evocation link.
}
\description{
This function constructs the binary adjacency matrix.
}
\examples{
data("preferences")
mtx = wordevok_binary_adjacency(noloop.pref)
mtx
}
\author{
Wesley Henrique Silva Pereira
}
