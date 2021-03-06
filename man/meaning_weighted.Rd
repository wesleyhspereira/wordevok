% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{meaning_weighted}
\alias{meaning_weighted}
\title{Measuring weights  to edges of personal network.}
\usage{
meaning_weighted(n, N, output = "Matrix")
}
\arguments{
\item{n}{Number of evocations expressed by a individual.}

\item{N}{Maximum number of evocations expressed allowed to an individual.}

\item{output}{Type of output to be generated. See details to more informations.}
}
\value{
\code{meaning_weighted} returns a set of numeric weights that will be a
ssigned to each evocation pair
of the individual vector.
}
\description{
Under Free Word Association Technique, this function measures the weights
assigned to each edge pair linking different evocations of the same individual.
}
\details{
Under the criteria of the Free Word Association Technique, we suppose that that
the evocations that are expressed by the same individual are mentally connected.
We propose a strategy to measure the weight of the connection for each pair of
these mental connections.
Details on the meaning_weighted can be found in
\url{http://wesleyhenriquesp.wixsite.com/rwordevok/meaning-weighted}.
}
\examples{
Murillo = c("Regression Analysis","Multivariate Statistics",
"General Statistiscs", "Experiment Planning", "Sampling")
n = length(Murillo)
N = n
mw = meaning_weighted(n,N)
colnames(mw) = Murillo
rownames(mw) = Murillo
mw
}
\author{
Wesley Henrique Silva Pereira
}
