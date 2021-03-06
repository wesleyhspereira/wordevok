% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{telp_wordevok_plot}
\alias{telp_wordevok_plot}
\title{Apresenta o gráfico de quadrantes de Vergés.}
\usage{
telp_wordevok_plot(twfile, ln = FALSE)
}
\arguments{
\item{twfile}{Uma dimensão de um objeto da classe \code{Multiple-TELP-wordevok}, gerado pela função \code{telp_wordevok}}

\item{ln}{Lógico. \code{TRUE} implica na utilização da escala logarítmica nos valores encontrados.
...}
}
\description{
Esta função apresenta o gráfico de quadrantes de Vergés.
}
\examples{
data("preferences")
mtx = wordevok_affinity_adjacency(noloop.pref)
g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
add.rownames = "names")
c = cluster_louvain(g, weights = E(g)$weight)
t = telp_wordevok(noloop.pref,c)
telp_wordevok_plot(t[[1]])
}
\author{
Wesley Henrique Silva Pereira
}
