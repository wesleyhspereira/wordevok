% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{telp_wordevok}
\alias{telp_wordevok}
\title{Armazena a frequência de evocação e a ordem média de evocação de cada evocação
única padronizada.}
\usage{
telp_wordevok(wordevok, groups)
}
\arguments{
\item{wordevok}{Um objeto da classe wordevok.}

\item{groups}{Uma lista contendo os membros da respectiva comunidade.
....}
}
\value{
\code{telp_wordevok} retorna uma lista que contém em cada dimensão sublistas
com um data.frame contendo as evocações únicas, ordem média de evocação e frequência
de evocação por comunidade.
...
}
\description{
....
}
\examples{
data("preferences")
mtx = wordevok_affinity_adjacency(noloop.pref)
g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
add.rownames = "names")
c = cluster_louvain(g, weights = E(g)$weight)
t = telp_wordevok(noloop.pref,c)
}
\author{
Wesley Henrique Silva Pereira
}
