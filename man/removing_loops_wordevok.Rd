% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{removing_loops_wordevok}
\alias{removing_loops_wordevok}
\title{Removal of loops of a wordevok object.}
\usage{
removing_loops_wordevok(wordevok)
}
\arguments{
\item{wordevok}{A wordevok class object.}
}
\value{
\code{removing_loops_wordevok} returns a wordevok object without loops.
}
\description{
This function removes the loops from an wordevok object.
}
\details{
Under the Free Word Association Technique, original evocations can vary from
single words to complete sentences. However, the analysis of this original
evocation set is challenging, since evocations are sometimes too specific.
Moreover, although certain evocations are syntactically distinct, they can be
semantically equivalent. In this case, qualitative standardization of the original
evocations is used, giving rise to unique standardized evocations. These evocations
can be found in the slot Evocations of wordevok object.

The standardization process creates unnatural loops in the evocation dataset:
the appearance of the same evocation twice for the same individual. In Graph Theory,
there are applications where such loops are acceptable. For the analysis of collective
meanings, however, this behavior does not make sense and represents a looped induced by
the standardization process. The best alternative to this problem is to manually try to
reclassify one of the repeated evocations in a different semantic class that agrees with
its original meaning. When manual reclassification is not possible, the
removing_loops_wordevok function removes from the wordevok object the repeated
evocation that has been assigned to the lower order of importance. In addition,
the function automatically rearranges the remaining evocations so that the NA
observations are placed to the rightest positions of the vector.
}
\examples{
#Creating a wordevok object:
data("preferences")
pref = as.wordevok(preferences,index = "Name")
pref
#Removing loops:
noloop.pref = removing_loops_wordevok(pref)
noloop.pref
}
\references{
Pereira, W. H. S. (2017). Representacao da estrutura do pensamento coletivo sobre
as enchentes do Rio Doce: conectando indivíduos afins através da teoria dos grafos.
Monografia (Graduacao em Estatistica) - Universidade Fedral de Minas Gerais.
}
\author{
Wesley Henrique Silva Pereira
}
