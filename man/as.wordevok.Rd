% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordevok.R
\name{as.wordevok}
\alias{as.wordevok}
\title{Creates an wordevok object.}
\usage{
as.wordevok(data, index = NULL, na = NULL)
}
\arguments{
\item{data}{A data.frame class object containing evocations. It may contain a column
identifying the individuals.}

\item{index}{String with the name of the column identifying the individuals.}

\item{na}{specifies how the missing information is informed in the original dataset.
For instance, if “na” is recorded as string (“Missing” or “NA” or “na” or
“null”), this argument should be specified accordingly: \code{na}=“Missing”}
}
\value{
\code{as.wordevok} returns a \code{wordevok} object. This class is based on a
list containing the following slots:

\code{index}: stores the identifier of the evocation vectors and
has the name of the string inserted in the argument \code{index} of the
function \code{as.wordevok}. If \code{index = NULL}, the slot is named "Index".

\code{Dataset}: stores the dataset organized according to the rules
established in \code{details}. Its rows are named according to the assigned
identifier.

\code{Evocations}: stores the vector of standardized single evocations
in the dataset.

\code{Loop}: stores the loops found (same evocation appearing two or
more times in the same row). Loops must be removed before using other
\code{wordevok} functions. Duplicates can be removed by the
\code{removing_loops_wordevok} function.
}
\description{
\code{as.wordevok} creates an object of the wordevok class from a data.frame object.
}
\details{
The \code{wordevok} class was developed to transform evocation data based on the
instrument used in the Free Words Association Technique (FWAT) into relational data
to be used in Social Network Analysis. Since evocations based on the FWAT instrument
are ranked by order of importance and qualitative coding and standardization of these
evocations may produce induced loops, it is likely that loop elimination causes blank
spaces between evocations for a person’s vector of evocations. The \code{wordevok} function
automatically corrects these blank spaces, moving evocations of lower order of
importance to the left of the evocation vector.

The \code{as.wordevok} function provides the option to index each respondent
by an identifier (\code{string} or \code{numeric}). This identifier must occupy a column
in the \code{data.frame}
object and its name must be indicated by a \code{string} in the \code{index} parameter
of the function. For instance: \code{index} = "ID". If \code{index = NULL}, the function automatically indexes
individuals by a \code{numeric} identifier.
}
\examples{
data("preferences")
pref = as.wordevok(preferences,index = "Name")
pref
}
\references{
Abric, J. C. (1994) Las representations sociales: aspects theoriques. In: Abric, J. C. Pratiques sociales et representations.
Paris: Presses Universitaires de France.

Abric, J. C. (2001) A structural approach to social representations.
In: Abric, J. C. Representations of the social: Bridging theoretical traditions.
Malden, UK: Blackwell Publishing.

Pereira, W. H. S. (2017). Representacao da estrutura do pensamento coletivo sobre
as enchentes do Rio Doce: conectando indivíduos afins através da teoria dos grafos.
Monografia (Graduacao em Estatistica) - Universidade Fedral de Minas Gerais.
}
\author{
Wesley H. S. Pereira (ICEx/UFMG), Denise Duarte (ICEx/UFMG), Rodrigo B. Ribeiro (IMPA), Gilvan Guedes (CEDEPLAR/UFMG).
}
