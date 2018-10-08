#' Creates an wordevok object.
#'
#' \code{as.wordevok} creates an object of the wordevok class from a data.frame object.
#'@param data A data.frame class object containing evocations. It may contain a column
#'identifying the individuals.
#'@param index String with the name of the column identifying the individuals.
#'@param na specifies how the missing information is informed in the original dataset.
#'For instance, if “na” is recorded as string (“Missing” or “NA” or “na” or
#'“null”), this argument should be specified accordingly: \code{na}=“Missing”
#'
#'@return \code{as.wordevok} returns a \code{wordevok} object. This class is based on a
#'list containing the following slots:
#'
#'\code{index}: stores the identifier of the evocation vectors and
#' has the name of the string inserted in the argument \code{index} of the
#' function \code{as.wordevok}. If \code{index = NULL}, the slot is named "Index".
#'
#'\code{Dataset}: stores the dataset organized according to the rules
#'established in \code{details}. Its rows are named according to the assigned
#'identifier.
#'
#'\code{Evocations}: stores the vector of standardized single evocations
#'in the dataset.
#'
#'\code{Loop}: stores the loops found (same evocation appearing two or
#'more times in the same row). Loops must be removed before using other
#'\code{wordevok} functions. Duplicates can be removed by the
#'\code{removing_loops_wordevok} function.
#'
#'@author Wesley H. S. Pereira (ICEx/UFMG), Denise Duarte (ICEx/UFMG), Rodrigo B. Ribeiro (IMPA), Gilvan Guedes (CEDEPLAR/UFMG).
#'@references Abric, J. C. (1994) Las representations sociales: aspects theoriques. In: Abric, J. C. Pratiques sociales et representations.
#' Paris: Presses Universitaires de France.
#'@references Abric, J. C. (2001) A structural approach to social representations.
#'In: Abric, J. C. Representations of the social: Bridging theoretical traditions.
#' Malden, UK: Blackwell Publishing.
#' @references Pereira, W. H. S. (2017). Representacao da estrutura do pensamento coletivo sobre
#' as enchentes do Rio Doce: conectando indivíduos afins através da teoria dos grafos.
#' Monografia (Graduacao em Estatistica) - Universidade Fedral de Minas Gerais.

#'@details
#'The \code{wordevok} class was developed to transform evocation data based on the
#'instrument used in the Free Words Association Technique (FWAT) into relational data
#'to be used in Social Network Analysis. Since evocations based on the FWAT instrument
#'are ranked by order of importance and qualitative coding and standardization of these
#'evocations may produce induced loops, it is likely that loop elimination causes blank
#'spaces between evocations for a person’s vector of evocations. The \code{wordevok} function
#'automatically corrects these blank spaces, moving evocations of lower order of
#'importance to the left of the evocation vector.
#'
#'The \code{as.wordevok} function provides the option to index each respondent
#'by an identifier (\code{string} or \code{numeric}). This identifier must occupy a column
#'in the \code{data.frame}
#'object and its name must be indicated by a \code{string} in the \code{index} parameter
#'of the function. For instance: \code{index} = "ID". If \code{index = NULL}, the function automatically indexes
#'individuals by a \code{numeric} identifier.
#'@examples
#'data("preferences")
#'pref = as.wordevok(preferences,index = "Name")
#'pref
#'@export

as.wordevok = function(data, index = NULL,na = NULL)
{
  if(!is.null(na))
  {
    data = apply(data,2,as.character)
    for(i in 1:dim(data)[2])
    {
      data[data[,i] == na,i] = NA
    }
    data = data.frame(data)
  }
  if(class(data)[1] != "data.frame")
  {
    stop ("Only objects of the data.frame class as supported by this function.")
  }
  if (!is.null(index)){
    idx = as.character(data[,match(index,colnames(data))])
    data = data[,-match(index,colnames(data))]
    name = colnames(data)
    data = data.frame(idx,data)
    colnames(data) = c(index,name)
  } else {
    index = "Index"
    idx = rep(1:dim(data)[1])
    name = colnames(data)
    data = data.frame(idx,data)
    colnames(data) = c("Index",name)
  }
  name = colnames(data)
  idx = as.character(data[,1])
  evocation = data[,-1]
  base = apply(evocation,2,as.character)
  if(dim(evocation)[1] == 1)
  {
    base = t(base)
  }
  for (i in 1:dim(evocation)[1])
  {
    exchange = which(is.na(base[i,]) == TRUE)
    if(length(exchange)!=0)
    {
      base[i,] = c(base[i,-exchange],base[i,exchange])
    }
  }
  out = matrix(ncol = 1+dim(evocation)[2],nrow = dim(evocation)[1])
  out[,1] = idx
  out[,2:(1+dim(evocation)[2])] = base
  colnames(out) = name
  out = as.data.frame(out)
  out = out[!is.na(out[,2]),]
  output = list()
  output[[index]] = out[,1]
  rownames(out) = out[,1]
  out = out[,-1]
  base = apply(out,2,as.character)
  if(dim(evocation)[1] == 1)
  {
    base = t(base)
  }
  rownames(base) = rownames(out)
  colnames(base) = colnames(out)
  evocations = list()
  for(i in 1:dim(base)[2])
  {
    evocations[[i]] = (base[,i])
  }
  evocations = Reduce(c,evocations)
  evocations = unique(evocations[!is.na(evocations)])
  loop = data.frame(NA,NA,NA,NA)
  colnames(loop) = c(index,"Evocation","Source","Target")
  iter = 0
  for(i in 1:dim(base)[1])
  {
    for(j in 1:(dim(base)[2]-1))
    {
      for(k in (j+1):(dim(base)[2]))
      {
        if(!is.na(base[i,j]) && !is.na(base[i,k]))
        {
          if(base[i,j] == base[i,k])
          {
            iter = iter + 1
            loop[iter,1] = rownames(base)[i]
            loop[iter,2] = base[i,j]
            loop[iter,3] = colnames(base)[j]
            loop[iter,4] = colnames(base)[k]
          }
        }
      }
    }
  }
  output$Dataset = out
  output$Evocations = sort(evocations)
  if(iter == 0)
  {
    output$Loop = NULL
  } else {
    output$Loop = as.data.frame(loop)
    warning("The wordevok file contains loops.")
  }
  class(output) = "wordevok"
  return(output)
}

#' Removal of loops of a wordevok object.
#'
#' This function removes the loops from an wordevok object.
#'@param wordevok A wordevok class object.
#'@return \code{removing_loops_wordevok} returns a wordevok object without loops.
#'@author Wesley Henrique Silva Pereira
#'@details
#'Under the Free Word Association Technique, original evocations can vary from
#'single words to complete sentences. However, the analysis of this original
#'evocation set is challenging, since evocations are sometimes too specific.
#'Moreover, although certain evocations are syntactically distinct, they can be
#'semantically equivalent. In this case, qualitative standardization of the original
#'evocations is used, giving rise to unique standardized evocations. These evocations
#'can be found in the slot Evocations of wordevok object.
#'
#'The standardization process creates unnatural loops in the evocation dataset:
#'the appearance of the same evocation twice for the same individual. In Graph Theory,
#'there are applications where such loops are acceptable. For the analysis of collective
#'meanings, however, this behavior does not make sense and represents a looped induced by
#'the standardization process. The best alternative to this problem is to manually try to
#'reclassify one of the repeated evocations in a different semantic class that agrees with
#'its original meaning. When manual reclassification is not possible, the
#'removing_loops_wordevok function removes from the wordevok object the repeated
#'evocation that has been assigned to the lower order of importance. In addition,
#'the function automatically rearranges the remaining evocations so that the NA
#'observations are placed to the rightest positions of the vector.
#'
#'@references Pereira, W. H. S. (2017). Representacao da estrutura do pensamento coletivo sobre
#' as enchentes do Rio Doce: conectando indivíduos afins através da teoria dos grafos.
#' Monografia (Graduacao em Estatistica) - Universidade Fedral de Minas Gerais.
#'@examples
#'#Creating a wordevok object:
#'data("preferences")
#'pref = as.wordevok(preferences,index = "Name")
#'pref
#'#Removing loops:
#'noloop.pref = removing_loops_wordevok(pref)
#'noloop.pref
#'@export

removing_loops_wordevok = function(wordevok)
{
  if(class(wordevok)[1] != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(is.null(wordevok$Loop))
  {
    return(wordevok)
  } else {
    loop = wordevok$Loop
    dataset = wordevok$Dataset
    base = apply(dataset,2,as.character)
    colnames(base) = colnames(dataset)
    rownames(base) = rownames(dataset)
    for (i in 1:dim(loop)[1])
    {
      base[match(loop[i,1],rownames(base)),loop[i,4]] = NA
    }
  }
  wordevok$Loop = NULL
  loop = data.frame(NA,NA,NA,NA)
  colnames(loop) = c(names(wordevok)[1],"Evocation","Source","Target")
  iter = 0
  for(i in 1:dim(base)[1])
  {
    for(j in 1:(dim(base)[2]-1))
    {
      for(k in (j+1):(dim(base)[2]))
      {
        if(!is.na(base[i,j]) && !is.na(base[i,k]))
        {
          if(base[i,j] == base[i,k])
          {
            iter = iter + 1
            loop[iter,1] = rownames(base)[i]
            loop[iter,2] = base[i,j]
            loop[iter,3] = colnames(base)[j]
            loop[iter,4] = colnames(base)[k]
          }
        }
      }
    }
  }
  if(iter == 0)
  {
    wordevok$Loop = NULL
  } else {
    warning("This dataset contains loops")
    wordevok$Loop = as.data.frame(loop)
  }
  for (i in 1:dim(base)[1])
  {
    exchange = which(is.na(base[i,]) == TRUE)
    if(length(exchange)!=0)
    {
      base[i,] = c(base[i,-exchange],base[i,exchange])
    }
  }
  wordevok$Dataset = as.data.frame(base)
  return(wordevok)
}

#' Measuring affinity within two individuals.
#'
#' Under Free Word Association Technique, this function measures the affinity
#' between the evoked vectors of two individuals.
#'@param subject1 Evoked vector of the individual 1.
#'@param subject2 Evoked vector of the individual 2.
#'@return \code{affinity} returns a numeric measure of affinity between the vectors
#'inputed.
#'@author Wesley H. S. Pereira (ICEx/UFMG), Denise Duarte (ICEx/UFMG), Rodrigo B. Ribeiro (IMPA), Gilvan Guedes (CEDEPLAR/UFMG).
#'@details
#'
#'Under the criteria of the Free Word Association Technique, we propose a
#'coefficient to calculate the affinity (similiraty) between two vectors of evocations.
#'Details on the affinity coefficient can be found in
#'\url{http://wesleyhenriquesp.wixsite.com/rwordevok/affinity}.
#'
#'The maximum number of evocations allowed to a vector is controlled by the length
#'of the vectors. Therefore, the vectors \code{subject1} and \code{subject2}
#'must have equal size. If an individual pronounces less than the maximum \code{N},
#'the vector must be structured so that n evocations pronounced are in the first \code{n}
#'positions of the vector according to their importance to the individual and the
#'remaining \code{N - n} positions must be filled with NA observations.
#'
#'@references Pereira, W. H. S. (2017).
#'@examples
#'#Creating evocation's vectors:
#'Murillo = c("Regression Analysis","Multivariate Statistics",
#'"General Statistiscs", "Experiment Planning", "Sampling")
#'Ingrid = c("Regression Analysis","Multivariate Statistics",
#'"Temporal Series","Statistical Quality Control", NA)
#'affinity(Murillo,Ingrid)
#'@export

affinity = function(subject1,subject2)
{
  if(length(subject1) != length(subject2))
  {
    stop("Arguments must be the same size.")
  }
  n1 = length(subject1[!is.na(subject1)])
  n2 = length(subject2[!is.na(subject2)])
  n = max(n1,n2)
  N = length(subject1)
  denominator = n*n*(n+1)
  cor1 = (N - n) * (N - n + 1)
  cor2 = N * (N + 1)
  omega_n = 1 - cor1/cor2
  affinity_coefficient = 0
  for (i in 1:n)
  {
    for(j in 1:n)
    {
      if(!is.na(subject1[i]) && !is.na(subject2[j]))
      {
        if(subject1[i] == subject2[j])
        {
          theta_i_j = 2*(n + 1) - (i + j)
          rho_i_j = n - abs(i-j)
          affinity_coefficient = affinity_coefficient + theta_i_j*rho_i_j
        }
      }
    }
  }
  beta_uv = affinity_coefficient/denominator
  affinity_coefficient = beta_uv*omega_n
  return(affinity_coefficient)
}

#' Measuring weights  to edges of personal network.
#'
#' Under Free Word Association Technique, this function measures the weights
#' assigned to each edge pair linking different evocations of the same individual.
#'@param n Number of evocations expressed by a individual.
#'@param N Maximum number of evocations expressed allowed to an individual.
#'@param output Type of output to be generated. See details to more informations.
#'@details
#'Under the criteria of the Free Word Association Technique, we suppose that that
#'the evocations that are expressed by the same individual are mentally connected.
#'We propose a strategy to measure the weight of the connection for each pair of
#' these mental connections.
#'Details on the meaning_weighted can be found in
#'\url{http://wesleyhenriquesp.wixsite.com/rwordevok/meaning-weighted}.
#'@author Wesley Henrique Silva Pereira
#'@return \code{meaning_weighted} returns a set of numeric weights that will be a
#'ssigned to each evocation pair
#'of the individual vector.
#'@examples
#'Murillo = c("Regression Analysis","Multivariate Statistics",
#'"General Statistiscs", "Experiment Planning", "Sampling")
#'n = length(Murillo)
#'N = n
#'mw = meaning_weighted(n,N)
#'colnames(mw) = Murillo
#'rownames(mw) = Murillo
#'mw
#'@export

meaning_weighted = function(n, N , output = "Matrix")
{
  iter = 0
  numerator = NULL
  denominator = 0
  if (output != "Matrix" & output != "Vector")
  {
    stop("This function is defined only to Vector or Matrix output.")
  }
  if(n < 2 || N < 2)
  {
    stop("This function is defined only to n and N greater than 2.")
  }
  if(n > N)
  {
    stop("This function is defined only to n lesser or equal than N.")
  }
  cor1 = (N - n)*(N - n + 1)
  cor2 = N*(N + 1)
  omega_n = 1 - cor1/cor2
  if (output == "Matrix")
  {
    numerator = matrix(ncol = N, nrow = N)
    colnames(numerator) = paste0("Evocation_",01:N)
    rownames(numerator) = paste0("Evocation_",01:N)
    numerator[,] = NA

    for (i in 1:(n-1))
    {
      for (j in (i+1):n){
        numerator[i,j] =  ((2^(2*n+2-i-j))-1)*(n - abs(i - j))
        denominator = denominator + numerator[i,j]
      }
    }
    meaning = (numerator/denominator)*omega_n
  }
  if (output == "Vector")
  {
    for (i in 1:(n-1))
    {
      for (j in (i+1):n){
        iter = iter + 1
        numerator[iter] =  ((2^(2*n+2-i-j))-1)*(n - abs(i - j))
        denominator = denominator + numerator[iter]
      }
    }
    meaning = (numerator/denominator)*omega_n
  }
  return(meaning)
}

#' Construction of the adjacency matrix by intersection.
#'
#' This function constructs the adjacency matrix by the number of
#' evocations in common.
#'@param wordevok A wordevok class object.
#'@author Wesley Henrique Silva Pereira
#'@return \code{wordevok_intersection_adjacency} returns the adjacency matrix in which
#'each line and each column are one
#'indiviual and each cell is the common evocation's number between the line's
#'individual and column's individual.
#'@examples
#'data("preferences")
#'mtx = wordevok_intersection_adjacency(noloop.pref)
#'mtx
#'@export

wordevok_intersection_adjacency = function(wordevok)
{
  if(class(wordevok) != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(!is.null(wordevok$Loop))
  {
    stop("To use this function, the wordevok object cannot contain loops.")
  }
  adjacency = matrix(nrow = dim(wordevok$Dataset)[1], ncol = dim(wordevok$Dataset)[1])
  colnames(adjacency) = as.character(wordevok[[1]])
  rownames(adjacency) = colnames(adjacency)
  diag(adjacency) = 0
  for(i in 1:(dim(wordevok$Dataset)[1]-1))
  {
    subject1 = wordevok$Dataset[i,]
    subject.1 = subject1[!is.na(subject1)]
    for(j in (i+1):(dim(wordevok$Dataset)[1]))
    {
      subject2 = wordevok$Dataset[j,]
      subject.2 = subject2[!is.na(subject2)]
      adjacency[i,j] = sum(subject.1 %in% subject.2)
      adjacency[j,i] = sum(subject.1 %in% subject.2)
    }
  }
  return(adjacency)
}

#' Construction of the adjacency matrix by affinity coefficient.
#'
#' This function constructs the adjacency matrix by affinity coefficient.
#'@param wordevok A object of wordevok class.
#'@return \code{wordevok_affinity_adjacency} returns the adjacency matrix in which each
#'line and each column are one indiviual
#' and each cell is the affinity coefficient between the line's individual and
#' column's individual.
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'mtx
#'@export

wordevok_affinity_adjacency = function(wordevok)
{
  if(class(wordevok) != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(!is.null(wordevok$Loop))
  {
    stop("To use this function, the wordevok object cannot contain loops.")
  }
  adjacency = matrix(nrow = dim(wordevok$Dataset)[1], ncol = dim(wordevok$Dataset)[1])
  colnames(adjacency) = as.character(wordevok[[1]])
  rownames(adjacency) = colnames(adjacency)
  adjacency[,] = 0
  for(i in 1:(dim(wordevok$Dataset)[1]-1))
  {
    subject1 = apply(wordevok$Dataset[i,],2,as.character)
    subject.1 = subject1[!is.na(subject1)]
    for(j in (i+1):(dim(wordevok$Dataset)[1]))
    {
      subject2 = apply(wordevok$Dataset[j,],2,as.character)
      subject.2 = subject2[!is.na(subject2)]
      if(sum(subject.1 %in% subject.2) > 0)
      {
        adjacency[i,j] = affinity(subject1,subject2)
        adjacency[j,i] = adjacency[i,j]
      }
    }
  }
  return(adjacency)
}

#' Construction of the binary adjacency matrix.
#'
#' This function constructs the binary adjacency matrix.
#'@param wordevok A wordevok class object.
#'@return \code{wordevok_intersection_adjacency} returns the adjacency matrix
#'in which
#'each line and each column are one indiviual
#'and each cell is the sum of all weigths assigned to that evocation link.
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_binary_adjacency(noloop.pref)
#'mtx
#'@export


wordevok_binary_adjacency = function(wordevok)
{
  if(class(wordevok) != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(!is.null(wordevok$Loop))
  {
    stop("To use this function, the wordevok object cannot contain loops.")
  }
  adjacency = matrix(nrow = dim(wordevok$Dataset)[1], ncol = dim(wordevok$Dataset)[1])
  colnames(adjacency) = as.character(wordevok[[1]])
  rownames(adjacency) = colnames(adjacency)
  adjacency[,] = 0
  for(i in 1:(dim(wordevok$Dataset)[1]-1))
  {
    subject1 = wordevok$Dataset[i,]
    subject.1 = subject1[!is.na(subject1)]
    for(j in (i+1):(dim(wordevok$Dataset)[1]))
    {
      subject2 = wordevok$Dataset[j,]
      subject.2 = subject2[!is.na(subject2)]
      if(sum(subject.1 %in% subject.2) > 0)
      {
        adjacency[i,j] = 1
        adjacency[j,i] = 1
      }
    }
  }
  return(adjacency)
}

#' Construction of the adjacency matrix by meaning weighted method.
#'
#' This function constructs the adjacency matrix by meaning weighted method.
#'@param wordevok A wordevok class object.
#'@author Wesley Henrique Silva Pereira
#'@return \code{wordevok_meaning_adjacency} returns the adjacency matrix in which each
#'line and each column are an evocation
#'and each cell takes on 1 if the line's individual and column's individual have
#'common evocations. Otherwise, the cells takes on 0.
#'@examples
#'data("preferences")
#'mtx = wordevok_meaning_adjacency(noloop.pref)
#'mtx
#'@export


wordevok_meaning_adjacency = function(wordevok)
{
  if(class(wordevok)[1] != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(!is.null(wordevok$Loop))
  {
    stop("To use this function, the wordevok object cannot contain loops.")
  }
  evocations = wordevok$Evocations
  meaning = wordevok_meaning_list(wordevok,"Simplify")
  adjacency = matrix(ncol = length(evocations),nrow = length(evocations))
  colnames(adjacency) = evocations
  rownames(adjacency) = evocations
  adjacency[,] = 0
  for(i in 1:dim(meaning)[1])
  {
    source = as.character(meaning[i,1])
    target = as.character(meaning[i,2])
    adjacency[source,target] = meaning[i,3]
    adjacency[target,source] = meaning[i,3]
  }
  return(adjacency)
}

#' Construction of the list of edges by meaning weighted method.
#'
#' This function constructs the list of edges by meaning weighted method.
#'@param wordevok A wordevok class object.
#'@param output Type of output to be generated. See details to more information.
#'@author Wesley Henrique Silva Pereira
#'@return \code{wordevok_meaning_list} returns a list containing the edges and their
#'respective weights.
#'@details
#'The output argument of the function x must be:
#'
#'- "Complete": the object contains the identifier of the individual
#'generating the edge, the vertices bordering the edge and the weights
#'assigned to the edges.
#'
#'- "Edgelist": the same elements of the "Complete" output, except the
#'identifier of the individual generator.
#'
#'- "Simplify": this output contains only the unique edges and their weights.
#'The weights associated with these represent the sum of the
#'weights of all your duplicates.
#'@examples
#'data("preferences")
#'lst = wordevok_meaning_list(noloop.pref)
#'lst
#'@export


wordevok_meaning_list = function(wordevok, output = "Complete")
{
  if(class(wordevok)[1] != "wordevok")
  {
    stop ("Only objects of the wordevok class as supported by this function.")
  }
  if(!is.null(wordevok$Loop))
  {
    stop("To use this function, the wordevok object cannot contain loops.")
  }
  if(output != "Complete" && output != "Edgelist" && output != "Simplify")
  {
    stop("This output type is not defined.")
  }
  dataset = apply(wordevok$Dataset,2,as.character)
  n = dim(dataset)[2]
  source = list()
  for (i in 1:n)
  {
    source[[i]] = as.data.frame(cbind(as.character(wordevok[[1]]),dataset[,i],i))
  }
  target = list()
  c = 0
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      c = c + 1
      target[[c]] = as.data.frame(cbind(source[[i]],source[[j]][,-1]))
    }
  }
  listing  = as.data.frame(rbindlist(target))
  listing  = as.data.frame(apply(listing,2,as.character))
  listing  = na.omit(listing)
  colnames(listing) = c(names(wordevok)[1],"Source","Source_ID","Target","Target_ID")
  listing = as.data.frame(apply(listing,2,as.character))
  ponder = split(listing,listing[,1])
  for (i in 1:length(ponder))
  {
    conect = max(as.numeric(as.character(ponder[[i]]$Target_ID)))
    ponder[[i]]$weight = meaning_weighted(conect,n,"Vector")
  }
  listing  = as.data.frame(rbindlist(ponder))
  if(output == "Complete")
  {
    return(listing[,c(names(wordevok)[1],"Source","Target","weight")])
  } else {
      listing = listing[,c("Source","Target","weight")]
      if(output == "Edgelist")
      {
        return(listing)
      } else {
          combining = t(combn(wordevok$Evocations,2))
          weigh = NULL
          for(i in 1:dim(combining)[1]){
            verse = sum(listing[listing$Source == combining[i,1]  &
                       listing$Target == combining[i,2],"weight"])
            inverse = sum(listing[listing$Source == combining[i,2] &
                        listing$Target == combining[i,1],"weight"])
            weigh[i] = verse + inverse
          }
          listing = as.data.frame(combining)
          listing$v3 = weigh
          listing = listing[listing$v3 > 0,]
          colnames(listing) = c("Source","Target","weight")
          return(listing)
      }
  }
}

#' Construction of the adjacency matrix by laplacian method.
#'
#' This function constructs the adjacency matrix by laplacian method.
#'@param wordevok A wordevok class object.
#'@param normalized Logical argument. If \code{normalized = TRUE}, the output will be
#'normalized.
#'@return \code{wordevok_laplacian_adjacency} returns the adjacency matrix in which each
#'line and each column are an
#'evocation and each cell takes on -1 if the line's individual and column's
#'individual have common evocations. The diagonal of the matrix equals the number of
#'edges incident on the vertex.
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_laplacian_adjacency(noloop.pref)
#'mtx
#'@export

wordevok_laplacian_adjacency = function(wordevok,normalized = FALSE)
{
  laplacian = wordevok_binary_adjacency(wordevok)
  dgn = as.numeric(rowSums(laplacian))
  laplacian = laplacian*(-1)
  for (i in 1:length(dgn))
  {
    laplacian[i,i] = dgn[i]
  }
  if (normalized)
  {
    for (i in 1:length(dgn))
    {
      if(dgn[i] != 0)
      {
        laplacian[i,] = laplacian[i,]/dgn[i]
      }
    }
    return(laplacian)
  } else {
    return(laplacian)
  }
}

#' Extrai uma lista contendo as evocações das comunidades do objeto wordevok.
#'
#' Essa função extrai, mediante a informação dos membros de cada comunidade
#' encontrada por um método de agrupamento em redes, uma lista contendo os
#' subconjuntos de evocações de cada uma das comunidades informadas.
#'@param wordevok Um objeto da classe wordevok.
#'@param groups Uma lista contendo os membros da respectiva comunidade.
#'....
#'@return \code{wordevok_comm_subsets} Retorna uma lista onde em cada uma das
#'dimensões são alocados os subconjuntos do objeto \code{wordevok} orginal
#'de acordo com as comunidades informadas.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#' add.rownames = "names")
#'c = cluster_louvain(g,weights = E(g)$weight)
#'sub = wordevok_comm_subsets(noloop.pref,c)
#'@export

wordevok_comm_subsets = function(wordevok,groups)
{
  n = length(groups)
  subsets = list()
  for (i in 1:n)
  {
    subsets[[i]] = wordevok$Dataset[match(groups[[i]],rownames(wordevok$Dataset)),]
  }
  names(subsets) = paste0("Group",1:n)
  return(subsets)
}

#' Extract the submatrix of adjacency matrix under founded communities.
#'
#' ....
#'@param adjacency A adjacency matrix.
#'@param groups Uma lista contendo os membros da respectiva comunidade.
#'@param laplacian ...
#'....
#'@return \code{wordevok_comm_submatrix} ...
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g,weights = E(g)$weight)
#'sub = wordevok_comm_submatrix(mtx,c,laplacian = FALSE)
#'@export

wordevok_comm_submatrix = function(adjacency,groups,laplacian = FALSE)
{
  n = length(groups)
  name = rownames(adjacency)
  submatrix = list()
  for (i in 1:n)
  {
    idx = match(groups[[i]],name)
    if (length(idx) > 1)
    {
      mtx = adjacency[idx,]
      mtx = mtx[,idx]
      submatrix[[i]] = mtx
    } else {
      mtx = adjacency[idx,]
      mtx = mtx[idx]
      submatrix[[i]] = mtx
    }
  }
  names(submatrix) = paste0("group",1:n)
  if(laplacian == FALSE)
  {
    return(submatrix)
  } else {
    for (i in 1:n)
    {
      mtx = submatrix[[i]]
      s = colSums(mtx)
      for (j in 1:dim(mtx)[1])
      {
        mtx[j,j] = mtx[j,j] - s[j]
      }
      submatrix[[i]] = mtx
    }
    return(submatrix)
  }
}

#' Transforma uma matriz de adjacências em uma lista de arestas.
#'
#' Esta função transforma a matriz de adjacências em uma lista de arestas.
#'@param adjacency Uma matriz de adjacencias.
#'...
#'@return \code{wordevok_adjacency_to_list} retorna uma lista de adjacências
#'conforme fornecida pela matriz de adjacências original. É importante ressaltar
#'que essa transformação ignora as limitações de se representar as conexões do
#'grafo através de uma lista de adjacência (como por exemplo a incapacidade de
#'represnetar vértices isolados). As modificações necessárias ficam
#'a cargo do usuário.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'lst = wordevok_adjacency_to_list(mtx)
#'@export

wordevok_adjacency_to_list = function(adjacency)
{
  name = rownames(adjacency)
  combination = as.data.frame(t(combn(name,2)))
  colnames(combination) = c("Source","Target")
  combination$weight = 0
  for (i in 1:dim(combination)[1])
  {
    combination$weight[i] = adjacency[as.character(combination$Source[i]),
                                      as.character(combination$Target[i])]
  }
  return(combination)
}

#' Cria subgrupos de uma lista de acordo com os grupos informados.
#'
#' Esta função cria sublistas de adjacência mediante à informação das
#' comunidades.
#'
#'@param list Uma lista de arestas.
#'@param groups Uma lista contendo os membros da respectiva comunidade.
#'....
#'@return \code{split_list_wordevok} retorna uma lista contendo em cada dimensão
#'as ligações por comunidade. Assim, a função ignora as ligações entre comunidades, deixando apenas as
#' ligações entre membros de uma mesma comunidade.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_meaning_list(noloop.pref)
#'g = graph_from_data_frame(mtx[,c(2,3)], directed = FALSE)
#'c = cluster_louvain(g)
#'lst = split_list_wordevok(mtx,c)
#'@export

split_list_wordevok = function(list,groups)
{
  evo = list()
  for(i in 1:length(groups))
  {
    step1 = split(list,list$Source)
    step2 = step1[groups[[i]]]
    step3 = data.frame(rbindlist(step2))
    step4 = split(step3,step3$Target)
    step5 = step4[groups[[i]]]
    evo[[i]] = data.frame(rbindlist(step5))
  }
  return(evo)
}

#' Armazena a frequência de evocação e a ordem média de evocação de cada evocação
#' única padronizada.
#'
#' ....
#'@param wordevok Um objeto da classe wordevok.
#'@param groups Uma lista contendo os membros da respectiva comunidade.
#'....
#'@return \code{telp_wordevok} retorna uma lista que contém em cada dimensão sublistas
#'com um data.frame contendo as evocações únicas, ordem média de evocação e frequência
#'de evocação por comunidade.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'@export

telp_wordevok = function(wordevok,groups)
{

  cp = wordevok_comm_subsets(wordevok,groups)


  for(i in 1:length(cp))
  {
    cp[[i]]$IDCLASS = rownames(cp[[i]])
  }


  lcp = list()
  for (i in 1:length(cp))
  {
    lcp[[i]] = as.wordevok(cp[[i]],"IDCLASS")
  }


  fcp = list()
  for (i in 1:length(lcp))
  {
    fcp[[i]] = table(as.matrix(lcp[[i]]$Dataset))
  }


  wcp = list()
  for (i in 1:length(fcp))
  {
    wcp[[i]] = rownames(fcp[[i]])
  }


  ome = list()
  for (i in 1:length(lcp))
  {
    ome[[i]] = rep(1:length(wcp[[i]]))
    for(j in 1:length(wcp[[i]]))
    {
      k =  dim(lcp[[i]]$Dataset)[2]
      eval(parse(text=paste0("ome[[",i,"]][",j,"] = round((",paste(paste0(1:k,
                                                                          "*","sum(na.omit(as.integer(lcp[[",i,
                                                                          "]]$Dataset[,",1:k,"]==wcp[[",i,"]][",j,
                                                                          "])))"),collapse = "+"),")/fcp[[",i,"]][",
                             j,"],2)")))
    }
  }



  resini = list()
  for(i in 1:length(ome))
  {
    resini[[i]] = cbind(wcp[[i]],fcp[[i]],ome[[i]])
    rownames(resini[[i]])=rep("",nrow(resini[[i]]))
  }


  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  cfmean = list()
  for(i in 1:length(resini))
  {
    cfmean[[i]] = rep(0,2)
    cfmean[[i]][1] = round(mean(as.numeric(resini[[i]][,2])),2)
    cfmean[[i]][2] = round(mean(as.numeric(resini[[i]][,3])),2)
  }

  cfmedian = list()
  for(i in 1:length(resini))
  {
    cfmedian[[i]] = rep(0,2)
    cfmedian[[i]][1] = round(median(as.numeric(resini[[i]][,2])),2)
    cfmedian[[i]][2] = round(median(as.numeric(resini[[i]][,3])),2)
  }


  cfmode = list()
  for(i in 1:length(resini))
  {
    cfmode[[i]] = rep(0,2)
    cfmode[[i]][1] = round(mode(as.numeric(resini[[i]][,2])),2)
    cfmode[[i]][2] = round(mode(as.numeric(resini[[i]][,3])),2)
  }

  output = list()
  for (i in 1:length(resini))
  {
    output[[i]] = list()
    output[[i]][[1]] = resini[[i]]
    output[[i]][[2]] = cfmean[[i]]
    output[[i]][[3]] = cfmedian[[i]]
    output[[i]][[4]] = cfmode[[i]]
    names(output[[i]]) = c("resini","cfmean","cfmedian","cfmode")
    class(output[[i]]) = "TELP-wordevok"
  }
  class(output) = "Multiple-TELP-wordevok"
  return(output)
}

#' Apresenta o gráfico de quadrantes de Vergés.
#'
#' Esta função apresenta o gráfico de quadrantes de Vergés.
#'
#'@param twfile Uma dimensão de um objeto da classe \code{Multiple-TELP-wordevok}, gerado pela função \code{telp_wordevok}
#'@param ln Lógico. \code{TRUE} implica na utilização da escala logarítmica nos valores encontrados.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'telp_wordevok_plot(t[[1]])
#'@export

telp_wordevok_plot = function(twfile,ln=FALSE)
{
  if(class(twfile) != "TELP-wordevok")
  {
    stop("Erro ainda n?o definido.")
  }
  if(ln)
  {
    plot(log(as.numeric(twfile$resini[,2])),
         log(as.numeric(twfile$resini[,3])),
         xlab="Logaritmo da Frequencia de Evocacao",
         ylab="Logaritmo da Ordem Media de Evocacao",
         col=rainbow(nrow(twfile$resini)),
         pch=19)

    abline(v=log(as.numeric(twfile$cfmean[1])),
           h=log(as.numeric(twfile$cfmean[2])), col = "gray20")
    abline(v=log(as.numeric(twfile$cfmedian[1])),
           h=log(as.numeric(twfile$cfmedian[2])), col = "red")
    abline(v=log(as.numeric(twfile$cfmode[1])),
           h=log(as.numeric(twfile$cfmode[2])), col = "purple")
  } else {

    plot((as.numeric(twfile$resini[,2])),
         (as.numeric(twfile$resini[,3])),
         xlab="Frequencia de Evocacao",
         ylab="Ordem Media de Evocacao",
         col=rainbow(nrow(twfile$resini)),
         pch=19)

    abline(v=(as.numeric(twfile$cfmean[1])),
           h=(as.numeric(twfile$cfmean[2])), col = "gray20")
    abline(v=(as.numeric(twfile$cfmedian[1])),
           h=(as.numeric(twfile$cfmedian[2])), col = "red")
    abline(v=(as.numeric(twfile$cfmode[1])),
           h=(as.numeric(twfile$cfmode[2])), col = "purple")
  }
}

#' Divide as evocações por seus respectivos quadrantes de Vergès segundo o método selecionado.
#'
#'
#'@param mtwfile Um objeto da classe \code{Multiple-TELP-wordevok}, gerado pela função \code{telp_wordevok}.
#'@param method Criterio de discriminação dos quadrantes.
#'....
#'@return \code{telp_wordevok_quad} retorna uma lista contendo os quadrantes de Vergès
#'da representação social. Cada dimensão do objeto resultante representa uma comunidade, e cada subdimensão
#'representa um quadrante daquela comunidade.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'twq = telp_wordevok_quad(t,"mean")
#'@export

telp_wordevok_quad = function(mtwfile,method = "mean")
{
  quad = list()
  if (method == "mean")
  {
    m = 2
  } else {
    if (method == "median")
    {
      m = 3
    } else {
      m = 4
    }
  }
  for(i in 1:length(mtwfile))
  {
    quad[[i]] = list()
    resini = mtwfile[[i]]$resini
    resini = data.frame(resini)
    quad[[i]][[1]] = data.frame(resini[(as.integer(resini[,2])>=mtwfile[[i]][[m]][1]&as.numeric(resini[,3])<mtwfile[[i]][[m]][2]),])
    colnames(quad[[i]][[1]]) = c("Evocation","Frequency","AEO")
    quad[[i]][[2]] = data.frame(resini[(as.integer(resini[,2])>=mtwfile[[i]][[m]][1]&as.numeric(resini[,3])>=mtwfile[[i]][[m]][2]),])
    colnames(quad[[i]][[2]]) = c("Evocation","Frequency","AEO")
    quad[[i]][[3]] = data.frame(resini[(as.integer(resini[,2])<mtwfile[[i]][[m]][1]&as.numeric(resini[,3])<mtwfile[[i]][[m]][2]),])
    colnames(quad[[i]][[3]]) = c("Evocation","Frequency","AEO")
    quad[[i]][[4]] = data.frame(resini[(as.integer(resini[,2])<mtwfile[[i]][[m]][1]&as.numeric(resini[,3])>=mtwfile[[i]][[m]][2]),])
    colnames(quad[[i]][[4]]) = c("Evocation","Frequency","AEO")
    names(quad[[i]]) = c("Quad1","Quad2","Quad3","Quad4")
  }
  class(quad) = "TELP-wordevok-quad"
  return(quad)
}

#' Calcula as coordenadas da relevância para os signficados dos quadrantes de
#' uma determinada comunidade.
#'
#'@param quad Uma subdimensão de um objeto da classe \code{TELP-wordevok-quad} gerado pela função
#'telp-wordevok-quad.
#'@param classe Um documento listando as classes às quais pertencem as evocações
#'únicas do wordevok. Veja os detalhes para mais informações.
#'....
#'@return \code{wordevok_quad_class} Retorna as coordenadas de relevância relativa
#'de um determinado significado para um quadrante de uma comunidade.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'twq = telp_wordevok_quad(t,"mean")
#'
#'# Macro-groups
#'
#'u = noloop.pref$Evocations
#'mg = c("Optional","Optional","Fifth","Fifth","Second","First","Optional",
#'       "Eigth","Fifth","Third","Second","Forth","Optional","Forth","Sixth",
#'       "Third","Optional","Forth","Optional","Sixth")
#'u = data.frame(cbind(u,mg))
#'wqc= wordevok_quad_class(twq[[1]][[1]],u)
#'@export

wordevok_quad_class = function(quad,classe)
{
  classe = apply(classe,2,as.character)
  if(dim(quad)[1] == 1)
  {
    quad = apply(quad,2,as.character)
    quad = t(quad)
  }else{
    quad = apply(quad,2,as.character)
  }
  zata = unique(as.character(as.character(classe[,2])))
  radgeral = as.data.frame(matrix(ncol = length(zata),nrow=1))
  colnames(radgeral) = zata
  radgeral[1,] = rep(0,length(zata))
  radgeral[2:3,] = rep(0,length(zata))
  for(i in 1:dim(classe)[1]){
    search = match(as.character(classe[i,1]),quad[,1])
    if(!is.na(search))
    {
      radcol = match(as.character(classe[i,2]),colnames(radgeral))
      radgeral[3,radcol] = radgeral[3,radcol] + log(as.numeric(quad[search,2]))*(1/as.numeric(quad[search,3]))
    }
  }
  radgeral[1,] = max(radgeral[3,])
  return(radgeral)
}

#' Calcula as coordenadas da relevância para os signficados dos quadrantes de
#' multicomunidades.
#'
#'@param twqfile Um objeto da classe \code{TELP-wordevok-quad} gerado pela função
#'\code{telp-wordevok-quad}.
#'@param classe Um documento listando as classes às quais pertencem as evocações
#'únicas do wordevok. Veja os detalhes para mais informações.
#'....
#'@return \code{wordevok_radar_attr} retorna uma lista de coordenadas de relevância relativa
#'de um determinado significado para cada quadrante de cada comunidade inserida.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'twq = telp_wordevok_quad(t,"mean")
#'
#'# Macro-groups
#'
#'u = noloop.pref$Evocations
#'mg = c("Optional","Optional","Fifth","Fifth","Second","First","Optional",
#'       "Eigth","Fifth","Third","Second","Forth","Optional","Forth","Sixth",
#'       "Third","Optional","Forth","Optional","Sixth")
#'u = data.frame(cbind(u,mg))
#'wra= wordevok_radar_attr(twq,u)
#'@export

wordevok_radar_attr = function(twqfile, classe)
{
  radar = list()
  for(i in 1:length(twqfile))
  {
    c = 1
    n = NULL
    radar[[i]] = list()
    for(j in 1:length(twqfile[[i]]))
    {
      if(dim(twqfile[[i]][[j]])[1] != 0)
      {
        radar[[i]][[c]] = wordevok_quad_class(twqfile[[i]][[j]],classe)
        n[c] = names(twqfile[[i]])[j]
        c = c + 1
      }
    }
    names(radar[[i]]) = n
  }
  class(radar) = "TELP-wordevok-radar"
  return(radar)
}

#' Adequa as coordenadas para a utilização dos radares criados usando o \code{ggplot2}.
#'
#'@param radar_attr Um objeto da classe \code{TELP-wordevok-radar} criado pela
#'função \core{wordevok_radar_attr}.
#'....
#'@return \code{wordevok_radar_gg} retorna uma lista onde cada dimensão representa
#'o núcleo central da comunidade em questão bem como as relevâncias relativas de
#'cada conceito na mesma.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#' add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'twq = telp_wordevok_quad(t,"mean")
#'
#'# Macro-groups
#'
#'u = noloop.pref$Evocations
#'mg = c("Optional","Optional","Fifth","Fifth","Second","First","Optional",
#'       "Eigth","Fifth","Third","Second","Forth","Optional","Forth","Sixth",
#'       "Third","Optional","Forth","Optional","Sixth")
#'u = data.frame(cbind(u,mg))
#'wra= wordevok_radar_attr(twq,u)
#'wrg= wordevok_radar_gg(wra)
#'@export

wordevok_radar_gg = function(radar_attr)
{
  core = list()
  for(i in 1:length(radar_attr))
  {
    base = t(radar_attr[[i]][[1]])
    new  = data.frame(Community = paste("Community", i),
                      Core = rownames(base),
                      Priority = base[,3]/max(base[,3]))

    setorder(new,"Core")
    new[nrow(new)+1, ] = new[1,]
    rownames(new) = 1:length(rownames(new))
    core[[i]] = new
  }
  names(core) = paste0("Community_",1:length(radar_attr))
  class(core) = "wordevok-radar-gg"
  return(core)
}

#' Gera os gráficos de radar para o pensamento coletivo das comunidades usando \code{ggplot2}.
#'
#' ....
#'@param radar_gg Um objeto da classe \code{wordovok-radar-gg} criado pela
#'função \core{wordevok_radar_gg}.
#'....
#'@return \code{wordevok_radar_plot} retorna uma lista onde cada dimensão representa
#'um gráfico de radar do núcleo central da comunidade em questão.
#' ...
#'@author Wesley Henrique Silva Pereira
#'@examples
#'data("preferences")
#'mtx = wordevok_affinity_adjacency(noloop.pref)
#'g = graph_from_adjacency_matrix(mtx,mode = "undirected", weighted = TRUE,
#'add.rownames = "names")
#'c = cluster_louvain(g, weights = E(g)$weight)
#'t = telp_wordevok(noloop.pref,c)
#'twq = telp_wordevok_quad(t,"mean")
#'
#'# Macro-groups
#'
#'u = noloop.pref$Evocations
#'mg = c("Optional","Optional","Fifth","Fifth","Second","First","Optional",
#'       "Eigth","Fifth","Third","Second","Forth","Optional","Forth","Sixth",
#'       "Third","Optional","Forth","Optional","Sixth")
#'u = data.frame(cbind(u,mg))
#'wra = wordevok_radar_attr(twq,u)
#'wrg = wordevok_radar_gg(wra)
#'wrp = wordevok_radar_plot(wrg)
#'wrp[[1]]
#'@export

wordevok_radar_plot = function(radar_gg)
{
  for(i in 1:length(radar_gg))
  {
    radar_gg[[i]] = ggplot(radar_gg[[i]], aes(x = Core, y = Priority)) +
      geom_polygon(aes(group = Community), color = "#006666", fill="#006666", size = 0.5, alpha=0.6) +
      ggtitle(as.character(radar_gg[[i]]$Community[1])) +
      xlab("") +
      ylab("") +
      ylim(-0.1,1) +
      guides(color = guide_legend(ncol=2)) +
      coord_radar() +
      guides(colour=guide_legend(nrow=4, byrow=TRUE), shape=guide_legend(nrow=1, byrow=TRUE)) +
      theme(
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.key = element_blank(), legend.title = element_blank(),
        legend.background = element_rect(color="#ffffff", fill="transparent"), ### neu check !!!
        panel.background = element_rect(fill = "white", colour = "white", size = 0.1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "#dddddd")
      )
  }
  return(radar_gg)
}


coord_radar <- function (theta = "x", start = 0, direction = 1)
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
  {
    "y"
  } else {
    "x"
  }
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, direction = sign(direction), is_linear = function(coord) TRUE)
}
