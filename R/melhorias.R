w_a_a = function(wordevok)
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
  wordevok$Dataset = apply(wordevok$Dataset,2,as.character)
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

w_a_a = function(wordevok)
{
  adjacency = matrix(nrow = dim(wordevok$Dataset)[1], ncol = dim(wordevok$Dataset)[1])
  colnames(adjacency) = as.character(wordevok[[1]])
  rownames(adjacency) = colnames(adjacency)
  adjacency[,] = 0
  z = apply(wordevok$Dataset,2,as.character)
  rownames(z) = rownames(wordevok$Dataset)
  
  while(nrow(z)>1)
  {
    v = z[1,,drop=FALSE]
    name = rownames(v)
    w = z[-1,,drop = FALSE]
    c = t(apply(w,1,match,v[!is.na(v)]))
    rownames(c) = rownames(w)
    s = apply(!is.na(c),1,sum)
    if(sum(s) == 0)
    {
      z = z[-1,,drop = FALSE]
      next
    } else {
      k = z[names(s[s>0]),,drop=FALSE]
      adjacency[name,names(s[s>0])] = apply(k,1,affinity,v)
      z = z[-1,,drop = FALSE]
    }
  }
  
  adjacency[lower.tri(adjacency)] = t(adjacency)[lower.tri(adjacency)]
  return(adjacency)
}

w_a_a = function(wordevok)
{
  adjacency = matrix(nrow = dim(wordevok$Dataset)[1], ncol = dim(wordevok$Dataset)[1])
  colnames(adjacency) = as.character(wordevok[[1]])
  rownames(adjacency) = colnames(adjacency)
  adjacency[,] = 0
  z = apply(wordevok$Dataset,2,as.character)
  rownames(z) = rownames(wordevok$Dataset)
  
  while(nrow(z)>1)
  {
    v = z[1,,drop=FALSE]
    name = rownames(v)
    w = z[-1,,drop = FALSE]
    c = t(apply(w,1,match,v[!is.na(v)]))
    rownames(c) = rownames(w)
    s = apply(!is.na(c),1,sum)
    if(sum(s) == 0)
    {
      z = z[-1,,drop = FALSE]
      next
    } else {
      k = z[names(s[s>0]),,drop=FALSE]
      adjacency[name,names(s[s>0])] = apply(k,1,affinity,v)
      z = z[-1,,drop = FALSE]
    }
  }
  
  adjacency[lower.tri(adjacency)] = t(adjacency)[lower.tri(adjacency)]
  return(adjacency)
}