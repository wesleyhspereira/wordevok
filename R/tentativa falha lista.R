wordevok = data

word = rownames(wordevok$Dataset)
c = t(combn(word,2))
z = cbind(s = match(c[,1],rownames(data$Dataset)),
      t = match(c[,2],rownames(data$Dataset)))

 

m = apply(t,1,k)

m = t[1,]

k = function(m)
{
  m = as.character(m)
  g = as.data.frame(
      as.matrix(
      as.data.frame(
      t(
      data.frame(x = t(data$Dataset[m[1],]), 
                 y = t(data$Dataset[m[2],])
                 )
                 )
                 )
                 )
                 )
  
  if(sum(as.matrix(g)[1,!is.na(as.matrix(g)[1,])] %in% 
         as.matrix(g)[2,!is.na(as.matrix(g)[2,])]) == 0)
  {
    return(0)
  } else {
    return(affinity(g))
  }
}

k(m)



affinity = function(comp)
{
  subject1 = comp[1,]
  subject2 = comp[2,]
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

f = function(p)
{
  n = wordevok$Dataset
  n$id = rownames(n)
  s = eval(parse(text = paste("filter(n, ", 
              paste0(colnames(wordevok$Dataset), "== '", 
              p ,"'",collapse = "|"),")")))
  if(nrow(s) > 1)
  {
    k = s$id
    if(nrow(s) == 2)
    {
      m = matrix(0,ncol = 2, nrow = 2)
      m[1,] = combn(k,2)
      m[2,] = combn(k,2)
      k = m
    } else{
      k = t(combn(k,2))
    }
    k = apply(k,2,as.numeric)
    k = as.data.frame(k,drop = F)
    return(k) 
  }
}
p = "Medo"

f(p)

p = wordevok$Evocations
t = sapply(p,f)

names(t) <- seq_along(t)
t = Filter(Negate(is.null), t)


v = sapply(g,t)
u = sapply(v,as.data.frame)


t = as.data.frame(rbindlist(t))
t = t[!duplicated(t),]

length(t)

combn(t[[1]],2)

c = apply(t,1,k)
dim(t)

comp = t
comp[1]
