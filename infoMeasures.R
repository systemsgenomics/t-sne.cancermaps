# Information theory measures for cluster evaluation

entropy = function(L, normalize = FALSE) {
  n = length(L)
  H = 0
  k = unique(L)
  for(i in k) H = H - (sum(L==i)/n) * log(sum(L==i)/n)
  if(normalize) H = H / log(n)
  if(is.na(H)) H = 0
  H
}

mutualinformation = function(L1,L2) {
  if(length(L1) != length(L2)) stop("Label lengths don't match.")
  n = length(L1)
  k1 = unique(L1)
  k2 = unique(L2)
  I = 0
  for(i in k1) {
    pi = sum(L1 == i) / n
    for(j in k2) {
      pij = sum(L1 == i & L2 == j) / n
      pj = sum(L2 == j) / n
      if (pij > 0) I = I + pij * log(pij/(pi*pj))
    }
  }
  I
}

NMI = function(L1,L2) {
  if(length(L1) != length(L2)) stop("Label lengths don't match.")
  I = mutualinformation(L1,L2)
  H1 = entropy(L1)
  H2 = entropy(L2)
  NMI = I / ((H1+H2)/2)
  NMI
}

variationinfo = function(L1,L2) {
  if(length(L1) != length(L2)) stop("Label lengths don't match.")
  I = mutualinformation(L1,L2)
  H1 = entropy(L1)
  H2 = entropy(L2)
  VI = H1 + H2 - 2*I
  VI
}
