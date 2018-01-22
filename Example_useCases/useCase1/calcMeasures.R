library(Rtsne) # BH-SNE implementation.
library(LPCM) # mean-shift clustering implementation.
library(foreach) # loop structure

# Function to run BH-SNE + mean-shift and calculate NMI
calcMeasures <- function(data, phenotype, sample_series, PCA = F) {
  
  if(PCA) { 
    dd = seq(20, 100, 5)
  } else {
    dd = seq(2.5, 50, 2.5)
    g = apply(data, 2, var)
    g = sort(g, decreasing = T)
    g = names(g)
  }
  
  RESULT <- foreach(i = dd, .combine = "rbind", .inorder = F) %dopar% {
    
    set.seed(1)
    
    # Run BH-SNE and mean-shift clustering.
    if(PCA) x = Rtsne(data, check_duplicates = F, pca = T, is_distance = F, initial_dims = i)
    if(!PCA) {
      n = round(ncol(data) * i / 100)
      G = g[1:n]
      D = data[,colnames(data)%in%G]
      x = Rtsne(D, check_duplicates = F, pca = F)
    }
    x = x$Y
    m1 = ms(x, h=h, scaled=F, thr=0.01, iter=500, plotms=0)
    X = data.frame(rownames(data), x, m1$cluster.label)
    
    # Modify clustering based on phenotype of samples inside cluster by majority vote.
    C = X[,4]
    ni = unique(C)
    for(f in ni) {
      C[C==f] = names(tail(sort(table(phenotype[C==f])),1))
    }
    
    # Calculate NMI between phenotypes and clusters.
    nMI = NMI(phenotype, C) # colorClass contains phenotype.
    
    # Calculate NMI between data series and clusters.
    nMI_lab = NMI(sample_series, X[,4]) # GSE is the dataseries ID.
    
    # Combine data to result matrix R.
    tmp = data.frame(cbind(i,length(ni),nMI,nMI_lab))
    colnames(tmp) = c("i", "no_clusters", "pNMI", "eNMI")
    
    tmp
  }
  
  RESULT
  
}
