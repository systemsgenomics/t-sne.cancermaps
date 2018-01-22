# Function for creating the cancer map.
# Input parameters:
# data: data matrix where rows represent samples and columns represent features.
# name: Parameter used in the naming of the output files.
# VAR: percentage of most variable features to retain in the data matrix for t-SNE.
# BW: Bandwidth parameter for mean-shfit clustering. Lower the bandwidth, more dense the clusters found.
# PATH_OUTPUT: Path where to save the output files.
CancerMap=function(data, name, VAR=NULL, BW=0.9, PATH_OUTPUT){
  
  if(!is.null(VAR)){
    # Get the j% most variating genes, feature selection.
    v = apply(data,2,var)
    v = sort(v,decreasing=T)
    v = names(v)
    j = VAR
    V = round(dim(data)[2] * j / 100)
    v = v[1:V]
    data = data[,colnames(data)%in%v]
    
    outname=paste0(PATH_OUTPUT, "cancermap_", name, "_", VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW, ".txt")
    outnamecent=paste0(PATH_OUTPUT, "cancermap_", name, "_", VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW, "_cluster_centroids.txt")
    
  }else{
    outname=paste0(PATH_OUTPUT, "cancermap_", name, "_", "BH-SNE_mean-shift_BW", BW, ".txt")
    outnamecent=paste0(PATH_OUTPUT, "cancermap_", name, "_", "BH-SNE_mean-shift_BW", BW, "_cluster_centroids.txt")
  }

  # Run BH-SNE + mean-shift clustering.
  set.seed(1) # Random number generator seed for reproducible results.
  x = Rtsne(data,perplexity = 30, check_duplicates = F, pca = F, is_distance = F) # BH-SNE.
  x = x$Y
  h = BW # Bandwidth parameter for clustering (for subsets of Hemap data we used h = 1.5).
  m1 = ms(x, h=h, scaled=F, thr=0.01, iter=500, plotms=0) # Mean-shift clustering.
  X = data.frame(rownames(data), x, m1$cluster.label) # Summarize results to a data.frame.
  colnames(X) = c("ID","x","y",paste0("cluster"))
  
  # Write embedding coordinates and the clustering result to a file.
  write.table(X, outname, quote = F, sep = "\t", row.names = F)

  # Write the cluster centroids(modes) to a file.
  write.table(data.frame(m1$cluster.center), outnamecent, row.names = F, sep = "\t", quote = F)
  return(X)
}
