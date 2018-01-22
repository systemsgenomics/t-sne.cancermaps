# Function for finding correlated genes in a cluster
# Parameters:
# clust: Cluster label from which cluster to find correlated genes.
# data: data matrix where rows contain genes and columns contain samples.
# anno: cancer map "coordinate file" where 4th row contains the clustering info.
Find_correlated_genes=function(clust, data, anno) {
  
  # Generate discrete vector for for clust.
  V1=as.numeric(anno[,4]==clust)

  # Calculate correlation and p-values of genes vs clust vector.
  genes=apply(data, 1, cor, V1, method="spearman")
  genes_pval=suppressWarnings(apply(data, 1, function(V2)cor.test(V2, V1, method="spearman")$p.value))
  
  # Generate result matrix.
  df=cbind(rownames(data), genes, genes_pval)
  df=df[order(as.numeric(df[,2]), as.numeric(df[,3]), decreasing = T),] # Order by correlation and p-value.
  df=df[!genes_pval>PVAL,] # Remove genes if p-value over threshold.
  
  # Report if no significant correlations found.
  if(dim(df)[2]<1){
    print(paste("no significantly correlated genes for:", clust))
    return(NULL)
  }
  
  # Get NUM_GENES most positively and most negatively correlated genes.
  top_genes=head(rownames(df), NUM_GENES)
  bottom_genes=tail(rownames(df), NUM_GENES)
   
  # Modify output and return genelists.
  out=list(top_genes, bottom_genes)
  names(out)=c(paste0("CLUSTER_", clust, "_UP"), paste0("CLUSTER_", clust, "_DOWN"))
  return(out)
}
