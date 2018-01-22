# *******************************************************************************************
#       FUNCTION: detect clusters on cancermap | 15/01/2018 | Petri Pölönen and Juha Mehtonen
# *******************************************************************************************

#********************************************************************************************
#************* RUNNING NOTES: run this script from the folder with the Rdata and
#*************                and R scripts. No need to modify this script!
#********************************************************************************************



#******************************** Load default parameters **************************************************

source("useCase3/load_scripts_set_parameters.R") # Loads default parameters and scripts

#******************************** Load data input data *****************************************************

source("useCase3/load_data_useCase3.R") # modify accordingly or use ROSS and Hemap ALL example data
#objects loaded: newdata, newdata_color_vector, matrix, hemap_color_vector

#********************************************* Step 1: Make cancermaps for both datasets ******************************************

# new data cancermap
print(paste0("Running TSNE for: ", NAME, "..."))
clust=CancerMap(data = t(newdata), name = NAME, VAR = VAR, BW = BW, PATH_OUTPUT = PATH_OUTPUT)
print(paste0("Running TSNE for: ", NAME, "... Done"))

# hemap cancermap
print(paste0("Running TSNE for: ", HEMAP, "..."))
clust_hemap=CancerMap(data = t(matrix), name = HEMAP, VAR=VAR, BW = BW_HEMAP, PATH_OUTPUT = PATH_OUTPUT)
print(paste0("Running TSNE for: ", HEMAP, "... Done"))

#********************************************* Make cancermap plot, color clusters and cytogenetics *******************************

print(paste0("plotting cancermap for: ", NAME, "..."))

# read cancermap data
X=read.delim(paste0(PATH_OUTPUT, "cancermap_", NAME, "_", VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW, ".txt"), header=T, stringsAsFactors=F)
peaks=read.delim(paste0(PATH_OUTPUT, "cancermap_", NAME, "_", VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW, "_cluster_centroids.txt"), header=T, stringsAsFactors=F)

Plot_cancermap_clusters(X, peaks, CLUSTER_CENTRE, NAME, SIZE, VAR, BW, NAME, PATH_OUTPUT) # Plot clusters with different colors

print(paste0("plotting cancermap for: ", NAME, "... Done"))

print(paste0("plotting cytogenetics to cancermap: ", NAME, "..."))

Plot_color_vector(X, NAME, SIZE, newdata_color_vector, NAME, PATH_OUTPUT, peaks) # Plot cytogenetics

print(paste0("plotting cytogenetics to cancermap: ", NAME, "... Done"))


#************************************** Step 2: Create genesets to identify cancermap clusters *************************************

print("making genesets")

#clusters
clusters=unique(X[,4])

# make gene sets
genesets=unlist(lapply(clusters, Find_correlated_genes, newdata, X), recursive=F)

print("genesets ready")

#****************************************************** run GSVA ************************************************************

# run GSVA
if(file.exists(paste0(PATH_OUTPUT, HEMAP, "_", NAME, "_", PATHW, "_GSVA.Rdata"))){
  load(paste0(PATH_OUTPUT, HEMAP, "_", NAME, "_", PATHW, "_GSVA.Rdata"))
  print("previously computed GSVA matrix loaded")
  print("if you want to re-compute, remove/rename file:")
  print(paste0(PATH_OUTPUT, HEMAP, "_", NAME, "_", PATHW, "_GSVA.Rdata"))
}else{
  print("GSVA running...")
  gsva_es <- gsva(as.matrix(matrix), method="gsva", genesets, mx.diff=F, tau=0.25, verbose=T, rnaseq=F, min.sz=5, max.sz=500)
  gsva_es=gsva_es$es.obs
  save(gsva_es, file=paste0(PATH_OUTPUT, HEMAP, "_", NAME, "_", PATHW, "_GSVA.Rdata"))
  print("GSVA done")
}

#********************************************* Step 3: Plot GSVA scores to hemap cancermap ************************************************************

X_hemap=read.delim(paste0(PATH_OUTPUT, "cancermap_", HEMAP, "_",VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW_HEMAP, ".txt"), header=T, stringsAsFactors=F)
peaks_hemap=read.delim(paste0(PATH_OUTPUT, "cancermap_", HEMAP, "_",VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW_HEMAP, "_cluster_centroids.txt"), header=T, stringsAsFactors=F)

feats=rownames(gsva_es)

print("making plots...")

# Plot GSVA scores to Hemap cancermap, obtain plotting data
p.all=lapply(feats, Plot_GSVA_scores, gsva_es, VALUE, SIZE, CLUSTER_CENTRE, X_hemap, peaks_hemap)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0(PATH_OUTPUT, HEMAP, "_", NAME, "_", PATHW, "_multipage.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)
print(paste0("done, file: ", paste0(PATH_OUTPUT, NAME, "_", PATHW, "_multipage.pdf")))

#********************************************* Plot cytogenetic data to hemap cancermap ************************************************************

print(paste0("plotting cytogenetics to cancermap: ", HEMAP, "..."))

Plot_color_vector(X_hemap, HEMAP, SIZE, hemap_color_vector, HEMAP, PATH_OUTPUT, peaks_hemap) # Plot cytogenetics

print(paste0("plotting cytogenetics to cancermap: ", HEMAP, "... Done"))
