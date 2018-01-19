# Function for plotting a generated cancer map and clustering.
# Parameters:
# X: "Coordinate file" of a generated cancer map, where
#    1st column is the ID a sample,
#    2nd and 3rd column are the x and y coordinates respectively, and
#    4th column is the cluster label.
# peaks: contains the cluster centroid coordinates.
# CLUSTER_CENTRE: Boolean variable whether to plot cluster centroids to cancer maps.
# NAME: Name to be used on the output file.
# SIZE: Variable defining point size in cancer maps.
# VAR: Percentage of most variable genes included when generating X. Used in naming the output file.
# BW: Bandwidth parameter for mean-shift clustering. Used in naming the output file.
# TITLE: Title for the plot.
# PATH_OUTPUT: path where to write output files.
Plot_cancermap_clusters=function(X, peaks, CLUSTER_CENTRE, NAME, SIZE, VAR, BW, TITLE, PATH_OUTPUT) {

  # Take cluster info from the 4th column in X.
  classes=as.character(X[,4])

  # Create color vector to be used in plotting.
  size=c(length(unique(classes)))
  COL<-brewer.pal(9,"Set1")
  col1 <- colorRampPalette(COL)(size)
  datCol=unlist(lapply(classes, function(c){ind=unique(classes)%in%c
                                            col1[ind]}))
  
  # Prepare input matrix for plotting.
  dat2show <- cbind(X$x, X$y)
  df=as.data.frame(dat2show)
  colnames(df) = c("X1","X2")
  
  # If for some reason we use a color vector which includes blanks,
  # create a color vector without them to be plotted at the front.
  # (Increases visibility of the "important" colored samples.)
  front=!datCol==""
  
  # Generate name
  MAPNAME=paste0("cancermap_", NAME, "_", VAR, "pct_genes_", "BH-SNE_mean-shift_BW", BW)
  
  # Call actual plotting function.
  p=drawFig(df, CLUSTER_CENTRE, datCol, front, TITLE, SIZE, peaks)

  # Write out print quality figure as PDF.
  ggsave(paste0(PATH_OUTPUT, MAPNAME, "_singlepage.pdf"), p, width = 210, height = 210, units = "mm", dpi=150) 
}
