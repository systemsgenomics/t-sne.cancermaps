# Function for plotting the GSVA results on top of the cancer map.
# Parameters:
# feat: Feature to plot. Row name from data_plot.
# data_plot: Matrix of GSVA scores where rows contain features (genesets) and columns contain samples.
# VALUE: GSVA score cutoff used for the GSVA run.
# SIZE: Variable defining point size in cancer maps.
# CLUSTER_CENTRE: Boolean variable whether to plot cluster centroids to cancer maps.
Plot_GSVA_scores=function(feat, data_plot, VALUE, SIZE, CLUSTER_CENTRE, coord, peaks){

  # Find specified feature(feat) from data_plot
  data=data_plot[rownames(data_plot)%in%feat,]
  
  # Transform matrix to numeric. 
  data=as.numeric(data)
  
  # Color vector for gradient colors from blue to red.
  rbPal <- colorRampPalette(c('blue','red'))
  
  # Adjust data for gradient colors.
  data=c(data, 1, -1) # adjust range
  datCol <- rbPal(10)[as.numeric(cut(data,breaks = 10))]
  datCol=datCol[-c(length(datCol)-1, length(datCol))]
  data=data[-c(length(data)-1, length(data))]
  
  # Samples below cutoff colored grey.
  datCol[abs(data)<VALUE]="grey75"
  front=abs(data)>VALUE
  
  # Prepare coordinate data for plotting.
  dat2show <- cbind(coord$x, coord$y)
  df=as.data.frame(dat2show)
  colnames(df) = c("X1","X2")
 
  # Generate plot title. 
  cutoff=paste("GSVA score >", VALUE)
  plotname=paste(feat, cutoff, sep="\n")
  
  # Call actual plotting function.
  drawFig(df, CLUSTER_CENTRE, datCol, front, plotname, SIZE, peaks) 
}
