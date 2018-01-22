# Function for plotting user-specified color vectors on to cancer maps.
# Parameters:
# X: "Coordinate file" of a generated cancer map, where
#    1st column is the ID a sample,
#    2nd and 3rd column are the x and y coordinates respectively, and
#    4th column is the cluster label.
# NAME: Name to be used on the output file.
# SIZE: Variable defining point size in cancer maps.
# color: user-specified color vector. Same length as number or samples in X.
# TITLE: Title for the plot.
# PATH_OUTPUT: path where to write output files.
Plot_color_vector=function(X, NAME, SIZE, color, TITLE, PATH_OUTPUT, peaks) {
  
  # Color vector
  datCol=color
  
  # Prepare input matrix for plotting.
  dat2show <- cbind(X$x, X$y)
  df=as.data.frame(dat2show)
  colnames(df) = c("X1","X2")

  # If for some reason we use a color vector which includes blanks,
  # create a color vector without them to be plotted at the front.
  # (Increases visibility of the "important" colored samples.)
  front=!datCol==""

  # Call actual plotting function.
  p=drawFig(df, CLUSTER_CENTRE, datCol, front, TITLE, SIZE, peaks)

  # Write out print quality figure as PDF.
  ggsave(paste0(PATH_OUTPUT, NAME, "_singlepage.pdf"), p, width = 210, height = 210, units = "mm", dpi=150)
}
