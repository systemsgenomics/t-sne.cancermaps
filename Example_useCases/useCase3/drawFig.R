# function for plotting
drawFig=function(df, addPeaks, col2show, front,title, SIZE, peaks){
  
  df=df[order(front),]
  col2show=col2show[order(front)]
  
  p <- ggplot(df, aes(X1, X2)) + geom_point(alpha=0.6, size=SIZE, colour=col2show) + 
    theme(text=element_text(size=8), title=element_text(size=6), axis.title=element_text(size=8))
  
  if(addPeaks){
    p <- p + geom_point(data=peaks, aes(X1, X2), colour = "red", size=4,fill=NA, shape=1)
    for (j in 1:dim(peaks)[1]){
      p <- p + annotate("text", label = as.character(j), x = peaks[j,1], y = peaks[j,2], size=2)
    }
  }
  p<- p + ggtitle(title) + xlab("y1") + ylab("y2")
  return(p)
}

