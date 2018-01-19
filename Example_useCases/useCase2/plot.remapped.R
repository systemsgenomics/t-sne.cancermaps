library(ggplot2)

plot.remapped=function(coord_A, coord_B, COLOR_A, COLOR_B, PATH_OUTPUT=getwd(), NAME="Remapped_samples", ALPHA_A=0.6, SIZE_A=3, ALPHA_B=0.8, SIZE_B=2){
  
  p = ggplot(coord_A, aes(x, y)) + geom_point(alpha=ALPHA_A, size=SIZE_A,
                                                         colour = COLOR_A)
  
  p = p + geom_point(data = data.frame(coord_B) ,alpha = ALPHA_B, aes(X1, X2),
                     colour = COLOR_B, size=SIZE_B)
  p = p + theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())
  ggsave(paste0(PATH_OUTPUT, NAME, "_singlepage.pdf"), p, width = 210, height = 210, units = "mm", dpi=150)

}
