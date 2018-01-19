# Function to get AML annotation from the master annotation table.
#
# This function simplifies the GSE information to one vector and
# greps the essential phenotype information to one vector.

getAMLannotation <- function(anno_aml) {
  
  # Simplify GSEs
  U = as.character(unique(anno_aml$altGSE))
  U = U[!(U%in%"na")]
  for (i in length(U):1) {
    g = unlist(strsplit(U[i],":"))
    U = c(U,g)
    U = U[-i]
  }
  U = unique(c(U,as.character(unique(anno_aml$GSE.identifier..experiment.))))
  gses = U
  for (i in grep(",",gses,value=T)) {
    s = unlist(strsplit(i,","))
    g = paste(s[1],s[2],sep="|")
    GSE = grepl(g,anno_aml$GSE.identifier..experiment.)
    altGSE = grepl(g,anno_aml$altGSE)
    if (sum(GSE) > 0) { anno_aml$GSE.identifier..experiment.[GSE] = i }
    if (sum(altGSE) > 0) { anno_aml$altGSE[altGSE] = i }
  }
  
  # Simplify cytogenetics annotation
  cytogenetics <- rep("other", nrow(anno_aml))
  t8_21 <- cytogenetics[grepl("8;21", anno_aml$Cytogenetics)] <- "t8;21"
  MLL_aml <- cytogenetics[grepl("MLL", anno_aml$Cytogenetics)] <- "MLL"
  t15_17 <- cytogenetics[grepl("15;17", anno_aml$Cytogenetics)] <- "t15;17"
  ab16 <- cytogenetics[grepl("idt", anno_aml$Cytogenetics) | grepl("inv\\(16\\)", anno_aml$Cytogenetics) | grepl("inv16", anno_aml$Cytogenetics)] <- "ab16"
  anno_aml$Cytogenetics <- cytogenetics
  
  anno_aml
}
