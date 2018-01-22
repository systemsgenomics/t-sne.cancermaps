# **********************************************************************************************************
#       FUNCTION: Evaluate NMI on cancer-map | 15/01/2018 | Juha Mehtonen and Petri Pölönen
# **********************************************************************************************************

#***********************************************************************************************************
#************* RUNNING NOTES: run this script from the folder with the Rdata and ***************************
#*************                and R scripts. No need to modify this script!      ***************************
#***********************************************************************************************************

# Hemap: Testing how the percentage of most variable genes retained or the number of principal components retained affect the cancer-map.
# NMI is used to quantify both phenotype (pNMI) and sample series (eNMI) based separation.

#************************************* Load input data *****************************************************

# Load gene expression data and annotations.
load("DATA/useCase1/HEMAP_data.Rdata") # 'data' variable.
anno = read.delim("DATA/useCase1/USE_final_anno_columnIDsorted_data9544_withColorClass.txt", stringsAsFactors = F)

# Subset AML
data <- t(matrix[,anno$colorClass=="AML"])
anno <- anno[anno$colorClass=="AML",]

# Use helper function to simplify GSE and phenotype information in the annotation
source("useCase1/getAMLannotation.R")
aml_anno <- getAMLannotation(anno)

# Create new output directory for results, if it doesn't already exist.
# By default creates the output directory under current working directory.
PATH_OUTPUT="output_useCase1/"
dir.create(file.path(getwd(), PATH_OUTPUT), showWarnings = FALSE, recursive=T)

#***************** Step 1: Calculate NMI for clusters based on phenotype and sample series *****************

# Load helper functions
source("useCase1/infoMeasures.R") # NMI calculation functions
source("useCase1/calcMeasures.R") # Function to automate NMI calculation

# Parallelize
library(doParallel)
nCores <- min(detectCores(), 20)
registerDoParallel(cores = nCores)

# Bandwidth parameter for mean-shift clustering.
h = 1.5

# Run calcMeasures with gene selection and with PCA
res.gs <- calcMeasures(data = data, phenotype = aml_anno$Cytogenetics, sample_series = aml_anno$GSE.identifier..experiment., PCA = F)
res.pca <- calcMeasures(data = data, phenotype = aml_anno$Cytogenetics, sample_series = aml_anno$GSE.identifier..experiment., PCA = T)

#******************************* Step 2: Plot NMI measures ******************************

# Modify data.frames for plotting
library(reshape2)
res.gs <- melt(data = res.gs, id.vars = 1:2, measure.vars = 3:4)
res.pca <- melt(data = res.pca, id.vars = 1:2, measure.vars = 3:4)

# Load library for plots
library(ggplot2)

# Generate plots
p.gs <- ggplot(res.gs, aes(x = i, y = value, color = variable)) + geom_line() + scale_colour_manual(values=c("blue", "red"))
p.gs <- p.gs + labs(title = "Gene selection", x = "% genes", y = "NMI")
p.pca <- ggplot(res.pca, aes(x = i, y = value, color = variable)) + geom_line() + scale_colour_manual(values=c("blue", "red"))
p.pca <- p.pca + labs(title = "PCA", x = "PCs", y = "NMI")

# Print to PDF
pdf(file.path(PATH_OUTPUT, "NMI_results.pdf"))
plot(p.gs)
plot(p.pca)
graphics.off()
