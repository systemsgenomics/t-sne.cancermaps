# Load needed packages.
library(reshape2)
library(gridExtra)
library(GSVA)
library(RColorBrewer)
library(ggplot2)
library(Rtsne)
library(LPCM)

# Load script for generating cancer maps.
source("useCase3/Cancermap.R")

# Load script for finding correlated genes and making cluster specific genesets.
source("useCase3/Find_correlated_genes.R")

# Load various scripts used in plotting the cancer maps and results.
source("useCase3/drawFig.R")
source("useCase3/Plot_color_vector.R")
source("useCase3/Plot_cancermap_clusters.R")
source("useCase3/Plot_GSVA_scores.R")

# Parameters used in generating genesets per cluster.
PVAL=0.05 # P-value threshold
NUM_GENES=20 # Number of genes to select per geneset.

# Parameters used in plotting.
VALUE=0.4 # GSVA score cutoff (default is 0.4 when NUM_GENES is 20).
CLUSTER_CENTRE=T # Boolean variable whether to plot cluster centroids to cancer maps.
SIZE=0.7 # Variable defining point size in cancer maps.

# Parameters for cancer map generation.
VAR=15 # Percentage of most variable genes to be included in the dataset when generating cancer map.
BW=0.9 # Bandwidth parameter for mean-shift clustering for new data, which typically has less samples and smaller bandwith
BW_HEMAP=1.5 # Bandwidth parameter for mean-shift clustering for hemap. Hemap bandwith was 2.5 and Hemap AML/ALL bandwith was 1.5

# Names to be used on the output file.
NAME="Ross_rma_u133a"
PATHW="clusters_top20"

# set hemap name for output
HEMAP="HEMAP_ALL"

# Create new output directory for results, if it doesn't already exist.
# By default creates the output directory under current working directory.
PATH_OUTPUT="output_useCase3/"
dir.create(file.path(getwd(), PATH_OUTPUT), showWarnings = FALSE, recursive=T)

print("Default parameters and scripts loaded.")

