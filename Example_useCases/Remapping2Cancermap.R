# *******************************************************************************************
#       FUNCTION: Remap new samples to TSNE map coordinate space | 15/01/2017 | Sergei Häyrynen, Petri Pölönen and Juha Mehtonen
# *******************************************************************************************

#********************************************************************************************
#************* RUNNING NOTES: run this script from the folder with the Rdata and
#*************                and R scripts. No need to modify this script!
#********************************************************************************************

#******************************** Load functions and libraries **************************************************

source("useCase2/run_remapping.R")
source("useCase2/plot.remapped.R")

#******************************** Load input data *****************************************************

# load hemap data, symbols as rows, gsm as column names:
data=get(load("DATA/useCase1/HEMAP_data.Rdata"))

# set name for output
NAME="HEMAP_ALL"

# Create new output directory for results, if it doesn't already exist.
# By default creates the output directory under current working directory.
PATH_OUTPUT="output_useCase2/"
dir.create(file.path(getwd(), PATH_OUTPUT), showWarnings = FALSE, recursive=T)

# TSNE coords to use:
coord = read.delim('DATA/useCase2/pre-B-ALL_T-ALL_15pct_genes_BHSNE_mean-shift.txt', stringsAsFactors = F)

# genes used in TSNE projection, can be used to filter:
top15var_genes = scan('DATA/useCase2/data9544_15pct_topVariating_genelist.txt',"a", quiet = T)

# subset to genes used in original projection and take only samples in TSNE projection:
data=data[rownames(data)%in%top15var_genes,colnames(data)%in%coord[,1]]

# Subset coordinate columns from the existing TSNE projection
coord_orgMap <- coord[,2:3]

# read newData, symbols as rows, gsm as column names
newData=get(load("DATA/useCase2/GSE49032_108_expressions.Rdata"))

# colors for remapped samples cytogenetics
col2show_remapped = scan('DATA/useCase2/GSE49032_108_cytogenetics_colors.txt',"a", quiet = T)

#******************************** run remapping algorithm, vizualize result *****************************************************

# input should be: original map data, map coordinates, VAR, new data to add
print(paste0("remapping samples to cancermap: ", NAME, "..."))
remapped_coord = run_remapping(originalData = data, map_coord = coord_orgMap, newData = newData, CORES=4)
print(paste0("remapping samples to cancermap: ", NAME, "..."))

# visualize remapped samples and cytogenetics on the original cancermap
col2show=rep("lightgrey", nrow(coord))

print(paste0("plotting remapped samples and cytogenetics to cancermap: ", NAME, "..."))
plot.remapped(coord_A = coord_orgMap, coord_B = remapped_coord, COLOR_A=col2show, COLOR_B=col2show_remapped, PATH_OUTPUT = PATH_OUTPUT, NAME = NAME)
print(paste0("plotting remapped samples and cytogenetics to cancermap: ", NAME, "... Done"))

