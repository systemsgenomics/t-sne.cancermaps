#********************************************************************************************
#************* RUNNING NOTES: See input data formats from user guide
#********************************************************************************************



# Load new dataset
#********************************************************************************************

#newdata=read.table("exprs.txt",  sep="\t",  header=T,  row.names=1,  stringsAsFactors=F) # also table is fine
load("DATA/useCase3/ROSS_ALL.Rdata") # Rdata for faster loading, object newdata

# Load color vectors of cytogenetics, mutations or clinical annotations to be used for plotting with New data cancer maps
# When plotting different datasets, user has to naturally provide different color vectors to match the samples in
# corresponding gene expression matrix.
newdata_color_vector=scan("DATA/useCase3/ROSS_cytogenetics_colors.txt", "cytogenetics", quiet=TRUE)

#********************************************************************************************


# load Hemap data                            
#********************************************************************************************


# Load Hemap input data matrix and subtype information.                                     
load("DATA/useCase1/HEMAP_data.Rdata")

# By default, Hemap ALL is used, and we choose only the ALL samples from the Hemap data matrix using a vector defining subtypes for samples.
# If user wants to use a different subtype, the 'subtypes' vector contains all the different subtypes in our data and can be
# used to select another subset. f.e. selecting Multiple Myeloma: matrix[,subtypes=="MM"].
matrix=matrix[,subtypes=="T-ALL"|subtypes=="pre-B-ALL"]

# use 15% variable genes defined from whole hemap dataset to reproduce the same ALL map coordinates
top15var_genes = scan('DATA/useCase2/data9544_15pct_topVariating_genelist.txt',"a", quiet = T) # use 15% of the most variable genes, based on whole hemap dataset
matrix=matrix[rownames(matrix)%in%top15var_genes,]

# Load color vector of cytogenetics, mutations or clinical annotations to be used for plotting with Hemap cancer maps.
# When plotting different datasets, user has to naturally provide different color vectors to match the samples in the
# corresponding gene expression matrix.
# By default, color vector used here is the Hemap ALL cytogenetics colors
hemap_color_vector=scan("DATA/useCase3/Hemap_cytogenetics_colors.txt", "cytogenetics", quiet=TRUE)

#********************************************************************************************

print("New data and Hemap data loaded")
