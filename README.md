# Scripts for t-SNE map generation, evaluation and sample remapping
Detailed User Guide: link. See use case examples for more information on usage in example dataset. 

- 1) To generate a t-SNE map from the required objects use the provided function CancerMap in **Cancermap.R** 

- 2) Quality control using **infoMeasures.R**. Evaluation of method and parameter choices. At this step it is wise to verify that the separation of the samples (e.g. clusters in a t-SNE map) reflects biology and not the data provider. We recommend to use established molecular subtypes for the disease in question. 

- 3) Remapping new samples to existing t-SNE map **run_remapping.R**. 

- 4) Plotting cancermaps and any color vector to t-SNE map using **Plot_cancermap_clusters.R, Plot_color_vector.R**

# Use Case Examples: t-SNE map optimization, remapping and data integration

![alt text](https://bioinformatics.uef.fi/~ppolonen/git_images/Cancermap_git.png)

More information: 
Mehtonen & Polonen et al. https://doi.org/10.1101/248096 

Use case 1: Find optimal parameters to use for t-SNE projection

Use case 2: Remap new samples to existing t-SNE map

Use case 3: Comparing an independent dataset with Hemap samples using t-SNE maps and gene set analysis

### Download example data:
The normalized data matrix is available for download from the web resource for further integrative analysis using standard R/Bioconductor software. The required objects are available as RData files and facilitate joint analysis of new datasets.
**Download data under cloned git folder Example_useCases/DATA**. This way, wrapper scripts work without modifying paths to read datasets in R.

https://drive.google.com/open?id=18i6EPHNDJyrAfZTLadKWfHQt3NxNd_za


### How to run example use cases:
Execute in R:
    
```R
source("Evaluate_cancermap_NMI.R") # Use Case 1
source("Remapping2Cancermap.R") # Use Case 2
source("Dataset_cancermap_comparison.R") # Use Case 3
```

### Main Inputs:
- **Use Case 1:** 
 * Hemap AML data set
 * Hemap AML data set id (GSEid) per sample

- **Use Case 2:** 
 * Hemap ALL data set
 * Hemap ALL t-SNE map coordinates
 * Remapping data set

- **Use Case 3:** 
 * Hemap ALL data set
 * Ross ALL data set
 * Cytogenetic annotations for both data sets

### Main Outputs:
- **Use Case 1:** 
 * NMI plot to find optimal parameters for t-SNE projection

- **Use Case 2:** 
 * Hemap ALL data set visualized together with new remapped samples
 * Coordinates for new samples

- **Use Case 3:** 
 * Visualize Hemap ALL clusters in ROSS ALL map (using gene set scores)
 * Visualize cytogenetic data in both datasets


### Libraries:

**R libraries:**
These R libraries are needed, for easy install, execute in R:
    
```R
# Packages, if not found, install automatically.
list.of.packages <- c("Rtsne", "LPCM", "foreach", "ggplot2", "data.table", "uuid", "reshape2", "gridExtra", "RColorBrewer", "doParallel", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

if(!"GSVA" %in% installed.packages()[,"Package"]) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GSVA")
}
```

**Remapping Algorithm**

cblas (and g++) is required!

First you need to compile the remapping program linking against cblas. Proceed to useCase2-directory and compile as follows. "Usually" the following will work:
```shell
g++ quadtree.cpp bh_tsne_remapping.cpp -o bh_tsne_remapping -O3 -I/usr/include -L/usr/lib64/atlas -lcblas
```
For CentOS/Fedora/RedHat use '-lsatlas' instead of '-lcblas'
```shell
g++ quadtree.cpp bh_tsne_remapping.cpp -o bh_tsne_remapping -O3 -I/usr/include -L/lib64/atlas -lsatlas
```
If you are using a conda environment, you can install f.e. openBLAS easily
```shell
conda install -c conda-forge openblas
```
Then compile
```shell
g++ quadtree.cpp bh_tsne_remapping.cpp -o bh_tsne_remapping -O3 -I/path/to/conda/env/include -L/path/to/conda/env/lib -lcblas
```
