# This script details the import and QC filtering of single-cell
# data from male and female Ancylostoma ceylanicum samples that
# had been pre-processed using CellRanger. After QC steps this
# script will then cluster each sample and identify marker 
# genes per cluster, then show how these 2 samples were integrated
# to support visualizations of the data and assist in assigning
# cell type per cluster
# NOTE: In our publication we selected resolution 0.1 for both the female and male samples, which resulted in 10 clusters and 13 clusters respectively.
#       But due to the stochastic nature of how doubletFinder subsamples the cells (based on the pN setting, for which we used 0.25, ie. 25% of cells)
#       the QC filtered set of cells used for all downstream analysis can vary slightly. This difference in the available 'singlet' cells can result in
#       different clustering at the same resolutions. The only way to confidently reproduce exact results from the paper is to use the RDS files we 
#       provide!




# -- Install packages as needed --
# Versions listed here are the ones used in the preparation of 
# our recent A.ceylanicum single-cell paper
install.packages("devtools")
install.packages("RColorBrewer")
install.packages("remotes")
remotes::install_version("Matrix", version = "1.6-1")
install.packages("https://cran.rstudio.com/src/contrib/matrixStats_1.2.0.tar.gz", repos=NULL, type="source")
devtools::install_version("Seurat", version = "4.4.0")
install_github("chris-mcginnis-ucsf/DoubletFinder")
install.packages("clustree")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mgcv") # Do not update packages if asked
BiocManager::install("MAST") # Do not update packages if asked




# -- Load libraries --
library(devtools)
library(RColorBrewer)
library(remotes)
library(Matrix)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(cowplot)
library(tidyr)
library(DoubletFinder)
library(clustree)
library(ggplot2)
library(patchwork)
library(mgcv)
library(MAST)
library(future)




# -- Set working directory --
setwd("/Users/jmartin/Library/CloudStorage/Box-Box/MitrevaLab/2023_Suman/scRNA_A.ceylanicum _Adult/assets_for_June24_paper/example_R_work")




# -- Import data from CellRanger --
# The Read10X function we use to import the data into a seurat object
# needs to be given the location of the feature_filtered_bc_matrix
# folder which will be in the 'outs' subdirectory of your CellRanger
# output. The feature_filtered_bc_matrix folder itself will contain
# 3 files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
# Path provided must be visible within the R working directory 
# set above

# Import sample counts data
AceyFemale.data <- Read10X(data.dir = "Data/female/filtered_feature_bc_matrix/")
AceyMale.data <- Read10X(data.dir = "Data/male/filtered_feature_bc_matrix/")

# Build seurat objects from counts data
# Initial filters applied here:
#   min.cells = 3  :only genes detected in >= 3 cells are kept
#   min.features = 100  :only keep cells with >= 100 genes detected
# Note that the A.ceylanicum gene names used in the reference provided to
# CellRanger had underscores ('_'), but Seurat does not allow that and
# so replaces them with dashes ('-') and gives warnings. This is not
# a problem, but is important to remember for downstream analysis
AceyFemale <- CreateSeuratObject(counts = AceyFemale.data, project = "AceyFemale", min.cells = 3, min.features = 100)
AceyMale <- CreateSeuratObject(counts = AceyMale.data, project = "AceyMale", min.cells = 3, min.features = 100)
# Clean up counts data which is no longer needed
rm(AceyFemale.data,AceyMale.data)




# -- Normalize per-sample counts using SCTransform --
AceyFemale <- SCTransform(AceyFemale, verbose = FALSE)
AceyMale <- SCTransform(AceyMale, verbose = FALSE)
# You may see warning messages at this step:
#   "In theta.ml(y = y, mu = fit$fitted) : iteration limit reached"
# This warning is thrown when a gene has very few non-zero observations
# and can be safely ignored (see developer post: 
# https://github.com/satijalab/seurat/issues/1378)




# -- Use the doubletFinder tool to filter out double cells --
# This should be done on each sample individually before
# integration

# Generate PCA dimensionality reduction and UMAP projection
# By default RunPCA will use 50 PC, but for analysis work
# we use only 20 PC which still captures almost all the 
# variability while being less compute intensive
AceyFemale <- RunPCA(AceyFemale, verbose = FALSE)
AceyFemale <- RunUMAP(AceyFemale, dims = 1:20)
AceyMale <- RunPCA(AceyMale, verbose = FALSE)
AceyMale <- RunUMAP(AceyMale, dims = 1:20)

# Determine optimal values for pN and pK which are needed
# arguments for doubleFinder
# pN is the proportion of cells that will be sampled per 
# iteration expressed as a decimel (ie. pN 0.25 means that 
# 25% of cells will be used). We will use pN = 0.25 which
# is what is used in the tutorial at:
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# pK defines the PC neighborhood size used to compute pANN
# ('proportion of artificial k nearest neighbors). This 
# value is estimated as follows
sweep.female.res.list <- paramSweep(AceyFemale, PC = 1:20, sct = TRUE, num.cores = 1) # If you have multiple cores available that can speed up this step
sweep.female.stats <- summarizeSweep(sweep.female.res.list, GT = FALSE)
bcmvn_female <- find.pK(sweep.female.stats)
sweep.male.res.list <- paramSweep(AceyMale, PC = 1:20, sct = TRUE, num.cores = 1) # If you have multiple cores available that can speed up this step
sweep.male.stats <- summarizeSweep(sweep.male.res.list, GT = FALSE)
bcmvn_male <- find.pK(sweep.male.stats)
# These commands will print the sweep table for each sample
# to the console. Review this table and choose the pK value
# that corresponds to the highest BCmetric value
bcmvn_female
bcmvn_male
# Note that because the set of cells being sampled (assuming 
# you are using pN < 1) the optimal value (and thus optimal pK)
# can vary between runs. This means the final set of cells
# identified as doublets can vary between runs

# To estimate the expected number of homotypic doublets (doublets
# composed of cells from the same cluster) we do need to setup
# a baseline clustering resolution. This will not be the final
# resolution used for analysis. We can use the default resolution
# for the FindClusters function of seurat (0.8) for this
AceyFemale <- FindNeighbors(AceyFemale, dims = 1:20)
AceyFemale <- FindClusters(AceyFemale, dims = 1:20, save.SNN = TRUE) # default resolution is 0.8
annotations_female <- AceyFemale@meta.data$SCT_snn_res.0.8
homotypic.prop.female <- modelHomotypic(annotations_female)
nExp_poi_female <- round(0.075 * length(AceyFemale@meta.data$SCT_snn_res.0.8)) # Initial estimate: ~7.5% of cells are doublets
nExp_poi_female.adj <- round(nExp_poi_female * (1 - homotypic.prop.female)) # adjust estimate for homotypic doublets
AceyMale <- FindNeighbors(AceyMale, dims = 1:20)
AceyMale <- FindClusters(AceyMale, dims = 1:20, save.SNN = TRUE) # default resolution is 0.8
annotations_male <- AceyMale@meta.data$SCT_snn_res.0.8
homotypic.prop.male <- modelHomotypic(annotations_male)
nExp_poi_male <- round(0.075 * length(AceyMale@meta.data$SCT_snn_res.0.8)) # Initial estimate: ~7.5% of cells are doublets
nExp_poi_male.adj <- round(nExp_poi_male * (1 - homotypic.prop.male)) # adjust estimate for homotypic doublets
# Run doubletFinder
AceyFemale <- doubletFinder(AceyFemale, PCs = 1:20, pN = 0.25, pK = 0.03, nExp = nExp_poi_female, reuse.pANN = FALSE, sct = TRUE)
AceyMale <- doubletFinder(AceyMale, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi_male, reuse.pANN = FALSE, sct = TRUE)
# Run doubletFinder a 2nd time to refine the DF classifications
# Use these commands to find the meta.data column name of the pANN values
head(AceyFemale@meta.data, n = 2) # in this case: pANN_0.25_0.03_364
head(AceyMale@meta.data, n = 2) # in this case: pANN_0.25_0.005_322
AceyFemale <- doubletFinder(AceyFemale, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi_female.adj, reuse.pANN = "pANN_0.25_0.03_364", sct = TRUE)
AceyMale <- doubletFinder(AceyMale, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi_male.adj, reuse.pANN = "pANN_0.25_0.005_322", sct = TRUE)

# Visualize doublets identified in each sample (optional)
# Use these commands to find the meta.data column name of the DF classifications from the 2nd iteration of doubletFinder
head(AceyFemale@meta.data, n = 2) # in this case: DF.classifications_0.25_0.06_335
head(AceyMale@meta.data, n = 2) # in this case: DF.classifications_0.25_0.29_300
DimPlot(AceyFemale, reduction = "umap", group.by = "DF.classifications_0.25_0.06_335")
VlnPlot(AceyFemale, features = "nFeature_RNA", pt.size = 0, group.by = "DF.classifications_0.25_0.06_335")
DimPlot(AceyMale, reduction = "umap", group.by = "DF.classifications_0.25_0.29_300")
VlnPlot(AceyMale, features = "nFeature_RNA", pt.size = 0, group.by = "DF.classifications_0.25_0.29_300")

# Filter out cells that are not singlets (ie. keep only cells classified as 'Singlet')
AceyFemale <- subset(AceyFemale, subset = DF.classifications_0.25_0.06_335 == "Singlet")
AceyMale <- subset(AceyMale, subset = DF.classifications_0.25_0.29_300 == "Singlet")

# Apply SCTransform on the filtered cells
AceyFemale <- SCTransform(AceyFemale, verbose = FALSE)
AceyMale <- SCTransform(AceyMale, verbose = FALSE)




# -- Identify marker genes and choose cluster resolution --
# We will try a range of cluster resolutions with the goal of
# finding one that splits the cells into clusters that
# represent specific cell types. If you have a biological
# expectation going in you can aim that that number of clusters
# Otherwise you will need to look at cluster specific marker 
# genes for the clusters across various cluster resolutions 
# and choose that one that seems most reasonable
AceyFemale <- RunPCA(AceyFemale, verbose = FALSE)
AceyFemale <- FindNeighbors(AceyFemale, dims = 1:20)
AceyFemale <- FindClusters(AceyFemale, dims = 1:20, resolution = c(0.01,0.1,0.5,1)) # We used 0.1 in paper, but you will usually need to try multiple resolutions
AceyFemale <- RunUMAP(AceyFemale, dims = 1:20)
AceyMale <- RunPCA(AceyMale, verbose = FALSE)
AceyMale <- FindNeighbors(AceyMale, dims = 1:20)
AceyMale <- FindClusters(AceyMale, dims = 1:20, resolution = c(0.01,0.1,0.5,1)) # We used 0.1 in paper, but you will usually need to try multiple resolutions
AceyMale <- RunUMAP(AceyMale, dims = 1:20)

# Visualize results
DimPlot(AceyFemale, reduction = "umap", group.by = "SCT_snn_res.0.01")
DimPlot(AceyFemale, reduction = "umap", group.by = "SCT_snn_res.0.1") # For object built in paper this resolution generates 10 clusters
DimPlot(AceyFemale, reduction = "umap", group.by = "SCT_snn_res.0.5")
DimPlot(AceyFemale, reduction = "umap", group.by = "SCT_snn_res.1")
DimPlot(AceyMale, reduction = "umap", group.by = "SCT_snn_res.0.01")
DimPlot(AceyMale, reduction = "umap", group.by = "SCT_snn_res.0.1") # For object built in paper this resolution generates 13 clusters
DimPlot(AceyMale, reduction = "umap", group.by = "SCT_snn_res.0.5")
DimPlot(AceyMale, reduction = "umap", group.by = "SCT_snn_res.1")
# NOTE: In our publication we selected resolution 0.1 for both the female and male samples, which resulted in 10 clusters and 13 clusters respectively.
#       But due to the stochastic nature of how doubletFinder sub samples the cells based on the pN setting, for which we used 0.25 (ie. 25% of cells),
#       the QC filtered set of cells used for downstream analysis can vary slightly. This difference in the available 'singlet' cells can result in
#       different clustering at the same resolutions. The only way to confidently reproduce exact results from the paper is to use the RDS files we 
#       provide

# Run FindAllMarkers on the resolutions to guide the choice of the
# final resolution to use, and what cell types are present. Here I 
# show how to generate marker genes per cluster at the resolution
# we chose as the 'correct' resolution. But in practice you may
# need to generated lists of marker genes per cluster at all the
# resolutions you are considering
AceyFemale_markers_res.0.1 <- FindAllMarkers(AceyFemale, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", group.by = "SCT_snn_res.0.1", test.use = "MAST")
write.table(AceyFemale_markers_res.0.1, "AceyFemale_markers_res.0.1.tsv", row.names = FALSE)
AceyMale_markers_res.0.1 <- FindAllMarkers(AceyMale, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", group.by = "SCT_snn_res.0.1", test.use = "MAST")
write.table(AceyMale_markers_res.0.1, "AceyMale_markers_res.0.1.tsv", row.names = FALSE)

# Before integration, it is useful to add meta.data columns to
# label cells with our selected resolution cluster number and/or
# male/female status. This helps later when we want to visualize
# the location of these specific groupings of cells in the context
# of the integrated dimension reduction
AceyFemale@meta.data$sample_gender <- "Female"
AceyFemale@meta.data$FemaleClusterNumber <- AceyFemale@meta.data$SCT_snn_res.0.1
AceyFemale@meta.data$MaleClusterNumber <- "Female"
AceyMale@meta.data$sample_gender <- "Male"
AceyMale@meta.data$MaleClusterNumber <- AceyMale@meta.data$SCT_snn_res.0.1
AceyMale@meta.data$FemaleClusterNumber <- "Male"




# -- Generate table of cell-counts per cluster per sample --
# After assigning cell type to cluster based on marker genes
# it is useful to report the number of cells detected per
# cluster(cell type) per each sample
AceyFemale_CellCount_perCluster <- table(AceyFemale@meta.data$SCT_snn_res.0.1)
write.table(AceyFemale_CellCount_perCluster, "AceyFemale_CellCount_perCluster_res.0.1.tsv", row.names =  FALSE)
AceyMale_CellCount_perCluster <- table(AceyMale@meta.data$SCT_snn_res.0.1)
write.table(AceyMale_CellCount_perCluster, "AceyMale_CellCount_perCluster_res.0.1.tsv", row.names =  FALSE)




# -- Integrate male and female objects --
# At this point we have assigned cell type to each male and female
# cluster based on marker genes, and often want to visualize cell type
# from one sample on top of the other sample. To do this we need to
# generate a dimension reduction that incorporates both samples. We
# need to integrate the sample objects to enable these visualizations

# Set number cores to use in 'plan' options
num_cores <- 1

# Define function to integrate Seurat objects
integrate_samples_from_SeuratObj <- function(samples) {
  for(i in samples) {
    # Find most variable and informative genes for each sample
    assign(i, FindVariableFeatures(get(i)))
  }
  # Select integration features and prepare for SCT integration
  int.features <- SelectIntegrationFeatures(object.list = mget(samples), nFeatures = 2000)
  int.list <- PrepSCTIntegration(object.list = mget(samples), anchor.features = int.features)
  # Increase memory limit and request multiple processors for FindIntegrationAnchors
  options(future.globals.maxSize = 4800000000)
  plan("multisession", workers = num_cores)
  # Find integration anchors
  int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features = int.features)
  # Finding the overlapping gene set for integration
  subgenes <- Reduce(intersect, lapply(mget(samples), row.names))
  # Increase memory limit for IntegrateData
  options(future.globals.maxSize = 6400000000)
  # Integrate data and perform dimensionality reduction
  all <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", features.to.integrate = subgenes)
  all <- RunPCA(all) # Run PCA
  all <- RunTSNE(all, dims = 1:20) # Run t-SNE using 20 PCs
  all <- RunUMAP(all, dims = 1:20) # Run UMAP using 20 PCs
  return(all)
}

# Make sure both objects are using SCT assay
DefaultAssay(AceyFemale) <- "SCT"
DefaultAssay(AceyMale) <- "SCT"

# Integrate the objects
list_Acey <- c("AceyFemale", "AceyMale")
AllAcey <- integrate_samples_from_SeuratObj(list_Acey)




# -- Subcluster a specified group of cells --
# This section shows how to subcet a specific grouping
# of cells into its own object, and then clustering
# and normalizing that subcluster. This is often done
# to explore a cluster whose cell type is know to see
# if those cells can be further clustered in a meaningful
# way. 
# NOTE: This example code shows how we subclustered
#       the cluster we had defined as male hypodermis.
#       If you are running this code from the raw
#       data provided it is possible that the cluster
#       number we identified as male hypodermis will
#       be different than the cluster's you define
#       (due to how doubletFinder subsamples data, and
#       may result in different identifications of 
#       doublet cells that are removed). In order to
#       reproduce our work you will need to load the
#       pre-existing .rds file that is provided and
#       start this section using that object

# Build your seurat object using our saved .rds file
# (this should reproduce subclustering work described
# in our publication)
AllAcey <- readRDS(file = "Acey_Integrated.rds")

# We will subcluster the cluster we identified as being
# Hypodermis in the male sample. First we need to build
# a subset object containing only these cells. In our
# local work (which you would reproduce if you are using
# the Acey_Integrated.rds shown above) Hypodermis is
# cluster #4
AceyMaleHypo <- subset(AllAcey, sample_gender == "Male" & MaleClusterNumber == 4)

# Normalize the subset object using SCTransform
AceyMaleHypo <- SCTransform(AceyMaleHypo, verbose = FALSE)

# Subcluster this subset object, apply dimension reduction
AceyMaleHypo <- RunPCA(AceyMaleHypo, verbose = FALSE)
AceyMaleHypo <- FindNeighbors(AceyMaleHypo, dims = 1:20)
AceyMaleHypo <- FindClusters(AceyMaleHypo, dims = 1:20, resolution = c(0.05, 0.1, 0.2, 0.4), save.SNN = TRUE) # Trying several resolutions with hope of finding clusters with biological meaning
AceyMaleHypo <- RunUMAP(AceyMaleHypo, dims = 1:20)

# Visualize the male hypodermis subcluster resolutions
DimPlot(AceyMaleHypo, reduction = "umap", group.by = "SCT_snn_res.0.05") + ggtitle("Male Hypo SCT_0.05")
DimPlot(AceyMaleHypo, reduction = "umap", group.by = "SCT_snn_res.0.1") + ggtitle("Male Hypo SCT_0.1")
DimPlot(AceyMaleHypo, reduction = "umap", group.by = "SCT_snn_res.0.2") + ggtitle("Male Hypo SCT_0.2")
DimPlot(AceyMaleHypo, reduction = "umap", group.by = "SCT_snn_res.0.4") + ggtitle("Male Hypo SCT_0.4")
# On review of these plots (assuming you used the provided .rds) 
# both the 0.1 and 0.2 resolutions produce cluster groups that 
# appear meaningful in terms of the UMAP plot. 

# Once you've chosen the cluster resolution you prefer, it
# is helpful to identify marker genes per each of these 
# subclusters. In this example we choose resolution 0.2
AceyMaleHypo_markers_res.0.2 <- FindAllMarkers(AceyMaleHypo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", group.by = "SCT_snn_res.0.2", test.use = "MAST")
write.table(AceyMaleHypo_markers_res.0.2, "AceyMaleHypo_markers_res.0.2.tsv", row.names = FALSE)




# -- Visualize specific set of cells on integrated plot --
# In general, to visualize specific groupings of cells within
# the context of all cells after integration, you can add 
# meta.data columns holding the information you want to 
# highlight per cell, along with some background value
# (eg. setting the FemaleClusterNumber to 'Male' for cells
# that are male and thus have no female cluster number).
# You can add these additional meta.data columns to the 
# original objects before integration, and those columns
# should propagate into the integrated object

# This is an example of highlighting male cells in the 
# integrated object's umap dimension reduction ... to
# set colors of data points (cells) you can modify the
# scale_color_manual 'values' list. Each category being
# plotted, as defined by the group.by argument in the
# DimPlot command, has a pair of values. The first value
# in each pair should match a value in the meta.data slot
# used in the group.by argument, and the second value is
# the color you want
ex_plot_1 <- DimPlot(AllAcey, reduction = "umap", group.by = "sample_gender")
ex_plot_1 <- ex_plot_1 + scale_color_manual(values = c("Male" = "red", "Female" = "gray91"))
ex_plot_1

# This is an example of highlighting the selected cluster number
# membership of all female cells in the integrated object's tsne
# dimension reduction
ex_plot_2 <- DimPlot(AllAcey, reduction = "tsne", group.by = "FemaleClusterNumber")
ex_plot_2 <- ex_plot_2 + scale_color_manual(values = c("Male" = "gray91", "0" = "slateblue", "1" = "thistle1", "2" = "tomato4", "3" = "wheat3", "4" = "palegreen4", "5" = "paleturquoise2", "6" = "orange2", "7" = "violetred2")) # Note, while test running this script there were 8 female cell clusters at chosen resolution, but this can change when steps are re-run (due to doubletFinder subsampling)
ex_plot_2
