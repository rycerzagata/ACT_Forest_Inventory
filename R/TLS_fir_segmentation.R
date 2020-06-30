"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create the CHM using TLS data. Derivates like treetops 
are computed. The last part of the script is tree segmentation using Dalponte approach.
"""

# Loading the required libraries
library(lidR)
library(raster)
library(sp)
library(rgl)
library(TreeLS)
library(rLiDAR)


## Setting working directory
setwd("~/ACT/ACT_Forest_Inventory")

## Douglas Fir species has a large size so it must be divided into 6 files
## create folder and unzip data
zipfile <- 'Data/ps01_TLS_douglas_fir.zip'
outdir <- 'Data/ps01_TLS_douglas_fir'
unzip(zipfile,exdir=outdir)

## load each data to process until segmentation
## read laz file of TLS data
tls_file <- "Data/ps01_TLS_douglas_fir/tls_fir_634_159.laz"
douglas <- readTLS(tls_file)

# Computing the DSM with the TLS dataset
DSM <- grid_canopy(douglas, res=1, p2r(0.2))


# Computing the DTM with the TLS dataset
DTM <- grid_terrain(douglas, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)

# We compute the CHM and remove one value which is below 0 (-0.005 m)
chm <- DSM - DTM
chm[is.na(chm)] <- 0

# Using focal statistics to smooth the CHM
chm <- focal(chm,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)


#### TREE SEGMENTATION - DALPONTE APPROACH ####

# Treetops detection using an algorithm based on a local maximum filter. 
# The function can be point cloud based but due to long processing time we use CHM.
ttops <- tree_detection(chm, lmf(4, 2))

# Normalize the point cloud using DTM
nlas <- lasnormalize(douglas, DTM)

# Individual tree segmentation based on the Dalponte and Coomes (2016) algorithm.
# The returned point cloud has a new extra byte attribute named treeID.
trees <- lastrees(nlas, dalponte2016(chm, ttops))
#plot(trees, color="treeID") 

# Extract every tree into a separate .laz file
# Give a different folder name for each dataset
dir.create( "Data/pr01_extracted_laz_TLS_fir")
dir.create( "Data/pr01_extracted_laz_TLS_fir/clipped_634_159")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("Data/pr01_extracted_laz_TLS_fir/clipped_634_159/tree", i, ".laz"))}




