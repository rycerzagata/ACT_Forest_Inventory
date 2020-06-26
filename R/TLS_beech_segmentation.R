"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create tree segmentation fitting using DALPONTE APPROACH.
The script should be run for each LAZ tile in the zip file.
"""

# Loading the required libraries
library(lidR)
library(raster)
library(sp)
library(rgl)
library(TreeLS)
library(rLiDAR)


#### CHM COMPUTATION ####
## Setting working directory
setwd("~/ACT/ACT_Forest_Inventory")

# We load and read the AHN3 file
zipfile <- 'Data/ps01_TLS_beech.zip'
outdir <- 'Data'
unzip(zipfile,exdir=outdir)

lasfile <- "Data/ps01_TLS_beech.laz"
beechLas <- readTLS(lasfile)

# Computing the DSM with the TLS dataset
DSM <- grid_canopy(beechLas, res=1, p2r(0.2))

# Computing the DTM with the TLS dataset
DTM <- grid_terrain(beechLas, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)

# We compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
CHM[is.na(CHM)] <- 0

# Using focal statistics to smooth the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)


#### TREE SEGMENTATION - DALPONTE APPROACH ####

# Treetops detection using an algorithm based on a local maximum filter. 
# The function can be point cloud based but due to long processing time we use CHM.
ttops <- tree_detection(CHM, lmf(4, 2))

# Normalize point cloud
nlas <- lasnormalize(beechLas, DTM)

# Individual tree segmentation based on the Dalponte and Coomes (2016) algorithm.
# The returned point cloud has a new extra byte attribute named treeID.
trees <- lastrees(nlas, dalponte2016(CHM, ttops))
plot(trees, color="treeID") 

# We extract every tree into a different .laz file
dir.create( "Data/pr01_extracted_laz_TLS_beech")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("Data/pr01_extracted_laz_TLS_beech/tree", i, ".laz"))}





