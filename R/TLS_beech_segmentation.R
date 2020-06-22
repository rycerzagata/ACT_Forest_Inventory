"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create tree segmentation fitting using DALPONTE APPROACH.
"""

# Loading the required libraries
library(lidR)
library(raster)
library(rgl)
library(rlas)
library(ForestTools)
library(itcSegment)
library(TreeLS)
library(EBImage)
library(rLiDAR)


#### CHM COMPUTATION ####
## Setting working directory
setwd("../ACT_Forest_Inventory")

# We load and read the AHN3 file
lasfile <- "Data/tls_beech_ground_aligned.laz"
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

#################################################
#### TREE SEGMENTATION - DALPONTE APPROACH ####
ttops <- tree_detection(CHM, lmf(4, 2))
crowns <- mcws(treetops = ttops, CHM=CHM, minHeight = 15, verbose=FALSE)
crownsPoly <- mcws(treetops = ttops, CHM=CHM, minHeight = 8, verbose=FALSE, format="polygons")

nlas <- lasnormalize(beechLas, DTM)
Vegpoints_norm <- nlas %>% lasfilter(Classification==1) 
trees <- lastrees(nlas, dalponte2016(CHM, ttops))
plot(trees, color="treeID") 

# We extract every tree into a different .laz file
dir.create( "Data/extracted_laz_tls_beech")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("Data/extracted_laz_tls_beech/tree", i, ".laz"))}





