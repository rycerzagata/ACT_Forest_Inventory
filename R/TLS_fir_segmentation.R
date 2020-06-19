"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create the CHM using TLS data. Derivates like treetops 
are computed. The last part of the script is tree segmentation using Dalponte approach.
"""

# Loading the required libraries
library(raster)
library(sp)
library(rgl)
library(lidR)
library(xlsx)
library(lidR)
library(colorRamps)
library(ggpubr)
library(rlas)
library(tiff)
library(ForestTools)
library(itcSegment)
library(TreeLS)

## Setting working directory


## read laz file of TLS data
tls_file <- "Data/LS_x=45_y=45_TLS_ground.laz"
tls <- readTLS(tls_file)

# Computing the DSM with the TLS dataset
DSM <- grid_canopy(tls, res=1, p2r(0.2))

# Computing the DTM with the TLS dataset
DTM <- grid_terrain(tls, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)

# We compute the CHM and remove one value which is below 0 (-0.005 m)
chm <- DSM - DTM
chm[is.na(chm)] <- 0


#################################################
#### TREE SEGMENTATION - DALPONTE APPROACH ####
ttops <- tree_detection(chm, lmf(4, 2))
crowns <- mcws(treetops = ttops, CHM=chm, minHeight = 15, verbose=FALSE)
crownsPoly <- mcws(treetops = ttops, CHM=chm, minHeight = 8, verbose=FALSE, format="polygons")

nlas <- lasnormalize(tls, DTM)
Vegpoints_norm <- nlas %>% lasfilter(Classification==1) 
trees <- lastrees(nlas, dalponte2016(chm, ttops))
plot(trees, color="treeID") 

# We extract every tree into a different .laz file
dir.create( "extracted_laz/x=45_y=45")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("extracted_laz/x=45_y=45/tree", i, ".laz"))}


