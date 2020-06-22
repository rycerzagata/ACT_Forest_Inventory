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
library(lidR)
library(colorRamps)
library(ggpubr)
library(rlas)
library(tiff)
library(ForestTools)
library(itcSegment)
library(TreeLS)

## Setting working directory
setwd("D:/01. Kuliah/02. Materi/Period 6/process/R_script")

## Douglas Fir species has a large size so it must divides into 6 files
## create folder and unzip data
zipfile <- 'D:/01. Kuliah/02. Materi/Period 6/process/R_script/Data/douglas_fir_data.zip'
outdir <- 'D:/01. Kuliah/02. Materi/Period 6/process/R_script/Data/douglas'
unzip(zipfile,exdir=outdir)

## load each data to process until segmentation
## read laz file of TLS data
tls_file <- "Data/douglas/tls_dfir_ground_cliped_agata_634_159.laz"
tls <- readTLS(tls_file)

## clip into the intersection part between UAV LiDAR and TLS
x <- c(176036, 176064, 176090, 176109,  176060, 176052)
y <- c(473695, 473725, 473723, 473679,  473657, 473657)

douglas <- lasclipPolygon(tls, x, y, inside=TRUE)

# Computing the DSM with the TLS dataset
DSM <- grid_canopy(douglas, res=1, p2r(0.2))


# Computing the DTM with the TLS dataset
DTM <- grid_terrain(douglas, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)

# We compute the CHM and remove one value which is below 0 (-0.005 m)
chm <- DSM - DTM
chm[is.na(chm)] <- 0

# Using focal statistics to smooth the CHM
chm <- focal(chm,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)

#################################################
#### TREE SEGMENTATION - DALPONTE APPROACH ####
ttops <- tree_detection(chm, lmf(4, 2))
crowns <- mcws(treetops = ttops, CHM=chm, minHeight = 15, verbose=FALSE)
crownsPoly <- mcws(treetops = ttops, CHM=chm, minHeight = 8, verbose=FALSE, format="polygons")

# Normalize the point cloud using DTM
nlas <- lasnormalize(douglas, DTM)

trees <- lastrees(nlas, dalponte2016(chm, ttops))
plot(trees, color="treeID") 

# Extract every tree into a separate .laz file
# Give a different folder name for each dataset
dir.create( "extracted_laz/clipped_634_159")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("extracted_laz/clipped_634_159/tree", i, ".laz"))}




