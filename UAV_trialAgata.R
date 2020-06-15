
"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create the CHM using Lidar data from UAV and AHN3 datasets. Derivates like treetops 
and crown area are computed. The last part of the script is tree segmentation and circumference fitting using RANSAC algorithm.
"""

# Loading the required libraries
library(rLiDAR)
library(lidR)
library(raster)
library(colorRamps)
library(sp)
library(rgl)
library(ggpubr)
library(rlas)
library(tiff)
library(ForestTools)
library(TreeLS)
library(EBImage)

readLAS<-lidR::readLAS
#### CHM COMPUTATION ####
## Setting working directory
#setwd("/Users/marariza/Downloads")

# We load and read the AHN3 file
AHN3_clip <- "/Users/HP/Documents/ACT/R/Data/AHN3_beech.laz"
AHN3 <- readLAS(AHN3_clip)

lasfile <- "/Users/HP/Documents/ACT/R/Data/UAV_withGround.laz"
beechLas <- readLAS(lasfile)
beechLas <- lasclipRectangle(beechLas, 176170, 473657, 176265, 473782)

# Computing the DSM with the AHN3 dataset
DSM <- grid_canopy(beechLas, res=1, p2r(0.2))
#plot(DSM, main="DSM", col=matlab.like2(50))

# Computing the DTM with the AHN3 dataset
DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
#plot(DTM, main="DTM", col=matlab.like2(50))

# We compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
CHM[is.na(CHM)] <- 0

# Using focal statistics to smooth the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)

#################################################
############## Alternative way of smoothing CHM
sCHM<-CHMsmoothing(CHM, ws=3)

fws<-3 # dimention 3x3
# Set the specified height above ground for the detection break
minht<-15
# Create the individual tree detection list and summarize it
loc<-FindTreesCHM(sCHM, fws, minht)
loc
maxcrown=10.0
# Set the exclusion parameter - A single value from 0 to 1 that represents the % of pixel
# exclusion. E.g. a value of 0.5 will exclude all of the pixels for a single tree that has
#a height value of less than 50% of the maximum height from the same tree. Default value is 0.3.
exclusion=0.1
# Compute canopy areas for the individual tree detections
canopy<-ForestCAS(sCHM, loc, maxcrown, exclusion)
boundaryTrees<-canopy[[1]]
canopyList<-canopy[[2]]
XY<-SpatialPoints(canopyList[,1:2])
#################################################

#### BASIC APPROACH - C CHM DERIVATES ####

# We use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After we plot it to check how the tree tops look like. 
lin <- function(x) {x* 0.06 + 1}
treetops <- vwf(CHM = CHM, winFun = lin, minHeight = 10)
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(treetops, col="black", pch = 20, cex=0.5, add=TRUE)

# We check the mean of the height of the detected tree tops 
mean(treetops$height)

# We compute the function MCWS function that implements the watershed algorithm. In this case, the argument
# minHeight refers to the lowest expected treetop. The result is a raster where each tree crown is 
# a unique cell value. 
crowns <- mcws(treetops = treetops, CHM=CHM, minHeight = 15, verbose=FALSE)
#plot(crowns, main="Detected tree crowns", col=sample(rainbow(50), length(unique(crowns[])),replace=TRUE), 
     #legend=FALSE, xaxt="n", yaxt="n")

# We do the same computation as before but changig the output format to polygons. It takes more processing
# time but polygons inherit the attributes of treetops as height. Also, crown area is computed for each polygon.
crownsPoly <- mcws(treetops = treetops, CHM=CHM, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2
#mean(crownsPoly$crownDiameter)
#mean(crownsPoly$crownArea)

sp_summarise(treetops)
sp_summarise(crownsPoly, variables=c("crownArea", "height"))


#### STEM SEGMENTATION ####
tls = tlsNormalize(beechLas)
# map the trees on a resampled point cloud so all trees have approximately the same point density
thin = tlsSample(tls, voxelize(0.01))
map = treeMap(thin, map.hough(hmin = 1, hmax = 2, max_radius = 0.3, min_density = 0.1, min_votes = 2))
tls = stemPoints(tls, map)
df = stemSegmentation(tls, sgmt.ransac.circle(n=10))
head(df)
tlsPlot(tls, df, map)


#### TREE SEGMENTATION - DALPONTE APPROACH ####

# Treetops detection
f <- function(x) { x * 0.07 + 3 }
ttops <- tree_detection(CHM, lmf(f))  # or lmf(4, 2)
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(ttops, col="black", pch = 20, cex=0.5, add=TRUE)

# Crowns
crowns <- mcws(treetops = ttops, CHM=CHM, minHeight = 15, verbose=FALSE)
crownsPoly <- mcws(treetops = ttops, CHM=CHM, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Calculate diameter of crowns and add to crown data
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2

# Normalize the point cloud using DTM
nlas <- lasnormalize(beechLas, DTM)
Vegpoints_norm <- nlas %>% lasfilter(Classification==1) 
trees <- lastrees(nlas, dalponte2016(CHM, ttops))
plot(trees, color="treeID") 

# We extract every tree into a different .laz file
dir.create( "/Users/HP/Documents/ACT/R/Data/extracted_trees")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("/Users/HP/Documents/ACT/R/Data/extracted_trees/tree", i, ".laz"))}

file1 <- "/Users/HP/Documents/ACT/R/Data/extracted_trees/tree 101 .laz"
tree <- readLAS(file1)
plot(tree)

# calculate tree metrics
metrics <- tree_metrics(trees, func = ~max(Z), attribute = "treeID")
