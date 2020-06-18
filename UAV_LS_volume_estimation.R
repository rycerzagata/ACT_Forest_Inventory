
"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script is made in order to create the CHM using Lidar data from UAV and AHN3 datasets. Derivates like treetops 
and crown area are computed. Next part of the script is tree segmentation and prediction of DBH using Random Forest 
algorithm. The prediction of tree volume is computed based on the model and results are validated.
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

#Setting working directory
setwd("/Users/marariza/Downloads")

#### CHM COMPUTATION ####
set.seed(2020)

# Load and read the AHN3 file
AHN3_clip <- "AHN3_beech.laz"
AHN3 <- readLAS(AHN3_clip)

# Load and clip the Laz file
lasfile <- "UAV_withGround.laz"
beechLas <- readLAS(lasfile)
beechLas <- lasclipRectangle(beechLas, 176170, 473657, 176265, 473782)

# Compute the DSM with the AHN3 dataset
DSM <- grid_canopy(beechLas, res=1, p2r(0.2))
#plot(DSM, main="DSM", col=matlab.like2(50))

# Compute the DTM with the AHN3 dataset
DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
#plot(DTM, main="DTM", col=matlab.like2(50))

# Compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
CHM[is.na(CHM)] <- 0

# Use focal statistics to smoothen the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)


#### TREE SEGMENTATION - DALPONTE APPROACH ####
set.seed(2020)

# Treetops detection using an algorithm based on a local maximum filter.
f <- function(x) { x * 0.08 + 2 }
ttops <- tree_detection(CHM, lmf(f))  # or lmf(4, 2) for las
#ttops <- tree_detection(beechLas, lmf(5)) # you can do it using las file too but it takes some time
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(ttops, col="black", pch = 20, cex=0.5, add=TRUE)

# Crowns detection using MCWS function that implements the watershed algorithm to produce a map of crowns as polygons.. 
# In this case, the argument minHeight refers to the lowest expected treetop.
crownsPoly <- mcws(treetops = ttops, CHM=CHM, minHeight = 15, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) * 2

# Normalize the point cloud using DTM and select the vegetation class
nlas <- lasnormalize(beechLas, DTM)
vegpoints <- nlas %>% lasfilter(Classification==1) 

# Individual tree segmentation based on the Dalponte and Coomes (2016) algorithm.
# The returned point cloud has a new extra byte attribute named treeID.
trees <- lastrees(vegpoints, dalponte2016(chm = CHM, treetops = ttops, ID = 'treeID'))
#plot(trees, color="treeID")

# Extract every tree into a separate .laz file
dir.create( "extracted_laz")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("extracted_laz/tree", i, ".laz"))}

blablabla

#### DBH PREDICTION ####
library(randomForest)
set.seed(2020)

# Create a dataframe out of the crown polygons with the chosen sample trees and introduce the DBH measured
# using the software Cloud Compare
# no double trees on plot,the stem must be visible, no understory covering stems, returns distributed in cylindrical shapes
sample_index <- c(9:10, 16:17, 19, 21, 23, 26:45, 50:51, 53, 55, 68:69, 75, 79, 82, 84, 88, 91, 95, 96, 99)
training <- as.data.frame(crownsPoly[sample_index,])
names(training) <- c("treeID", "height", "crownArea", "crownDiameter")
training$DBH <- c(0.426, 0.45, 0.402, 0.594, 0.391, 0.814, 0.507, 0.391, 0.387, 0.359, 0.487, 0.511, 
                  0.377, 0.362, 0.328, 0.412, 0.427, 0.35, 0.447, 0.389, 0.325, 0.442, 0.389, 0.428, 
                  0.432, 0.39, 0.481, 0.376, 0.505, 0.525, 0.454, 0.47, 0.359, 0.514, 0.516, 0.461, 
                  0.308, 0.74, 0.479, 0.411, 0.31, 0.252)

# Create the test dataset 
test <- as.data.frame(crownsPoly[-sample_index,])
names(test) <- c("treeID", "height", "crownArea", "crownDiameter")

# Train a Random Forest model to predict DBH
model <- randomForest(DBH ~ height + crownArea, 
                      data = training, importance= TRUE, proximity = TRUE, ntree=500)

# Predict the DBH of the test dataset
test$DBH <- predict(model, test)

# Merge both datasets
full_dataset <- rbind(training, test)

# Compute the standing volume with the DBH and the height. Source: https://silvafennica.fi/pdf/smf004.pdf
full_dataset$standing_volume <- ((0.049/100)*(full_dataset$DBH^1.78189)*(full_dataset$height)^1.08345)*1000

totalVolume <- sum(as.matrix(full_dataset$standing_volume))
emptyArea <- 240                                     # area of empty spaces in the forest
totalArea <- raster::area(AHN3) - emptyArea
m3ha <- totalVolume/(totalArea/10000) 
m3ha


write.csv(full_dataset,"/Users/marariza/Downloads/UAV-LS-results.csv", row.names = TRUE)


#### VALIDATION ####

# Mean squared error of the RF model
MSE <- mean(model$mse[1:500])

TLS_dataset <- read.csv("TLS_beech.csv", header=TRUE, sep = ",")
UAV_LS_dataset <- read.csv("UAV-LS-results.csv", header=TRUE, sep = ",")

# Compare if there is a significant difference between the DBH from TLS and UAV-LS
comparison_RGB_LiDAR_DBH <- t.test(TLS_dataset$DBH, UAV_LS_dataset$DBH, 
                                   paired = FALSE, alternative = "two.sided")

# Compare if there is a significant difference between the standing volume from TLS and UAV-LS
comparison_RGB_LiDAR <- t.test(TLS_dataset$standing_volume, UAV_LS_dataset$standing_volume, 
                               paired = FALSE, alternative = "two.sided")




