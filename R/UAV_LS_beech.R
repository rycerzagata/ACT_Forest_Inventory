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
library(sp)
library(rgl)
library(rlas)
library(ForestTools)
library(TreeLS)
library(randomForest)
library(EBImage)
library(Metrics)
library(raster)
library(colorRamps)

readLAS<-lidR::readLAS

#Setting working directory
setwd("~/ACT/ACT_Forest_Inventory")

#### CHM COMPUTATION ####
set.seed(2020)

# Load and read the AHN3 file
AHN3_clip <- "Data/ps02_AHN3.laz"
AHN3 <- readLAS(AHN3_clip)

# Load and clip the Laz file
lasfile <- "Data/ps02_UAV_LS.laz"
beechLas <- readLAS(lasfile)
x <-  c(176254, 176185, 176167, 176236)
y <- c(473741, 473712, 473754, 473783)
AHN3 <- lasclipPolygon(AHN3, x, y, inside=TRUE)
beechLas <- lasclipPolygon(beechLas, x, y, inside=TRUE)
                           
# Compute the DSM with the AHN3 dataset
DSM <- grid_canopy(beechLas, res=1, p2r(0.2))
#plot(DSM, main="DSM", col=matlab.like2(50))

# Compute the DTM with the AHN3 dataset
DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
DTM[is.na(DTM)] <- 0
#plot(DTM, main="DTM", col=matlab.like2(50))

# Compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
CHM[is.na(CHM)] <- 0

# Use focal statistics to smoothen the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)
CHM[is.na(CHM)] <- 0


#### TREE SEGMENTATION - DALPONTE APPROACH ####
set.seed(2020)

# Treetops detection using an algorithm based on a local maximum filter.
f <- function(x) { x * 0.08 + 2 }
ttops <- tree_detection(CHM, lmf(f))  
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
dir.create( "Data/pr02_extracted_laz_beech")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("Data/pr02_extracted_laz_beech/tree", i, ".laz"))}


#### DBH PREDICTION ####
set.seed(2020)

# Create a dataframe out of the crown polygons with the chosen sample trees and introduce the DBH measured
# using the software Cloud Compare. There are some rules for choosing the right trees:
# no double trees on plot,the stem must be visible, no understory covering stems, returns distributed in cylindrical 
# shapes, trees distributed across a wide range of DBH (5-50 cm) and geographically distributed throughout the area.
sample_index <- c( 1, 4, 6, 8, 10, 15, 18, 20, 23, 25, 30, 38 )
training <- as.data.frame(crownsPoly[sample_index,])
names(training) <- c("treeID", "height", "crownArea", "crownDiameter")
training$DBH <- c( 0.61, 0.55, 0.4, 0.57, 0.74, 0.53, 0.54, 0.56, 0.7, 0.89, 0.32, 0.56 )

# Create the test dataset 
test <- as.data.frame(crownsPoly[-sample_index,])
names(test) <- c("treeID", "height", "crownArea", "crownDiameter")

# Train a Random Forest model to predict DBH
model <- randomForest(DBH ~ height + crownArea, 
                      data = training, importance= TRUE, proximity = TRUE, ntree = 500)

# Predict the DBH of the test dataset
test$DBH <- predict(model, test)

# Merge both datasets
full_dataset <- rbind(training, test)

# Compute the standing volume with the DBH and the height. Source: https://silvafennica.fi/pdf/smf004.pdf
full_dataset$standing_volume <- ((0.049)*((full_dataset$DBH*100)^1.78189)*(full_dataset$height)^1.08345)/1000

totalVolume <- sum(as.matrix(full_dataset$standing_volume))
emptyArea <- 200                            # area of empty spaces in the forest measured with polygons in ArcgIS/QGIS
totalArea <- raster::area(AHN3) - emptyArea
m3ha <- totalVolume/(totalArea/10000) 
m3ha

write.csv(full_dataset,"Data/pr02_UAV_LS_beech.csv", row.names = TRUE)


#### VALIDATION ####

# Mean squared error of the RF model
MSE <- mean(model$mse[1:500])
plot(model,  main="MSE of a RF model")

# Read the CSVs generated above and in the TLS scripts
TLS_dataset <- read.csv("Data/pr01_TLS_beech_valid.csv", header=TRUE, sep = ",")
UAV_dataset <- read.csv("Data/pr02_UAV_LS_beech.csv", header=TRUE, sep = ",")


# Compute some statistics of both datasets
trees_UAV <- nrow(UAV_dataset)
trees_TLS <- nrow(TLS_dataset)
trees_ha_UAV <- trees_UAV/(totalArea/10000)
trees_ha_TLS <- trees_TLS/(totalArea/10000)
mean_height_UAV <- mean(UAV_dataset$height)
mean_height_TLS <- mean(TLS_dataset$height)
mean_DBH_UAV <- mean(UAV_dataset$DBH)
mean_DBH_TLS <- mean(TLS_dataset$DBH)
m3ha_TLS <- sum(TLS_dataset$standing_volume)/(totalArea/10000)

# Compare if there is a significant difference between the DBH from TLS and UAV-LS
t_test_DBH <- t.test(TLS_dataset$DBH, UAV_dataset$DBH, 
                                   paired = FALSE, alternative = "two.sided")

mean_volume_UAV <- mean(UAV_dataset$standing_volume)
mean_volume_TLS <- mean(TLS_dataset$standing_volume)

# Compare if there is a significant difference between the standing volume from TLS and UAV-LS
t_test_volume <- t.test(TLS_dataset$standing_volume, UAV_dataset$standing_volume, 
                               paired = FALSE, alternative = "two.sided")


validation_results <- data.frame("Dataset" = c("UAV", "TLS"), "Number of trees" = c(trees_UAV, trees_TLS), 
                                 "Trees per ha" = c(trees_ha_UAV, trees_ha_TLS), "Mean height (m)" = c(mean_height_UAV, mean_height_TLS),
                                 "Mean DBH (m)" = c(mean_DBH_UAV, mean_DBH_TLS), "t-test DBH" = t_test_DBH$p.value,
                                 "CI DBH" =t_test_DBH$conf.int,"St error DBH" = t_test_DBH$stderr,
                                 "Mean volume (m3)" = c(mean_volume_UAV, mean_volume_TLS),"t-test volume" = t_test_volume$p.value,
                                 "CI volume"=t_test_volume$conf.int,"St error volume" = t_test_volume$stderr, "m3 per ha" = c(m3ha, m3ha_TLS))

# Compute the RMSE of DBH and standing volume
rmse(TLS_dataset$DBH, UAV_dataset$DBH)
rmse(TLS_dataset$standing_volume, UAV_dataset$standing_volume)
  
# Export the results in an excel
write.table(validation_results, "Data/pr02_UAV_LS_beech_valid.csv", row.names = TRUE)
