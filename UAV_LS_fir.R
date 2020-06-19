"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script can be used to compute the standing volume for Douglas fir from AHN3 and UAV-LS data.
At the end of the script some validation with TLS data is performed. 

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
AHN3_clip <- "AHN3.laz"
AHN3 <- readLAS(AHN3_clip)

# Load and clip the Laz file
lasfile <- "UAV_withGround.laz"
douglasLas <- readLAS(lasfile)
x <- c(176036, 176064, 176090, 176109,  176060, 176052)
y <- c(473695, 473725, 473723, 473679,  473657, 473657)
AHN3 <- lasclipPolygon(AHN3, x, y, inside=TRUE)
douglasLas <- lasclipPolygon(douglasLas, x, y, inside=TRUE)

# Compute the DSM with the AHN3 dataset
DSM <- grid_canopy(douglasLas, res=1, p2r(0.2))
#plot(DSM, main="DSM", col=matlab.like2(50))

# Compute the DTM with the AHN3 dataset
DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
#plot(DTM, main="DTM", col=matlab.like2(50))

# Compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
#CHM[is.na(CHM)] <- 0

# Use focal statistics to smoothen the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)

#### TREE SEGMENTATION - DALPONTE APPROACH ####
set.seed(2020)

# Treetops detection using an algorithm based on a local maximum filter.
f <- function(x) { x * 0.08 + 2 }
ttops <- tree_detection(CHM, lmf(f))  
#ttops <- tree_detection(douglasLas, lmf(5)) # you can do it using las file too but it takes some time
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
nlas <- lasnormalize(douglasLas, DTM)
vegpoints <- nlas %>% lasfilter(Classification==1) 

# Individual tree segmentation based on the Dalponte and Coomes (2016) algorithm.
# The returned point cloud has a new extra byte attribute named treeID.
trees <- lastrees(vegpoints, dalponte2016(chm = CHM, treetops = ttops, ID = 'treeID'))
#plot(trees, color="treeID")

# Extract every tree into a separate .laz file
dir.create( "extracted_laz_douglas")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("extracted_laz_douglas/tree", i, ".laz"))}


#### DBH PREDICTION ####
library(randomForest)
set.seed(2020)

# Create a dataframe out of the crown polygons with the chosen sample trees and introduce the DBH measured
# using the software Cloud Compare
sample_index <- c(3, 4, 9, 11, 15:17, 19, 23, 37:41, 45, 51, 52)
training <- as.data.frame(crownsPoly[sample_index,])
names(training) <- c("treeID", "height", "crownArea", "crownDiameter")
training$DBH <- c(0.42, 0.37, 0.32, 0.436, 0.32, 0.456, 0.43, 0.25, 0.41, 0.482, 0.455, 0.42, 0.41, 0.39, 0.44, 0.432, 0.41)

# Create the test dataset 
test <- as.data.frame(crownsPoly[-sample_index,])
names(test) <- c("treeID", "height", "crownArea", "crownDiameter")

# Train a Random Forest model to predict DBH
model <- randomForest(DBH ~ height + crownArea, 
                      data = training, importance= TRUE, proximity = TRUE, ntree=200)

# Predict the DBH of the test dataset
test$DBH <- predict(model, test)

# Merge both datasets
full_dataset <- rbind(training, test)

# Compute the standing volume with the DBH and the height. Source: https://silvafennica.fi/pdf/smf004.pdf
full_dataset$standing_volume <- (((full_dataset$DBH*100)^1.90053)*(full_dataset$height^0.80726)*exp(-2.43151))/1000

totalVolume <- sum(as.matrix(full_dataset$standing_volume))
emptyArea <- 230                    # area of empty spaces in the forest measured with polygons in ArcgIS/QGIS
totalArea <- raster::area(AHN3) - emptyArea
m3ha <- totalVolume/(totalArea/10000) 
m3ha

write.csv(full_dataset,"/Users/marariza/Downloads/UAV-LS-results_douglas.csv", row.names = TRUE)

#####################################################################################################
# VALIDATION

# Read the excels generated above and in the TLS scripts
TLS_dataset <- read.csv("TLS_fir.csv", header=TRUE, sep = ",")
UAV_LS_dataset <- read.csv("UAV-LS-results_douglas.csv", header=TRUE, sep = ",")

# Compute some statistics of both datasets
trees_UAV <- nrow(UAV_LS_dataset)
trees_TLS <- nrow(TLS_dataset)
trees_ha_UAV <- trees_UAV/(totalArea/10000)
trees_ha_TLS <- trees_TLS/(totalArea/10000)
mean_height_UAV <- mean(UAV_LS_dataset$height)
mean_height_TLS <- mean(TLS_dataset$height)
mean_DBH_UAV <- mean(UAV_LS_dataset$DBH)
mean_DBH_TLS <- mean(TLS_dataset$DBH)
m3ha_TLS <- sum(TLS_dataset$standing_volume)/(totalArea/10000)

# Compare if there is a significant difference between the DBH from TLS and UAV-LS
t_test_DBH <- t.test(TLS_dataset$DBH, UAV_LS_dataset$DBH, 
                     paired = FALSE, alternative = "two.sided")

mean_volume_UAV <- mean(UAV_LS_dataset$standing_volume)
mean_volume_TLS <- mean(TLS_dataset$standing_volume)

# Compare if there is a significant difference between the standing volume from TLS and UAV-LS
t_test_volume <- t.test(TLS_dataset$standing_volume, UAV_LS_dataset$standing_volume, 
                        paired = FALSE, alternative = "two.sided")


validation_results <- data.frame("Dataset" = c("UAV", "TLS"), "Number of trees" = c(trees_UAV, trees_TLS), 
                                 "Trees per ha" = c(trees_ha_UAV, trees_ha_TLS), "Mean height" = c(mean_height_UAV, mean_height_TLS),
                                 "Mean DBH" = c(mean_DBH_UAV, mean_DBH_TLS), "t-test DBH" = t_test_DBH$p.value,
                                 "CI DBH" =t_test_DBH$conf.int,"St error DBH" = t_test_DBH$stderr,
                                 "Mean volume" = c(mean_volume_UAV, mean_volume_TLS),"t-test volume" = t_test_volume$p.value,
                                 "CI volume"=t_test_volume$conf.int,"St error volume" = t_test_volume$stderr, "m^3 per ha" = c(m3ha, m3ha_TLS))


# Export the results in an excel
write.table(validation_results, "/Users/marariza/Downloads/validation_results_LS_fir.csv", row.names = TRUE)




