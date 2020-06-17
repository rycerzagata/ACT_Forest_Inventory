"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
Description: This script can be used to compute a canopy height model from AHN3 and UAV RGB data and the CHM derivates. 
These are later used to estimate the DBH and volume of trees. The results are validated with TLS data.
"""

# Loading the required libraries
library(lidR)
library(raster)
library(colorRamps)
library(sp)
library(rgl)
library(ggpubr)
library(rlas)
library(tiff)
library(itcSegment)
library(ForestTools)

readLAS<-lidR::readLAS

## Setting working directory
setwd("~/ACT_Forest_Inventory")

#### CHM COMPUTATION ####
# We load and read the beech AHN3 file
AHN3beech <- "Data/AHN3_beech.laz"
AHN3_beech <- readLAS(AHN3beech)

# We compute two DTMs modifying just one parameter (keep lowest)
DTM1 <- grid_terrain(AHN3_beech, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
DTM2 <- grid_terrain(AHN3_beech, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest=TRUE)

# We plot both DTMs to observe the differences
plot(DTM1, main="DTM1", col=matlab.like2(50))
plot(DTM2, main="DTM2", col=matlab.like2(50))

# We compute the differences between DTM and create and histogram to check if they are close to 0
Difference <- DTM2 - DTM1
plot(Difference, main="Difference", col=matlab.like2(50))
hist(Difference)
# As the values are almost the same (all around 0 in the difference), we will work with the first DTM (DTM1)

# We load the DSM created with photogrammetry from RGB images, we project it to RD New and resample to 
# have the same resolution as the DTM
RGB <- "Data/DEM_speulderbos_georef_ar.tif"
DSM_no_proj <- raster(RGB)
DSM_proj <- projectRaster(DSM_no_proj, crs = crs(DTM1))
DSM <- resample(DSM_proj, DTM1, method = "bilinear")

# Compute the CHM substracting also 40.68, which is the difference height between the DSM and the DTM
# in the area on the left (U-shaped area) and plot if needed.
CHM <- DSM - DTM1 - 40.68
# plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")

#Remove NA values
CHM[is.na(CHM)] <- 0

# Use the focal statistics to smooth the CHM
CHM_smooth <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)
plot(CHM_smooth)

#### CHM DERIVATES COMPUTATION ####
# Use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After this we plot it to check how the tree tops look like. 
lin <- function(x) { x * 0.02 + 0.5 }
treetops <- vwf(CHM = CHM_smooth, winFun = lin, minHeight = 20)
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(treetops, col="black", pch = 20, cex=0.5, add=TRUE)

# Check the mean of the height of the detected tree tops 
mean(treetops$height)

# Use the function MCWS function that implements the watershed algorithm to produce a map of crowns. 
# In this case, the argument minHeight refers to the lowest expected treetop. The result is a raster 
# where each tree crown is a unique cell value. 
crowns <- mcws(treetops = treetops, CHM=CHM_smooth, minHeight = 15, verbose=FALSE)
plot(crowns, main="Detected tree crowns", col=sample(rainbow(50), length(unique(crowns[])),replace=TRUE), 
     legend=FALSE, xaxt="n", yaxt="n")

# Use the  the MCWS function that implements the watershed algorithm to produce a map of crowns as polygons. It takes more processing
# time but polygons inherit the attributes of treetops as height. Also, crown area is computed for each polygon.
crownsPoly <- mcws(treetops = treetops, CHM=CHM_smooth, minHeight = 15, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2
mean(crownsPoly$crownDiameter)
mean(crownsPoly$crownArea)

sp_summarise(treetops)
sp_summarise(crownsPoly, variables=c("crownArea", "height"))

#### DBH AND TREE VOLUME ESTIMATION ####
# Compute the DBH based on the tree height in m and crown diameter in m^2, adapted to the type of biome.
crownsPoly$DBH <- dbh(H=crownsPoly$height, CA = crownsPoly$crownDiameter, biome=20)/100

# Compute the standing volume in m^3 from the DBH in cm and the height in m. Source: https://silvafennica.fi/pdf/smf004.pdf
crownsPoly$standing_volume <- (0.049)*(crownsPoly$DBH^1.78189)*(crownsPoly$height)^1.08345
hist(crownsPoly$standing_volume)

# Save the results to a csv file.
dataset <- as.data.frame(crownsPoly)
write.csv(dataset,"./Data/photogrammetry.csv", row.names = TRUE)

# Compute the total tree volume in m^3
totalVolume <- sum(as.matrix(crownsPoly$standing_volume))
emptyArea <- 873                                       # area of empty spaces in the forest
totalArea <- raster::area(AHN3_beech) - emptyArea      # area of forest
m3ha <- totalVolume/(totalArea/10000)                  # total tree volume in m^3 per hectare
m3ha

#### VALIDATION USING LiDAR DATA ####
photogrammetry_dataset <- read.csv("photogrammetry.csv", header=TRUE, sep = ",")
LiDAR_dataset <- read.csv("LiDAR.csv", header=TRUE, sep = ",")

comparison_RGB_LiDAR <- t.test(photogrammetry_dataset$standing_volume, LiDAR_dataset$standing_volume, 
                               paired = FALSE, alternative = "two.sided")

comparison_RGB_LiDAR_DBH <- t.test(photogrammetry_dataset$DBH, LiDAR_dataset$DBH, 
                               paired = FALSE, alternative = "two.sided")
library(Metrics)
rmse(photogrammetry_dataset$standing_volume, LiDAR_dataset$standing_volume)





