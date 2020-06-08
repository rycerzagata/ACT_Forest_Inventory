"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
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

## Setting working directory
setwd("/Users/marariza/Downloads")

# We load and read the AHN3 file
AHN3_clip <- "AHN3.laz"
AHN3 <- readLAS(AHN3_clip)

# We compute two DTMs modifying just one parameter (keep lowest)
DTM1 <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
DTM2 <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest=TRUE)

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
RGB <- "DEM_speulderbos_georef_ar.tif"
DSM_no_proj <- raster(RGB)
DSM_proj <- projectRaster(DSM_no_proj, crs = crs(DTM1))
DSM <- resample(DSM_proj, DTM1, method = "bilinear")

# We compute the CHM substracting also 40.68, which is the difference height between the DSM and the DTM
# in the area on the left (U-shaped area).  
CHM <- DSM - DTM1 - 40.68
plot(CHM, main="CHM", col=matlab.like2(50))

##########################################################################################################

plot(ground)
rground <- (grid_metrics(ground, res=1, mean(Z)))
plot(rground, col=matlab.like2(50))
minHeight = grid_metrics(AHN3_clip, res=1, min(Z)) 
maxHeight = grid_metrics(AHN3_clip, res=1, max(Z))
plot(maxHeight, main="Maximum Height (m)", col=matlab.like2(50))
plot(minHeight, main="Minimum Height (m)", col=matlab.like2(50))
rHeightdiff <- maxHeight - minHeight 
plot(rHeightdiff, main="Object Height", col=matlab.like2(50))
DTM <- grid_terrain(AHN3_clip, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE) 
plot(DTM, main="DTM", col=matlab.like2(50))
CHM<-grid_canopy(AHN3_clip, res=1, p2r(0.2))
plot(CHM, main="CHM", col=matlab.like2(50))

