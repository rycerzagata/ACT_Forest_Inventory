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

## Setting working directory
setwd("/Users/marariza/Downloads")

# Loading the downloaded datasets (32fn2 from AHN3)
AHN3_laz <- "Data/C_32FN2.laz"
AHN3 <- catalog(AHN3_laz)

# Checking the current extent
extent(AHN3)

# Clipping the extent to the study area (https://www.gpscoordinaten.nl/converteer-rd-coordinaten.php)
# The website above can be used to convert global coordinates to RD New
AHN3_clip <- lasclipRectangle(AHN3, 175981, 473657, 176265, 473782)

# Checking the new extent
extent(AHN3_clip)

# We reproject it to RD New (it should already RD new, just to make sure)
epsg(AHN3_clip) = 28992

# Plot it to see the result
plot(AHN3_clip)
rgl.close()

# Save the clipped AHN3 for future use
writeLAS(AHN3_clip, "AHN3.laz")
writeLAS(AHN3_clip, "AHN3.las")

# Select area with beech forest and save it in separate files
ahn3laz <- "/Users/HP/Documents/ACT/R/Data/AHN3.laz"
AHN3 <- readLAS(ahn3laz)
AHN3_beech <- lasclipRectangle(AHN3, 176170, 473657, 176265, 473782)
#plot(AHN_beech)
writeLAS(AHN3_beech, "/Users/HP/Documents/ACT/R/Data/AHN3_beech.laz")
writeLAS(AHN3_beech, "/Users/HP/Documents/ACT/R/Data/AHN3_beech.las")

# Select area with Douglas Fir forest and save it in separate files
AHN3_fir <- lasclipRectangle(AHN3, 176040, 473657, 176100, 473710)
#plot(AHN3_fir)
writeLAS(AHN3_fir, "/Users/HP/Documents/ACT/R/Data/AHN3_fir.laz")
writeLAS(AHN3_fir, "/Users/HP/Documents/ACT/R/Data/AHN3_fir.las")


######################################################################################################



