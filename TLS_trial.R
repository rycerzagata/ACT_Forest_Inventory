"""
ACT group 10
Remote Sensing and GIS Integration 2020
Title: Forest Inventory through UAV based remote sensing
"""

# Loading the required libraries
library(raster)
library(sp)
library(rgl)
library (lidR)
library (TreeLS)

## Setting working directory

## read laz file of TLS data
TLS <- "Data/LS_x=105_y=105_TLS.laz"
TLS_a <- readLAS(TLS)
summary(TLS_a)
plot (TLS_a)

setTLS(TLS_a, colNames = NULL)

# normalize the point cloud
tls = tlsNormalize(TLS_a)
# map the trees on a resampled point cloud so all trees have approximately the same point density
thin = tlsSample(tls, voxelize(0.02))
map = treeMap(thin, map.hough(min_density = 0.05))
tls = stemPoints(tls, map)
df = stemSegmentation(tls, sgmt.ransac.circle(n=10))
head(df)
tlsPlot(tls, df, map)