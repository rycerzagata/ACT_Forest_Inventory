
##

## Determine which file to load: 

library(lidR)
library(raster)
library(colorRamps)
library(sp)
library(rgl)
library(ggpubr)
library(rlas)
library(tiff)
library(ForestTools)
## Setting working directory
setwd("/Users/marariza/Downloads")


#load data
UAVfile <-"UAV_withGround.laz"
UAV = readLAS (UAVfile)
epsg(UAV) = 28992
UAV_beech <- lasclipRectangle(UAV, 176170, 473657, 176265, 473782)

# We load and read the AHN3 file
#AHN3_clip <- "AHN3.laz"
#AHN3 <- readLAS(AHN3_clip)
# Computing the DTM with the AHN3 dataset
#DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
#plot(DTM, main="DTM", col=matlab.like2(50))

# We load and read the beech AHN3 file
AHN3beech <- "AHN3_beech.laz"
AHN3_beech <- readLAS(AHN3beech)
DTM <- grid_terrain(AHN3_beech, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
plot(DTM, main="DTM", col=matlab.like2(50))

#creat DSM
DSM<-grid_canopy(UAV_beech, res=1, p2r(0.2))

#OBJECT Height
CHM <- DSM - DTM

#Remove negetative data
CHM[CHM<0] <- NA

# Using focal statistics to smooth the CHM
CHM_smooth <- focal(CHM,w=matrix(1/9, nc=3, nr=3), fun=mean, na.rm=TRUE)
calc(CHM_smooth, function(x) {x*10})


par(mfrow=c(1,2))
plot(CHM, main="Regular")
plot(CHM_smooth, main="Smooth by 3")


# We use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After we plot it to check how the tree tops look like. 
lin <- function(x) {x*0.09+2}
treetops <- vwf(CHM = CHM_smooth, winFun = lin, minHeight = 15)

plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(treetops, col="black", pch = 20, cex=0.5, add=TRUE)

# We compute the function MCWS function that implements the watershed algorithm. In this case, the argument
# minHeight refers to the lowest expected treetop. The result is a raster where each tree crown is 
# a unique cell value. 
crowns <- mcws(treetops = treetops, CHM=CHM, minHeight = 15, verbose=FALSE)
plot(crowns, main="Detected tree crowns", col=sample(rainbow(50), length(unique(crowns[])),replace=TRUE), 
     legend=FALSE, xaxt="n", yaxt="n")

# We do the same computation as before but changig the output format to polygons. It takes more processing
# time but polygons inherit the attributes of treetops as height. Also, crown area is computed for each polygon.
crownsPoly <- mcws(treetops = treetops, CHM=CHM, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2
mean(crownsPoly$crownDiameter)
mean(crownsPoly$crownArea)

sp_summarise(treetops)
sp_summarise(crownsPoly, variables=c("crownArea", "height"))

#####################################################################################################

# Point density
density <- grid_density(UAV_beech, res=1)
plot(density)

# Normalize laz to correct height of all points for the terrain height
nlas <- lasnormalize(UAV_beech, DTM)
#plot(nlas)

# Select stems/crowns segments
DBH_slice <-  nlas %>% lasfilter(Z>0.8 & Z<1.5) ## slice around Breast height 
crown_slice <- nlas %>% lasfilter(Z>10)
plot(DBH_slice, color="Classification")

# Select all vegetation and other objects
vegpoints_norm <- nlas %>% lasfilter(Classification==1)
trees <- lastrees(vegpoints_norm, dalponte2016(CHM, treetops))
plot(trees, color="treeID")



stem_segm <- sgmt.ransac.circle(conf=0.99,inliers = 0.7,n=10, tol=0.025)

