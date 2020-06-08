setwd("E:/ACT DATA/2017_SpeulderbosLS")
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

#load data
UAVfile <-"UAV_withGround.laz"
UAV = readLAS (UAVfile)
epsg(UAV) = 28992

#DTM creation 
DTM <- grid_terrain(UAV, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)  ## create a DTM using IDW interpolation

#creat DSM
DSM<-grid_canopy(UAV, res=1, p2r(0.2))

#OBJECT Height
CHM <- DSM - DTM

#Remove negetative data
CHM[CHM<0] <- NA

# We use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After we plot it to check how the tree tops look like. 
lin <- function(x) {x*0.05+0.6}
treetops <- vwf(CHM = CHM, winFun = lin, minHeight = 15)

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


##DO same thing with beech subsit
#clip out beech region
beech <- lasclipRectangle(UAV,176170, 473657, 176265, 473782)

#DTM creation 
DTM_beech <- grid_terrain(beech, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)  ## create a DTM using IDW interpolation

#creat DSM
DSM_beech <- grid_canopy(beech, res=1, p2r(0.2))

#OBJECT Height
CHM_beech <- DSM_beech - DTM_beech

#Remove negetative data
CHM_beech[CHM<0] <- NA

# We use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After we plot it to check how the tree tops look like. 
lin_beech <- function(x) {x*0.05+0.6}
treetops_beech <- vwf(CHM = CHM_beech, winFun = lin, minHeight = 15)

plot(CHM_beech, main="CHM_beech", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(treetops_beech, col="black", pch = 20, cex=0.5, add=TRUE)

# We compute the function MCWS function that implements the watershed algorithm. In this case, the argument
# minHeight refers to the lowest expected treetop. The result is a raster where each tree crown is 
# a unique cell value. 
crowns_beech <- mcws(treetops = treetops_beech, CHM=CHM_beech, minHeight = 15, verbose=FALSE)
plot(crowns, main="Detected Beech crowns", col=sample(rainbow(50), length(unique(crowns[])),replace=TRUE), 
     legend=FALSE, xaxt="n", yaxt="n")

# We do the same computation as before but changig the output format to polygons. It takes more processing
# time but polygons inherit the attributes of treetops as height. Also, crown area is computed for each polygon.
crownsPoly_beech <- mcws(treetops = treetops_beech, CHM=CHM_beech, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM_beech, main="CHM_beech", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly_beech, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly_beech[["crownDiameter"]] <- sqrt(crownsPoly_beech$crownArea/pi) *2
mean(crownsPoly_beech$crownDiameter)
mean(crownsPoly_beech$crownArea)

sp_summarise(treetops_beech)
sp_summarise(crownsPoly_beech, variables=c("crownArea", "height"))

