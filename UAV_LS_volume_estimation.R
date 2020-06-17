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
#Setting working directory
setwd("/Users/marariza/Downloads")

set.seed(2020)
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


#### TREE SEGMENTATION - DALPONTE APPROACH ####
set.seed(2020)
# Treetops detection
f <- function(x) { x * 0.07 + 3 }
ttops <- tree_detection(CHM, lmf(f))  # or lmf(4, 2)
#ttops <- tree_detection(beechLas, lmf(5))
#x <- plot(beechLas)
#add_treetops3d(x, ttops)
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(ttops, col="black", pch = 20, cex=0.5, add=TRUE)

# Crowns
crownsPoly <- mcws(treetops = ttops, CHM=CHM, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Calculate diameter of crowns and add to crown data
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2

# Normalize the point cloud using DTM
nlas <- lasnormalize(beechLas, DTM)
Vegpoints_norm <- nlas %>% lasfilter(Classification==1) 
trees <- lastrees(nlas, dalponte2016(chm = CHM, treetops = ttops, ID = 'treeID'))
#plot(trees, color="treeID")

# We extract every tree into a different .laz file
dir.create( "/Users/HP/Documents/ACT/R/Data/extracted_trees")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("/Users/HP/Documents/ACT/R/Data/extracted_trees/tree", i, ".laz"))}



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





