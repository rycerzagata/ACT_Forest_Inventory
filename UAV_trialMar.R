
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
library(itcSegment)
library(TreeLS)

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

#Remove NA values
CHM[is.na(CHM)] <- 0

# Using focal statistics to smooth the CHM
CHM_smooth <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)
plot(CHM_smooth)

# We use the Variable Window Filter (VWF) to detect dominant tree tops. We use a linear function used in 
# forestry and set the minimum height of trees at 10, but those variables can be modified. 
# After we plot it to check how the tree tops look like. 
lin <- function(x) {x*0.07+3}
treetops <- vwf(CHM = CHM_smooth, winFun = lin, minHeight = 15)

plot(CHM_smooth, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(treetops, col="black", pch = 20, cex=0.5, add=TRUE)

# We compute the function MCWS function that implements the watershed algorithm. In this case, the argument
# minHeight refers to the lowest expected treetop. The result is a raster where each tree crown is 
# a unique cell value. 
crowns <- mcws(treetops = treetops, CHM=CHM_smooth, minHeight = 15, verbose=FALSE)
plot(crowns, main="Detected tree crowns", col=sample(rainbow(50), length(unique(crowns[])),replace=TRUE), 
     legend=FALSE, xaxt="n", yaxt="n")

# We do the same computation as before but changig the output format to polygons. It takes more processing
# time but polygons inherit the attributes of treetops as height. Also, crown area is computed for each polygon.
crownsPoly <- mcws(treetops = treetops, CHM=CHM_smooth, minHeight = 8, verbose=FALSE, format="polygons")
plot(CHM_smooth, main="CHM", col=matlab.like2(50), xaxt="n", yaxt="n")
plot(crownsPoly, border="black", lwd=0.5, add=TRUE)

# Assuming each crown has a roughly circular shape,the crown area is used to compute its average circular diameter.
crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly$crownArea/pi) *2
mean(crownsPoly$crownDiameter)
mean(crownsPoly$crownArea)

# Make a summary of computed tree parameters
sp_summarise(treetops)
sp_summarise(crownsPoly, variables=c("crownArea", "height"))

#####################################################################################################
library(randomForest)
set.seed(5)

# Create a dataframe out of the polygons
training <- as.data.frame(crownsPoly)

# Creating Train and Test samples
split <- runif(nrow(training),0,1)>0.3
train_samples <- subset(training, split)
test_samples <- subset(training, !split)


#train_ind <- sample(seq_len(nrow(training)), size = 50)
#train_Agata <- training[train_ind,]
#val_Agata <- training[-train_ind,]

# Checking the number of row of each splitted dataset
nrow(train_samples)
nrow(test_samples)

# Train a Random Forest model to predict DBH
train_samples$DBH <- rnorm(nrow(train_samples), mean=0.4, sd=0.1)
model <- randomForest(DBH ~ height + crownArea, 
                      data = train_samples, importance= TRUE, proximity = TRUE, ntree=200)

# Predict the DBH of the test dataset
test_samples$DBH <- predict(model, test_samples)
test_samples

# Merge both datasets
full_dataset <- rbind(train_samples, test_samples)

# Compute the standing volume with the DBH and the height. Source: https://silvafennica.fi/pdf/smf004.pdf
full_dataset$standing_volume <- ((0.049/100)*(full_dataset$DBH^1.78189)*(full_dataset$height)^1.08345)*1000

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
#plot(DBH_slice, color="Classification")

# Select all vegetation and other objects
vegpoints_norm <- nlas %>% lasfilter(Classification==1)

# Individual tree segmentation 
trees <- lastrees(vegpoints_norm, silva2016(CHM_smooth, treetops, max_cr_factor = 0.8, exclusion = 0.1, ID = "treeID"))
plot(trees, color="treeID")

# Check the number of trees we have
(max(trees@data$treeID, na.rm=TRUE))

# Normalize UAV point cloud dataset and perform stem segmentation with RANSAC 
#tls = tlsNormalize(UAV_beech, keepGround = FALSE)
#thin = tlsSample(tls, voxelize(0.01))
#map = treeMap(thin, map.hough(hmin = 1, hmax = 2, max_radius = 0.3, min_density = 0.01, min_votes = 2))
#tls2 = stemPoints(tls, map)
#df = stemSegmentation(tls2, sgmt.ransac.circle(n=10))
#tlsPlot(tls, df, map)

# We extract every tree into a different .laz file
dir.create( "extracted_laz")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("extracted_laz/tree", i, ".laz"))}

##################################################################################################



