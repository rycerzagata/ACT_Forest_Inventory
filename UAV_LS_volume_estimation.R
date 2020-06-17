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
AHN3_clip <- "/Users/HP/Documents/ACT/R/Data/AHN3_beech.laz"
AHN3 <- readLAS(AHN3_clip)

# Load and clip the Laz file
lasfile <- "/Users/HP/Documents/ACT/R/Data/UAV_withGround.laz"
beechLas <- readLAS(lasfile)
beechLas <- lasclipRectangle(beechLas, 176170, 473657, 176265, 473782)

# Compute the DSM with the AHN3 dataset
DSM <- grid_canopy(beechLas, res=1, p2r(0.2))
#plot(DSM, main="DSM", col=matlab.like2(50))

# Compute the DTM with the AHN3 dataset
DTM <- grid_terrain(AHN3, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)
#plot(DTM, main="DTM", col=matlab.like2(50))

# Compute the CHM and remove one value which is below 0 (-0.005 m)
CHM <- DSM - DTM
CHM[is.na(CHM)] <- 0

# Use focal statistics to smoothen the CHM
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)


#### TREE SEGMENTATION - DALPONTE APPROACH ####
set.seed(2020)

# Treetops detection using an algorithm based on a local maximum filter.
f <- function(x) { x * 0.08 + 2 }
ttops <- tree_detection(CHM, lmf(f))  # or lmf(4, 2) for las
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
dir.create( "/Users/HP/Documents/ACT/R/Data/extracted_trees")
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  writeLAS(tree, paste("/Users/HP/Documents/ACT/R/Data/extracted_trees/tree", i, ".laz"))}


#### DBH PREDICTION ####
library(randomForest)
set.seed(2020)

# Create a dataframe out of the crown polygons
training <- as.data.frame(crownsPoly)
names(training) <- c("treeID", "height", "crownArea", "crownDiameter")

# Create Training, Test and Validation samples
train_ind <- sample(seq_len(nrow(training)), size = 70)
train_samples <- training[train_ind,]
test_samples <- training[-train_ind,]
valid_ind <- sample(seq_len(nrow(train_samples)), size = 20)
valid_samples <- train_samples[valid_ind,]
train_samples <- train_samples[-valid_ind,]

# Check the number of row of each splitted dataset
nrow(train_samples)
nrow(test_samples)
nrow(valid_samples)

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

total_volume <- sum(as.matrix(full_dataset$standing_volume))
total_area <- raster::area(beechLas)
m3ha <- totalVolume/(total_area/10000)
m3ha

write.csv(full_dataset,"/Users/marariza/Downloads/UAV-LS-results.csv", row.names = TRUE)

