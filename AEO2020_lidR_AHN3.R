############################################################################
#### script AEO_LiDAR_AHN3
#### V1.0 2018-03-06
#### Harm Bartholomeus
#### 
#### AEO exercise on AHN3 lidar processing 
#### 
#### Instructions: Go through the script step by step and execute the R-code
#### (select code and Ctrl-R). Answer the questions in a seperate document.
#### If required you can also copy the figures to this answer document.
############################################################################


# For this part of the exercise you will work with pointcloud data which are collected for the 
# creation of the Actueel Hoogtemodel Nederland (up-to-date elevation model of the Netherlands), 
# version 3 (AHN3). Those data are collected every 6-7 years by means of Airborne Laser Scanning 
# to create a Digital Terrain Model with a spatial resolution of 5m. Next to the rasterized DTM 
# the pointcloud data are provided. 
# # 
# # The pointclouds are delivered in .las/.laz format where all individual points have been 
# classified into the classes: ground, building, water, artificial object or other 
# (ref: https://www.pdok.nl/nl/ahn3-downloads). Next to that, some additional attributes are 
# recorded for each point. The .las format is a binary format which is used to store and exchange 
# LiDAR data. A .laz file is a zipped/compressed version, which reduces the file size to ~10% 
# without loss of information.  
# # 
# # Since a few years the data are openly available through https://www.pdok.nl/nl/producten/pdok-downloads/atomfeeds 
# (AHN1 and AHN2) and https://www.pdok.nl/nl/ahn3-downloads (AHN3) . AHN1 and AHN2 data cover the whole of 
# the Netherlands, AHN3 has now been released for part of the Netherlands (status early 2018)
# # 
# # For the exercise the Tile: 39EZ2, covering the area around Rhenen was downloaded, which has 
# a file size of 2.5 GB in .laz format. Since this is a bit too much to process on the computers
# used during the course tiles of 200*200 m were created using the LAStile function in Lastools
# (https://rapidlasso.com/lastools/). This software is also available for use during the practicals! 
# See download section on blackboard, including the software license. Much of the functionality we show 
# in the rest of the practical can be done in LAStools as well, where LAStools is better capable to deal 
# with very large datasets. 
#
# For information on the functions check the lidR documentation on https://cran.r-project.org/web/packages/lidR/lidR.pdf
#
# And... when scripting... Google is your friend!

#### INSTALL PACKAGES ####
# in the PC rooms this can't be done with the script, so please do it manually:
# Go to: Tools - Install Packages
# As a repositry select: CRAN/CRANextra
# As a library path select: "C:/Program Files/R/R-3.4.3/library"
# you can install multiple packages at once, so just copy the names given below: 
# raster, lidR, colorRamps, sp, rgl
### pressing Install will start the installation which may take about 5 mins (coffee!)

## install.packages("lidR", "raster", "colorRamps", "sp", "rgl") ## normally running this works as well, but not in the pc rooms

#### PARAMETER SETTING #### 
library(lidR)
library(raster)
library(colorRamps)
library(sp)
library(rgl)

## set working directory. ADAPT IT TO YOUR OWN DATA DIRECTORY!!
setwd("/Users/marariza/Dropbox/MSc/Advanced Earth Observation/Practical_7/AEO_LIDAR_DATA_2020")
##

## Determine which file to load: 

LASfile <-"AHN3/tile_rhenen_west_small.laz" ## set the name of the file to be loaded

las= readLAS (LASfile) ## open the .las file, which will be named "las"
epsg(las) = 28992  ## set RD coordinate system


#### CREATE A SUBSET IF YOUR COMPUTER NEEDS IT:
## remove the ## at the beginning of the line below and run the line. This will give you 
## a small subset with most features in it as well

##las = lasclipRectangle(las, 166820, 441763, 166950, 441863)


#-- on the right panel you can see which data/variables are loaded in the memory of R and which values they have. 

#-- let's first check what the spatial extent of the file is: 

extent(las)

## Q1: What is the actual spatial extent of this tile?

#--> let's visualize the pointcloud - by default it will be colored using the Z value. 

plot(las) ## a 3D viewer will be opened with the point cloud displayed
rgl.close() ## to close the viewer again. Do so regularly to prevent memory issues. 

#### INVESTIGATING THE POINT CLOUD ####

#--> click on the blue circle with the white triangle next in front of the las file on the right. This will
# show all attributes that are linked to this dataset. 

## Q2: Which attributes are available for every single point?

#--> Let's display the points with some of the other attributes. You may have to close some viewers if your computer get's slow
plot(las, color="Intensity", colorPalette= heat.colors(50), trim=300)
rgl.close() 
## Q3: we don't know exactly which type of scanner and which  wavelength was used to acquire the AHN3 data, but can you explain the intensity pattern?

plot(las, color="ReturnNumber")
## Q4: Investigate, describe and explain the ReturnNumber distribtution. Look in the header of the file (where you checked the attributes as well). What is the distribution of 1st-2nd.... nth retutns?
rgl.close() 
#--> feel free to have a look at some of the other attributes as well. 

#--> Let's look at the classification

plot(las, color="Classification")

#--> these are the values present in the as classification

unique(las@data$Classification)

### There are fixed classification codes for las files. Use google to find which classes are present in this file

## Q5: how many and which classes are present? Give a table with the numbers and their meaning

rgl.close()  ## just close the 3d plots to keep the computers going

#--> now let's select the points which are classified as ground points only, there is a special function for this:
ground <- lasfilterground(las)

plot(ground) ## to plot all points classified as ground points

rground <- (grid_metrics(ground, res=1, mean(Z)))## to convert the ground points to a raster format 

plot(rground, col=matlab.like2(50))

#--> If you would have a raw point file you could create such a classification using LasClassify in Lastools. 

## if you want to select points defined by other attributes you can query your data: 

buildings <- las %>% lasfilter(Classification==6)

plot(buildings) ## visualize the points classified as buildings

## Q6: Do all buildings have points on all the side walls? Explain your observation. If you are doubting, run the line below. 

plot(buildings, colour="ScanAngle") ## visualize the points classified as buildings

## if you want to filter the data based on other attributes you can use this function
other <- las %>% lasfilter(Classification==1) ## this contains all vegetation and other objects (lamp posts etc)
plot(other, color="ReturnNumber")

#### CALCULTING LIDAR METRICS ####
## you can calculate values of the different attributes per grid cell. The lidR package uses grid_metrics for that, which returns the results as a raster object. 

## let's first determine some basic rasters
minHeight = grid_metrics(las, res=1, min(Z)) ## returns a 1m raster with the minimum Z value per grid cell
maxHeight = grid_metrics(las, res=1, max(Z)) ## returns a 1m raster with the maximum Z value per grid cell

plot(maxHeight, main="Maximum Height (m)", col=matlab.like2(50))
plot(minHeight, main="Minimum Height (m)", col=matlab.like2(50))

## from these two rasters you can calculate the difference in minimum and maximum height
rHeightdiff <- maxHeight - minHeight ## calculate the heightdifference
plot(rHeightdiff, main="Object Height", col=matlab.like2(50))

## Q7: what do you notice for e.g. the buildings when you look at the object height? Can you explain this?

##  Let's evaluate the point density distribution
density <- grid_density(las, res=1) ## calculate the number of points per surface area
plot(density, main="Point Density", col= matlab.like2(50)) ## to plot a map with the point density
hist(density, breaks=50, main="Histogram of point density per m2", xlab="Point Density per m2") ## create a frequency histogram 

## Q8: How many points do you have more or less per m2? 
## Q9: What does the variation in point density depend on?

### for creating a DEM you mainly need a lot of points on the ground. Let's check that density
ground_dens <- grid_density(ground, res=1) ## calculate density of gorund points
plot(ground_dens, main="Ground Points Density", col= matlab.like2(50))
hist(ground_dens, breaks=50, main="Histogram of Ground point density per m2", xlab="Point Density per m2") ## create a frequency histogram 

## Q10: How many ground points do you have more or less per m2? 

## let's have a look at the relative point density: 
rel_dens <- ground_dens/density
plot(rel_dens, zlim=c(0,1), col= matlab.like2(50), main="Relative Point Density")

## Q11: Why does the relative density show spatial patterns? What causes this spatial variation in Ground point density?

#### DTM CREATION  ####
## To create a DTM we can use the grid_terrain function in LidR. 
## To see how it works check the lidR documentation on https://cran.r-project.org/web/packages/lidR/lidR.pdf

## Q12: which interpolation methods area to create a DTM are available?

DTM <- grid_terrain(las, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest = FALSE)  ## create a DTM using IDW interpolation
plot(DTM, main="DTM", col=matlab.like2(50))

## Feel free to try out different interpolation techniques!!

### for the terrain interpolation you have to option to show interpolated values, or keep the lowest point per grid cell. 

## Let's check what the difference in output is if you use this function
DTMkl <- grid_terrain(las, res=1, algorithm = knnidw(k=6L, p = 2), keep_lowest=TRUE) ## interpolate DTM keeping lowest points 
plot(DTMkl, main="DTM - keep lowest point", col=matlab.like2(50))
kl_diff <- DTM - DTMkl ## calculate the difference between the two rasters
plot(kl_diff, col=matlab.like2(50), main="DTM Difference using different interpolation options", zlim=c(0,0.5)) ## create a plot that shows the difference

## Q13: Explain the spatial pattern you notice 

#### CHM CREATIONl ####
## Off course we can also look at the highest point, which is often done to create a canopy height model 
CHM<-grid_canopy(las, res=1, p2r(0.2))
plot(CHM, main="CHM", col=matlab.like2(50))
CHM2<-grid_canopy(las, res=1, dsmtin(0))
plot(CHM, main="CHM", col=matlab.like2(50))
dif = CHM - CHM2
plot(dif, col=matlab.like2(50), main="CHM Difference using different interpolation options", zlim=c(0,0.5)) ## create a plot that s

## Q14: does the CHM give the absolute height, or the height relative to the terrain?

## Q15: There are different algorithms available to create the Canopy Height model. Which one is used here? 
## Which other algorithms are available? Try at least 2 algorithms and calculate the difference. 
## Describe which methods you used, show the difference map and describe Where do you notice the major differences. 

#### Simple Terrain correction ####
## To calculate the height relative to the terrain you can substract the DTM from the CHM. 
Object_height <- CHM - DTMkl

plot(Object_height, main="Object Height", col=matlab.like2(50))
click(Object_height) ## now you can click in the viewer to read the raster values
## Q16: how high is the highest building? And how high are the trees? 
## press ESC to stop the click funtion

#### NORMALIZE DATA FOR TERRAIN HEIGHT####
### What you did above is just one way to calculate the object height. However, for some calculations it is better
### to correct the height of all points for the terrain height. This can be done by normalizing the z-values. 

nlas <- lasnormalize(las,DTM) ## here you correct all z-values for the corresponding terrain height.
## be aware that the DTM is used as a reference. You know there are different choices when creating that! 
## this will influence your result. Also the DTM resolution will influence the outcomes, especially in steep terrain. 
plot(nlas)  ## you can see that the terrain has become as flat as the Netherlands are often thought to be

hist(nlas@data$Z, main="Histogram of normalized Z values", xlab="Height (m)", ylab="Frequency") ## plot a histogram of the normalized Z values

## Using a query we can create a slices at a certain height above ground level

DBH_slice <-  nlas %>% lasfilter(Z>1 & Z<1.5) ## slice around Breast height 

Crowns_slice <-  nlas %>% lasfilter(Z>10) ## or for the tree crowns and higher parts of buildings

plot(DBH_slice, color="Classification")
plot(Crowns_slice, color="Classification")

#### CANOPY PROFILES ####
#### For vegetation analysis often vertical profiles are created. 
#### First, let's select a few vegetated areas with 10m diameter. 

sub1 <- lasclipCircle(nlas,166862,441800, 10 )
plot(sub1)
sub2<- lasclipCircle(nlas,166730,441851, 10 )
plot(sub2)
sub3<- lasclipCircle(nlas,166909,441836, 10 )
plot(sub3)
sub4<- lasclipCircle(nlas,166900,441790, 10 )
plot(sub4)

## For each subset we can calculate the gap fraction at different heights above the terrain
gapFraction1 <- (gap_fraction_profile(sub1@data$Z))
gapFraction2 <- (gap_fraction_profile(sub2@data$Z))
gapFraction3 <- (gap_fraction_profile(sub3@data$Z))
gapFraction4 <- (gap_fraction_profile(sub4@data$Z))

## And put this in a nice plot
plot(gapFraction1$z ~ gapFraction1$gf, type="l", xlim=c(0.8,1), ylim=c(0,20), main="Gap Fraction Profiles", xlab="Gap Fraction", ylab="Height [m]")
points(gapFraction2$z ~ gapFraction2$gf, type="l",col=2, lty=2)
points(gapFraction3$z ~ gapFraction3$gf, type="l",col=3, lty=3)
points(gapFraction4$z ~ gapFraction4$gf, type="l",col=4, lty=4)
legend(0.8,20, legend=c("1", "2", "3","4" ), col=c(1:4), lty=1:4, bty="n")

## Q17: describe the differences between the profiles. Which plot has the most dense understorey? Which plot is the most "open"?

#### INDIVIDUAL TREE SEGMENTATION ####
## Sometimes we want to extract individual trees from the individual point cloud
###For this different methods are available, please see the "lastrees" function in the documentation

Vegpoints_norm <- nlas %>% lasfilter(Classification==1) ## this contains all vegetation and other objects (lamp posts etc)

# We use the method of Li 2012 to extract individual trees.
# This may take some time, so grab a coffee or tea. 
trees <- lastrees(Vegpoints_norm, li2012(R=5, speed_up=10, hmin=5))  ## this function creates a new pointcloud where evry point has an extra attribute, being the treeID

plot(trees, color="treeID") ## visualise the tree segmentation

## Q18: investigate the tree segmentation visually and describe the quality. 

## Q19: what do the R, speed_up and hmin settings determine?

### to check how many trees are detected we can see what the max treeID is

(max(trees@data$treeID, na.rm=TRUE))  ## give the maximum treeID

## Q20: how many trees are detected?

## select a random tree for visulation: 
tree <- trees %>% lasfilter(treeID ==10)
plot(tree)

### Q21: how did the segmentation work out? What are good/bad points you notice. 

## let's make a nice figure for all trees and extract some parameters of each tree. 
## These are: height, x, y locations and the convex hull of the crown
## This does take quite some time! The results are best if you run it for all trees, but feel free to
## make graphs and calculate the parameters for only part of the trees. 
## In that case replace line 1 of the for-loop by: 
## for (i in 1:50){
## or any other number instead of 50... 

dir.create( "extracted_trees_Vegpoints") ## create a new directory to write the outputs to 

d=NULL ## create an empty file to write some statistical data per tree
for (i in 1:max(trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- trees %>% lasfilter(treeID==i, Classification==1)
  height <- max(tree@data$Z)
  meanx <- mean(tree@data$X)
  meanY <- mean(tree@data$Y)
  int <- mean(tree@data$Intensity)
  count <- as.numeric(length(tree@data$Z)) ## count the number of points for each tree
  ## calculate the crown area by fitting a convex hull around it
  canopy <- tree %>% lasfilter(Z>5)
  X <- as.vector(canopy@data$X)
  Y <- as.vector(canopy@data$Y)
  XY<- cbind(X,Y)
  hpts<-chull(XY)
  hpts<- c(hpts, hpts[1])
  crown.coords <- XY[hpts,]
  crown.poly <- Polygon(crown.coords, hole=F)
  crown.area <- crown.poly@area
  
  ## calculate PCA for XY values. The ratio between the length of the main axes is used as a controlvalue to check if it is really a tree
  PCA <- prcomp(XY, scale=F)
  PCA.ev.ratio <- PCA$sdev[1]/PCA$sdev[2]
  score <- ifelse(PCA.ev.ratio>1.2,0,1) ## so if the ratio is >1.2, it is probably not a proper tree. 
  
  # plot(XY)
  ### write the data to a table
  d = rbind(d, data.frame(i,meanx,meanY, height, crown.area,PCA.ev.ratio,score, count))

## make a nice figure for every individual tree  
png(paste("extracted_trees_Vegpoints/tree", i, ".png"), height=800, width=400) 
par(mfrow=c(2,2))
plot(tree@data$Y~tree@data$X, asp=1, pch=19, col= tree@data$NumberOfReturns, main="Topview - NumberOfReturns", xlab="X-coordinate", ylab="Y-coordinate")
plot(tree@data$Z~tree@data$X, asp=1, pch=19, col=tree@data$ReturnNumber, main ="Sideview - Returnnumber", xlab="X-coordinate", ylab="Height (m)")
plot(tree@data$Z~tree@data$Y, asp=1, pch=19, col=tree@data$Classification, main="Sideview - Classification", xlab="Y-coordinate", ylab="Height (m)")
hist(tree@data$Z, main="Height Profile", xlab="Elevation (m)")
dev.off()
}

### See if there is any relation between the height of the tree and the crown size
plot(d$height~d$crown.area, xlab="crown area (m2)", ylab="tree height (m)")  ### 

## Q22: is there a relation between tree height and crown dimension?

#### "URBAN" TREEs ####
#### Run the same analysis for another area, where trees grow more seperated from each other. 
## this part of the script just does the analysis. Reuse parts of the script above if you want to visualize it. 

LaSfile_urban_trees <- "AHN3/tile_urban_167800_441400.laz"

ulas= readLAS (LaSfile_urban_trees)  ## let's first load a copy of the data
extent(ulas)

#### CREATE A SUBSET IF YOUR COMPUTER NEEDS IT:
## remove the ## at the beginning of the line below and run the line. This will give you 
## a small subset with most features in it as well
## subset_ulas =lasclipRectangle(ulas, 167800, 441400, 167900, 441500)


DTMu <- grid_terrain(ulas, res=1, algorithm = knnidw(k=2L, p = 2), keep_lowest = FALSE) ## create a DTM using IDW interpolation, ignore the error message
ulas<- lasnormalize(ulas,DTMu) ## normalize the Z values
Vegpoints_norm_u <- ulas %>% lasfilter(Classification==1) ## this contains all vegetation and other objects (lamp posts etc)

# We use the method of Li 2012 to extract individual trees
####lastrees(Vegpoints_norm_u, "li2012", R=7.5, hmin=5)  ## this function adds an extra attribute to each point, being the treeID
urban_trees <- lastrees(Vegpoints_norm_u, li2012(R=10, speed_up=5, hmin=5))  ## this function adds an extra attribute to each point, being the treeID

plot(urban_trees, color="treeID") ## visualise the tree segmentation

## Q23: investigate the tree segmentation visually and describe the quality. 

### to check how many trees are detected we can see what the max treeID is

(max(urban_trees@data$treeID, na.rm=TRUE))  ## give the maximum treeID

## Q24: how many trees are detected?

## let's make a nice figure for each tree and extract some parameters of the tree
## These are: height, x, y locations and the convex hull of the crown

dir.create( "extracted_trees_Vegpoints_urban") ## create a new directory to write the outputs to 

d_u=NULL ## create an empty file to write some statistical data per tree

## feel free to decrease the number of trees which you use in you analysis again!

for (i in 1:max(urban_trees@data$treeID, na.rm=TRUE)){
  print(i)
  tree <- urban_trees %>% lasfilter(treeID==i, Classification==1)
  height <- max(tree@data$Z)
  meanx <- mean(tree@data$X)
  meanY <- mean(tree@data$Y)
  int <- mean(tree@data$Intensity)
  count <- as.numeric(length(tree@data$Z)) ## count the number of points for each tree
  ## calculate the crown area by fitting a convex hull around it
  canopy <- tree %>% lasfilter(Z>5)
  X <- as.vector(canopy@data$X)
  Y <- as.vector(canopy@data$Y)
  XY<- cbind(X,Y)
  hpts<-chull(XY)
  hpts<- c(hpts, hpts[1])
  crown.coords <- XY[hpts,]
  crown.poly <- Polygon(crown.coords, hole=F)
  crown.area <- crown.poly@area
  
  ## calculate PCA for XY values. The ratio between the length of the main axes is used as a controlvalue to check if it is really a tree
  PCA <- prcomp(XY, scale=F)
  PCA.ev.ratio <- PCA$sdev[1]/PCA$sdev[2]
  score <- ifelse(PCA.ev.ratio>1.2,0,1) ## so if the ratio is >2, it is not a round tree when you look from above. 
  
  # plot(XY)
  ### write the data to a table
  d_u = rbind(d_u, data.frame(i,meanx,meanY, height, crown.area,PCA.ev.ratio, score, count))
  
  
  ## make a nice figure for every individual tree  
  png(paste("extracted_trees_Vegpoints_urban/tree", i, ".png"), height=800, width=400) 
  par(mfrow=c(2,2))
  plot(tree@data$Y~tree@data$X, asp=1, pch=19, col= tree@data$NumberOfReturns, main="Topview - NumberOfReturns", xlab="X-coordinate", ylab="Y-coordinate")
  plot(tree@data$Z~tree@data$X, asp=1, pch=19, col=tree@data$ReturnNumber, main ="Sideview - Returnnumber", xlab="X-coordinate", ylab="Height (m)")
  plot(tree@data$Z~tree@data$Y, asp=1, pch=19, col=tree@data$Classification, main="Sideview - Classification", xlab="Y-coordinate", ylab="Height (m)")
  hist(tree@data$Z, main="Height Profile", xlab="Elevation (m)")
  dev.off()
}

### let's remove all trees that don't have an elongated crown according to the PCA-axes and contain less than 500 points to clean up the dataset a little:
d.sel <- d[which(d$PCA.ev.ratio<1.3 & d$count>500),]
d_u.sel <- d_u[which(d_u$PCA.ev.ratio<1.3 & d_u$count>100),]

## Let's plot the selected trees of both areas in one graph:
plot(d.sel$height~d.sel$crown.area,ylim=c(0,25), xlim=c(0,400), xlab="crown area (m2)", ylab="tree height (m)", main="Relation tree height vs crown area - large round trees")  ### 
points(d_u.sel$height~d_u.sel$crown.area, col=2, pch=2)  ### 
legend(300,4, legend=c("forest trees", "urban trees" ), col=c(1:2),pch=c(1:2), bty="n")

## Q25: is there a relation between tree height and crown dimension in the more urban area and is it different from the forest dataset?
## Is this an ecologically logical relation? Which improvements would you need in terms of data collection and processing to prove this?

#### THE END ####
### please hand in your answers in blackboard, together with the exercise on TLS ###

