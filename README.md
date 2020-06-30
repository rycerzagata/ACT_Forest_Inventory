# Forest Inventory Project Repository

## Contents:
1. Preprocessing of the data
2. Creating a canopy height model from AHN3 and the UAV/TLS data
3. Segmentation of trees from UAV-LS/TLS data
4. Prediction of tree volume from UAV-LS derived tree metrics using RandomForest model
5. Prediction of tree volume from UAV RGB using empirical equations
6. Validation of volume estimations from both UAV-LS and UAV RGB using TLS data

## Input data for Speulderbos area
Explanation of abbreviations: ps - project source, valid - validation dataset, ps01 - the project source for first stage of the project - data preparation, ps02 - the project source for the second stage of the project - volume estimation.
The data should be downloaded to a local device and stored in a folder *ACT_Forest_Inventory*, subfolder *Data*.The paths in scripts are adjusted to file names presented below. The scripts should automatically recognize the paths after setting up the correct working directory. Link to the data: https://filesender.surf.nl/?s=download&token=e30fbda8-bd14-4df9-990b-7d63c66e99c8
1. ***ps01_TLS_beech.laz, ps01_TLS_fir.laz***   
TLS point cloud datasets, note that for fir the dataset might be divided into a couple of files for the purpose of computing limitations so the scripts should be run for each file.  
2. ***ps02_AHN3.laz***  
AHN3 puntenwolk (kaartblad 32fn2). The downloaded dataset is cropped to the extent of the study area (175981, 473657, 176265, 473782) to reduce its size (from 4 GB to 4,6 MB). Source: http://esrinlcontent.maps.arcgis.com/apps/Embed/index.html?appid=a3dfa5a818174aa787392e461c80f781 
3. ***ps02_UAV_LS.laz***  
UAV LS point cloud dataset clipped to the Speulderbos study site (175981, 473657, 176265, 473782).  
4. ***ps02_UAV_RGB.tif***  
Digital Surface Model created from photogrammetry data cropped to the study area (175981, 473657, 176265, 473782).  

## Instructions  
1.  ***TLS_beech_segmentation.R and TLS_fir_segmentation.R***  
This script is made in order to create the CHM using TLS data. Derivates like treetops are computed. The last part of the script is tree segmentation using Dalponte's approach. TLS script should be run before the others, since the output of this script is needed to perform the validation for both RGB and UAV-LS scripts.  

2. ***RGB_beech.R and RGB_fir.R***   
The scripts can be used to estimate the tree volume using CHM based on orthophoto and AHN3 and empirical equations for beech/ Douglas fir forest.

3. ***UAV_LS_beech.R and UAV_LS_fir.R***   
The scripts can be used to compute tree volume using CHM derivates produced from UAV point cloud and AHN3 as well as the products of Dalponte's tree segmentation algorithm. The volume is predicted by RandomForest model with manually measured tree metrics used as training data.


