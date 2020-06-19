# Forest Inventory Project Repository

## Contents:
1. AHN3 preprocessing
2. Creating a canopy height model from AHN3 and UAB LS/ UAV RGB
3. Segmentation of trees from UAV-LS/TLS data and computation of tree volume
4. Processing the RGB data and tree volume computation

## Input data for Speulderbos area
1. AHN3 puntenwolk (kaartblad 32fn2). Source: http://esrinlcontent.maps.arcgis.com/apps/Embed/index.html?appid=a3dfa5a818174aa787392e461c80f781
2. UAV LS point cloud dataset
3. TLS point clud dataset
4. UAV RGB orthophoto

## Instructions
1. ***AHN3_preprocessing.R***   
This script can be used for cropping/reprojecting and saving AHNs.
2. ***CHM_from RGB.R, CHM_from_AHN3.R, CHM_from_UAV_LS.R***   
The scripts than can be used to compute canopy height model (CHM) using three types of input data - UAV RGB imagegery, AHN3 point cloud   or UAV LS point cloud. Additionally tree metrics based on CHM are computed (crown areas, crown diameters and height of trees).
3. ***RGB_volume_estimation.R***  
The script estimates the DBH and volume of each tree based on CHM derivates and empirical models.
4. ***UAV_LS_volume_estimation.R***  
The script estimates the DBH and volume of each tree using a RandomForest model taking CHM derivates and manually collected training data as input.
5. ***TLS_beech_segmentation.R,TLS_fir_segmentation.R***  
This script is made in order to create the CHM using TLS data. Derivates like treetops are computed. The last part of the script is tree segmentation using Dalponte approach.
