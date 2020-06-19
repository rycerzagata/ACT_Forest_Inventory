# Forest Inventory Project Repository

## Contents:
1. AHN3 preprocessing
2. Creating a canopy height model from AHN3 and UAB LS/UAV RGB
3. Segmentation of trees from UAV-LS/TLS data and computation of tree volume
4. Processing the RGB data and tree volume estimation
5. Validation of volume estimations using TLS data

## Input data for Speulderbos area
1. AHN3 puntenwolk (kaartblad 32fn2). Source: http://esrinlcontent.maps.arcgis.com/apps/Embed/index.html?appid=a3dfa5a818174aa787392e461c80f781
2. UAV LS point cloud dataset
3. UAV RGB orthophoto
4. TLS point clud dataset

## Instructions
1. ***RGB_beech.R and RGB_fir.R***   
The scripts can be used to estimate the tree volume using CHM based on orthophoto and AHN3 and empirical equations for beech/ Douglas fir forest.
2. ***UAV_LS_beech.R and UAV_LS_fir.R***   
The scripts than can be used to compute tree volume using CHM derivates produced from UAV point cloud and AHN3 as well as the products of Dalponte's tree segmentation algorithm. The volume is predicted by RandomForest model with manually measured tree metrics used as training data.
3. ***TLS_beech_segmentation.R and TLS_fir_segmentation.R***  
This script is made in order to create the CHM using TLS data. Derivates like treetops are computed. The last part of the script is tree segmentation using Dalponte approach.
