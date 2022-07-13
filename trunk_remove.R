library(nabor)
library(dplyr)
library(tmap)
library(data.table)
library(terra)
library(rlas)
library(raster)
library(lidR)
library(e1071)
library(ggplot2)
library(plotly)
library(geometry)
library(MASS)
library(sampSurf)
library(rgdal)
library(sf)
library(Morpho)
library(TreeLS)

#==================================================================================

# Create a LAS header
lasdata = data.frame(X = c(0.001, 0.001, 0.001),
                     Y = c(0.001, 0.001, 0.001),
                     Z = c(0.001, 0.001, 0.001),
                     R = c(0.001, 0.001, 0.001),
                     G = c(0.001, 0.001, 0.001),
                     B = c(0.001, 0.001, 0.001),
                     I = c(0.001, 0.001, 0.001),
                     gpstime = c(0L, 0L, 0L),
                     Intensity = c(0L, 0L, 0L),
                     ReturnNumber = c(0L, 0L, 0L),
                     NumberOfReturns = c(0L, 0L, 0L),
                     ScanDirectionFlag = c(0L, 0L, 0L),
                     EdgeOfFlightline = c(0L, 0L, 0L),
                     Classification = c(0L, 0L, 0L),
                     ScanAngleRank = c(0L, 0L, 0L),
                     UserData = c(0L, 0L, 0L),
                     PointSourceID = c(0L, 0L, 0L))

lhead = header_create(lasdata)


inDir <- "D:/ptx/out/for_comparison/rotate/"
outDir <- "D:/ptx/out/for_comparison/no_trees3/"

files <-list.files(pattern=".las", inDir)   

for(i in 148:length(files)){
  
  file <- files[i]
  tlsIn <- readTLS(paste0(inDir, file))
  
  #Prep LAS
  tls2 <- nnFilter2(tlsIn)
  las_norm = tlsNormalize(tls2, keep_ground = T)
  las_norm <- filter_poi(las_norm, Z >= 0L)
  thin = tlsSample(las_norm, smp.voxelize(0.02))
  
  tls <- las_norm
  
  # Tree Segmentation  
  all_metrics = fastPointMetrics.available()
  my_metrics = all_metrics[c(16, 11)]
  Stems = fastPointMetrics(thin, ptm.knn(50), my_metrics)
  
  tlsfilter <- filter_poi(Stems, Verticality > 80, Verticality < 95)
  tlsfilter <- filter_poi(tlsfilter, Eigentropy < .03)
  map = treeMap(tlsfilter, map.hough(min_h = 2, max_h= 4, min_votes = 1), merge = 0)
  tls = treePoints(tls, map, trp.crop())
  tls = stemPoints(tls, stm.hough(h_step = 0.2, h_base = c(0.05, 0.25), min_votes = 1))
  tls@data[Stem == T, Classification := 20] 
  justfuels <- filter_poi(tls, Classification == 0 | Classification == 2)
  writeLAS(justfuels, paste0(outDir, file))
}

# Understory Segmentation of Terrestrial Laser Scanning. 
# Aim: Script to remove the trees from point cloud files in order to only account for fuels. 
# Authors: Guillén, L.A. using information adapted from Tall Timbers script. 
# 2020-12-04 
# Arguments: 

# load libraries -----
library(TreeLS)

# Declare functions ----

understory_segmentation_function<-function(list,source_path,new_path){ #open function
  #Aim: Function to obtain only understory vegetation and fuel bed LAS files
  #Arguments: list:files that are preprocessed, source_path: the path to the files, new_path:folder where to save new files. 
  #Output: writes a new file to the disk that is with understory vegetation.
  n<-length(list) #read the number of files
  i<-1 # initiate the value of counter
  while(i<=n){# open file/plot loop
    tls <- readTLS(paste0(source_path,list[i]))#read the ".las" file
    map <- treeMap(filter_poi(tls,Classification==0), # Map the point that form trees based on the stem geometry
                   map.hough(min_h = 2, max_h  = 3, h_step = 0.5, min_votes = 1,# and proximity.
                             pixel_size = 0.05, min_density = 0.01),merge=0)#  
    tls <- treePoints(tls, map, trp.crop())# obtain the point of the trees previously mapped
    tls <- stemPoints(tls, stm.hough())# obtain the points that form the stems
    seg <- stemSegmentation(tls, sgt.ransac.cylinder(n = 10))# This process created cylinders on the stems. #check usefulness. 
    tls@data[Stem == T, Classification := 20] # assign value "20" to the points that are stems. 
    juststem  <- filter_poi(tls, Classification == 20L)# create point cloud of only the tree stems. 
    justfuels <- filter_poi(tls, Classification == 0) # create point cloud of only fuels. 
    veght<- filter_poi(justfuels, Z > 0L, Z < 3L,TreeID=!0)# filter point cloud to only understory (not ground and lower than 3m)
    writeLAS(veght,file=paste0(new_path,list[i]))# write/create new file of only understory
    i<-i+1 # increase counter
  }#close file/plot loop
}#close function

# Applying functions ----

files<-list.files(path = "C:/LiDAR/USFWS/LAG_processing/voxelized/",pattern=".las") # read files. 

understory_segmentation_function(files,source_path ="C:/LiDAR/USFWS/LAG_processing/voxelized/", 
                                 new_path="C:/LiDAR/USFWS/LAG_processing/understory/") #apply function.



