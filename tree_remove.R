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

#Set in and out directories   
inDir <- "D:/ptx/out/for_comparison/rotate/"
outDir <- "D:/ptx/out/for_comparison/no_trees2/"

#List input files
files <-list.files(pattern=".las", inDir)   

for(i in 14:length(files)){
  
  #Grab file name
  file <- files[i]
  #Read in file
  tls <- readTLS(paste0(inDir, file))
  
  #Prep LAS
  tls <- nnFilter2(tls)
  
  #Find trees and stems
  map <- treeMap(filter_poi(tls,Classification==0), 
                 map.hough(min_h = 2, max_h  = 3, h_step = 0.5, min_votes = 1,
                           pixel_size = 0.05, min_density = 0.01),merge=0)
  tls <- treePoints(tls, map, trp.crop())
  tls <- stemPoints(tls, stm.hough())
  tls@data[Stem == T, Classification := 20] 
  juststem  <- filter_poi(tls, Classification == 20L)
  justfuels <- filter_poi(tls, Classification == 0)
  
  #Filter out just fuels and save
  veght<- filter_poi(justfuels, Z > 0L, Z < 3L,TreeID=!0)
  writeLAS(veght, paste0(outDir, file))
}