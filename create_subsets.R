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

#Set input and output directories
inDir <- "D:/ptx/out/for_comparison/no_trees3/"
outDir <- "D:/ptx/out/for_comparison/all_subsets/"

r1 <- 4
r2 <- 6

files <-list.files(pattern=".las", inDir)   


for(i in 1:length(files)){
  
  file <- files[i]
  lasIn <- readLAS(paste0(inDir, file))
  clipIt(lasIn, paste0(outDir, "clip/", file), 0 , 10)
  lasIn <- readLAS(paste0(outDir, "clip/", file))
  clipIt(lasIn, paste0(outDir, "donut/", file), r1 , r2)
  sN <- clip_circle(lasIn, 0, 5.5, 1)
  sS <- clip_circle(lasIn, 0, -4.5, 1)
  sW <- clip_circle(lasIn, -5.5, 0, 1)
  sE <- clip_circle(lasIn, 4.5, 0, 1)
  writeLAS(sN, paste0(outDir, "siteN/", file))
  writeLAS(sS, paste0(outDir, "siteS/", file))
  writeLAS(sW, paste0(outDir, "siteW/", file))
  writeLAS(sE, paste0(outDir, "siteE/", file))
  
  tN <- clip_transect(lasIn, c(0,.5), c(0,10), 1)
  tS <- clip_transect(lasIn, c(0,-.5), c(0,-10), 1)
  tW <- clip_transect(lasIn, c(-.5,0), c(-10,0), 1)
  tE <- clip_transect(lasIn, c(.5,0), c(10,0), 1)
  
  writeLAS(tN, paste0(outDir, "transectsScratch/N_", file))
  writeLAS(tS, paste0(outDir, "transectsScratch/S_", file))
  writeLAS(tW, paste0(outDir, "transectsScratch/W_", file))
  writeLAS(tE, paste0(outDir, "transectsScratch/E_", file))
  filelist <- c(paste0(outDir, "transectsScratch/N_", file), 
                paste0(outDir, "transectsScratch/S_", file), 
                paste0(outDir, "transectsScratch/W_", file), 
                paste0(outDir, "transectsScratch/E_", file))
  allT <- rlas::read.las(filelist)
  write.las(paste0(outDir, "transects/", file), lhead, allT)
  print(paste0("Completed for ", file))
}