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

#Read in and subset list of plots in which placard was found
theCoords <- read.csv("D:\\ptx\\TLS_Placard_Coords.csv")
theCoords2 <- theCoords[,1:4]
names(theCoords2) <- c("Name", "X", "Y", "Z")
theCoords3 <- theCoords2 %>% filter(complete.cases(.))

#Define ground classification algorithm
mycsf <- csf(TRUE, 1, 1, time_step = 1)

#Set directories
inDir <- "D:/ptx/in/"
outDir <- "D:/ptx/out/for_comparison/"

for(t in 1:nrow(theCoords3)){
  #Grab file name
  curRow <- theCoords3[t,]
  thePTX <- paste0(inDir, curRow$Name, ".ptx")
  
  #Convert PTX to LAS
  ptxToLas(thePTX, paste0(outDir, "las/", curRow$Name, ".las"), lhead)
  
  #Read in PTX
  tls <- readLAS(paste0(outDir, "las/", curRow$Name, ".las"))
  
  #Rotate to North align and clip
  tlsRot <- tlsRotate(tls=tls, x=curRow$X, y=curRow$Y, z=curRow$Z)
  tlsClp <- clip_circle(tlsRot, xcenter=0, ycenter=0, radius=10)
  
  #Classify and nromalize
  tlsCls <- classify_ground(tlsClp, mycsf)
  tlsNorm = tlsNormalize(tlsCls, min_res = 0.05, keep_ground = T)
  
  #Save result
  writeLAS(tlsNorm, paste0(outDir, "rotate/", curRow$Name, ".las"))
}

