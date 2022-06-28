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

#Read in placard coords table
theCoords <- read.csv("D:\\ptx\\TLS_Placard_Coords.csv")
theCoords2 <- theCoords[,1:4]
names(theCoords2) <- c("Name", "X", "Y", "Z")
theCoords3 <- theCoords2 %>% filter(complete.cases(.))

#Set input and output directories
inDir <- "D:/ptx_data/"
outDir <- "D:/ptx/out/"

circleRad <- 10
r1 <- 4.5
r2 <- 6.5
sampDist <- 5.5

ptxIn <- paste0(inDir, "2021_0927_1557.ptx")
ptxToLas(ptxIn, paste0(outDir, "las/", "ptx1.las"), lhead)
lasIn <- readLAS(paste0(outDir, "las/", "ptx1.las"))
lasRot <- tlsRotate(lasIn, -8, 6, 2)
writeLAS(lasRot, paste0(outDir, "rotate/", "ptx1.las"))
lasIn <- readLAS(paste0(outDir, "rotate/", "ptx1.las"))
clipIt(lasIn, paste0(outDir, "plot/", "ptx1.las"), 0 , 10)
lasIn <- readLAS(paste0(outDir, "plot/", "ptx1.las"))
mycsf <- csf(TRUE, 1, 1, time_step = 1)
lasCls <- classify_ground(lasIn, mycsf)
lasNorm <- normalize_height(lasCls,tin())
writeLAS(lasNorm, paste0(outDir, "normalized/", "ptx1.las"))
lasIn <- readLAS(paste0(outDir, "normalized/", "ptx1.las"))

clipIt(lasIn, paste0(outDir, "donut/", "ptx1.las"), r1 , r2)

sN <- clip_circle(lasIn, 0, 5.5, 1)
sS <- clip_circle(lasIn, 0, -5.5, 1)
sW <- clip_circle(lasIn, -5.5, 0, 1)
sE <- clip_circle(lasIn, 5.5, 0, 1)
writeLAS(sN, paste0(outDir, "sites/", "ptx1_N.las"))
writeLAS(sS, paste0(outDir, "sites/", "ptx1_S.las"))
writeLAS(sW, paste0(outDir, "sites/", "ptx1_W.las"))
writeLAS(sE, paste0(outDir, "sites/", "ptx1_E.las"))

tN <- clip_transect(lasIn, c(0,0), c(0,10), .5)
tS <- clip_transect(lasIn, c(0,0), c(0,-10), .5)
tW <- clip_transect(lasIn, c(0,0), c(-10,0), .5)
tE <- clip_transect(lasIn, c(0,0), c(10,0), .5)

writeLAS(tN, paste0(outDir, "transect/", "ptx1_N.las"))
writeLAS(tS, paste0(outDir, "transect/", "ptx1_S.las"))
writeLAS(tW, paste0(outDir, "transect/", "ptx1_W.las"))
writeLAS(tE, paste0(outDir, "transect/", "ptx1_E.las"))



