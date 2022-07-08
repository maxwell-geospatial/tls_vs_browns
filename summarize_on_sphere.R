library(rlas)
library(dplyr)
library(ggplot2)
library(plotly)
library(geometry)
library(MASS)
library(raster)
library(sampSurf)
library(rgdal)
library(raster)
library(lidR)
library(tmap)
library(sf)
library(Morpho)


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


#Define radius of sphere
radii <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)

#Declaring constants ----
height = 2 #height of the scanner 
thetaDiv=.5 
phiDiv=.5 
grdDist=.25 
r_blank <- raster(ncol= 360/thetaDiv, nrow = 180/phiDiv, 
                  xmn=-180, xmx=180,ymn=-90,ymx=90) 


cld_fun<-function(r,tlsData){
  #Aim: Calculate coordinates where each ray would intersect the sphere
  # create a dataframe to calculate distances to the scanner
  cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height))
  cld$a <- cld$x^2 + cld$y^2 + cld$z^2
  cld$t1 <- -sqrt(0^2-(4*cld$a*(-r^2)))/(2*cld$a)
  cld$t2 <- sqrt(0^2-(4*cld$a*(-r^2)))/(2*cld$a)
  cld$x1 <- cld$x*cld$t1
  cld$y1 <- cld$y*cld$t1
  cld$z1 <- cld$z*cld$t1
  cld$x2 <- cld$x*cld$t2
  cld$y2 <- cld$y*cld$t2
  cld$z2 <- cld$z*cld$t2
  cld$slp <- cld$z/(sqrt(cld$x^2+ cld$y^2))
  cld$slp1 <- cld$z1/(sqrt(cld$x1^2+ cld$y1^2))
  cld$slp2 <- cld$z2/(sqrt(cld$x2^2+ cld$y2^2))
  
  #Find correct coordinate based on slope
  cld<-cld %>% mutate(xs = case_when(sign(slp)==sign(slp1)~ x1,sign(slp)!=sign(slp1)~x2),
                      ys = case_when(sign(slp)==sign(slp1)~ y1,sign(slp)!=sign(slp1)~y2),
                      zs = case_when(sign(slp)==sign(slp1)~ z1,sign(slp)!=sign(slp1)~z2))
  
}

sphere_func<-function(cld=cld){
  #Aim: calculate spherical coordinates. 
  #Calculate Euclidean distance from sensor to point
  cld$dist_act <- sqrt(cld$x^2 + cld$y^2 + cld$z^2)
  #Calculate Euclidean distance from sensor to sphere/point intersection
  cld$dist_sphere <- sqrt(cld$xs^2 + cld$ys^2 + cld$zs^2)
  
  #Filter out points that actually pass through the sphere
  pass <- cld %>% filter(dist_act >= dist_sphere)

  #Convert from Cartesian to spherical coordinates
  s <- as.data.frame(cart2sph(x=cld$xs, 
                              y=cld$ys, 
                              z=cld$zs))
  p <- as.data.frame(cart2sph(x=pass$xs, 
                              y=pass$ys, 
                              z=pass$zs))
  
  s2 <- s[,c(1,2)]
  p2 <- p[,c(1,2)]
  
  s2$theta <- s$theta*(180/pi)
  s2$phi <- s$phi*(180/pi)
  p2$theta <- p$theta*(180/pi)
  p2$phi <- p$phi*(180/pi)
  
  s3 <- st_multipoint(x=as.matrix(s2), dim="XY") # make spatial object
  p3 <- st_multipoint(x=as.matrix(p2), dim="XY") # make spatial object
  
  #create lists and finalize function
  multipoint_list<-list(s3,p3)
  names(multipoint_list)<-c("s3","p3")
  return(multipoint_list)
} 

pointCount <- function(r, pts){
  # Idea source: https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r
  # Aim: Calculate the number of points in a raster cell grid. 
  
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

on_sphere_multi_v2 <- function(tlsData, radii, height, thetaDiv=1, phiDiv=1){
  # Aim: Transforms a point cloud into a raster showing the points that passed or not pass a certain radius, 
  # depending on the layer of the raster.
  for(ra in radii){ # Begin for radius loop 
    if(ra == radii[1]){ # Logic statement of 1st radius
      r <- ra #read radius
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates 
      
      gridS <- pointCount(r_blank,list_point$s3) #create point object of points on the sphere
      gridP <- pointCount(r_blank, list_point$p3)#create point object of points that pass the sphere
      
      null1 <- raster::calc(gridS, function(x){x[x<=0]<-NA; return(x)}) #If count in cell is zero, recode to NA since no returns were recorded (canopy gap)
      null2 <- null1 > 0 ## Make a 0/1 mask for NA vs. Value
      
      out <- (1-(gridP/gridS))*null2 #create raster for first radius
    }else if(ra == radii[2]){# logical statement 2nd radius
      r<-ra # read radius
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates 
      
      gridPre <- gridS - gridP # make temporal raster of difference. This represents number of points remaining once returns from prior bins are subtracted.
      
      gridS <- pointCount(r_blank, list_point$s3) - gridPre # create point object of points on the sphere minus the previous sphere
      gridP <- pointCount(r_blank, list_point$p3) #create point object of points that pass the sphere
      
      mask1 <- out # read previous layer. 
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)}) #If proportion of returns in prior bin is 1 then recode to NA since no returns passed to current bin
      mask3 <- mask2 >= 0 #Make a 0/1 mask for NA vs. cells with values
      
      out <- (1-(gridP/gridS)) # make faster for 2nd radius
      out2 <- stack(mask1, out) # stack raster grids
    }else{ # logic statement n_th radius
      r<-ra # read data
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates
      
      gridPre <- gridS - gridP# make temporal raster of difference.
      
      gridS <- pointCount(r_blank, list_point$s3) - gridPre # create point object of points on the sphere minus the previous sphere
      gridP <- pointCount(r_blank, list_point$p3) #create point object of points that pass the sphere
      
      mask1 <- out # read raster of previous radius.
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)}) ##If proportion of returns in prior bin is 1 then recode to NA since no returns passed to current bin
      mask3 <- mask2 >= 0 #Make a 0/1 mask for NA vs. cells with values
      
      out <- (1-(gridP/gridS))*mask3 # create raster of n radius. 
      out2 <- stack(out2, out) # stack raster grids of n radius 
    }
  }#Close radius loop
  names(out2) <- as.character(radii)# name layers
  return(out2) # set function output. 
} 

# Applying functions ----

test <- on_sphere_multi_v2(tlsData=tlsData, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

writeRaster(test, "D:/sphere_vox1111.tif") # save data. 

summarizeSpheres <- function(rast, bins){
  rast[is.na(test)]<- 9999
  asPnts <- as.data.frame(rasterToPoints(rast))
  asPnts2 <- asPnts[,c(1:3)]
  asPnts2$r <- .5
  names(asPnts2) <- c("t", "p", "data", "r")
  for(i in 2:bins){
    asPntsx <- asPnts[,c(1:2, i+2)] 
    asPntsx$r <- i - .5
    names(asPntsx) <- c("t", "p", "data", "r")
    asPnts2 <- bind_rows(asPnts2, asPntsx)
  }
  
  toCart <- as.data.frame(sph2cart(theta=asPnts2$t, 
                              phi=asPnts2$p, 
                              r=asPnts2$r))
  asPnts2 <- bind_cols(asPnts2, toCart)
  asPnts2$z <- asPnts2$z + 2
  
  lng <- asPnts2 %>% filter((z > 0 & z <= 20) & (x >= -10 & x <= 10) & (y >= -10 & y <= 10))
  l1 <- lng %>% filter(z > 0 & z <= 0.5) #Subset all non-ground samples between 0.001 and 0.5 m
  l2 <- lng %>% filter(z > 0.5 & z <= 1) #Subset all non-ground samples between 0.5 and 1 m
  l3 <- lng %>% filter(z > 1 & z <= 1.5) #Subset all non-ground samples between 1 and 1.5 m
  l4 <- lng %>% filter(z > 1.5 & z <= 2) #Subset all non-ground samples between 1.5 and 2 m
  l5 <- lng %>% filter(z > 2) #Subset all non-ground samples above 2 m
  
  lng_na_per <- nrow(lng %>% filter(data == 9999))/nrow(lng)
  lng_zero_per <- nrow(lng %>% filter(data == 0))/nrow(lng)
  lng_prop_per <- 1 - lng_na_per - lng_zero_per
  
  l1_na_per <- nrow(l1 %>% filter(data == 9999))/nrow(l1)
  l1_zero_per <- nrow(l1 %>% filter(data == 0))/nrow(l1)
  l1_prop_per <- 1 - l1_na_per - l1_zero_per
  
  l2_na_per <- nrow(l2 %>% filter(data == 9999))/nrow(l2)
  l2_zero_per <- nrow(l2 %>% filter(data == 0))/nrow(l2)
  l2_prop_per <- 1 - l2_na_per - l2_zero_per
  
  l3_na_per <- nrow(l3 %>% filter(data == 9999))/nrow(l3)
  l3_zero_per <- nrow(l3 %>% filter(data == 0))/nrow(l3)
  l3_prop_per <- 1 - l3_na_per - l3_zero_per
  
  l4_na_per <- nrow(l4 %>% filter(data == 9999))/nrow(l4)
  l4_zero_per <- nrow(l4 %>% filter(data == 0))/nrow(l4)
  l4_prop_per <- 1 - l4_na_per - l4_zero_per
  
  l5_na_per <- nrow(l5 %>% filter(data == 9999))/nrow(l5)
  l5_zero_per <- nrow(l5 %>% filter(data == 0))/nrow(l5)
  l5_prop_per <- 1 - l5_na_per - l5_zero_per
  
  sphere_metrics <- as.data.frame(cbind(lng_na_per,lng_zero_per,lng_prop_per,l1_na_per,l1_zero_per,l1_prop_per,l2_na_per,l2_zero_per,l2_prop_per,l3_na_per,l3_zero_per,l3_prop_per,l4_na_per,l4_zero_per,l4_prop_per,l5_na_per,l5_zero_per,l5_prop_per))
  
  #return(sphere_metrics)
  return(asPnts2)
}

setwd("D:/cylinder_sets")
las1 <- readLAS("set_0_frames_1_to_1.las")
test <- on_sphere_multi_v2(tlsData=las1, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)
summary1 <- summarizeSpheres(test,26)

names(summary1) <- c("R", "G", "Intensity", "B", "X", "Y", "Z")
summary1$id <- as.numeric(row.names(summary1))
summary1 <- mutate(summary1, dist = sqrt(X^2+Y^2))
summary2 <- summary1[,c(5:7,3)]
summary2$Intensity <- summary2$Intensity*100
summary2$Intensity <- as.integer(summary2$Intensity)
summary2 <- summary2 %>% mutate(Intensity2 = ifelse(Intensity == 999900, 101, Intensity))
summary3 <- summary2 %>% filter((X<10 & X >-10) & (Y<10 & Y >-10) & (Z > -2 & Z <20))
summary4 <- summary3 %>% dplyr::select(X, Y, Z, Intensity2)
names(summary4) <- c("X", "Y", "Z", "Intensity")
summary4$Intensity <- as.integer(summary4$Intensity)
write.las("D:/sphere_points2222.las", lhead, summary4)
lasGrid <- readLAS("D:/sphere_points2222.las")
l1 <- filter_poi(lasGrid, Intensity <101 & Intensity > 0) 
writeLAS(l1, "D:/occlusions2222.las")


xseq <- seq(-10, 10, by=.2) 
yseq <- seq(-10, 10, by=.2) 
zseq <- seq(-2, 20, by=.2) 

regPnts <- expand.grid(xseq, yseq, zseq)
names(regPnts) <- c("X", "Y", "Z")
regPnts <- mutate(regPnts, dist = sqrt(X^2+Y^2))
regPntsClip <- regPnts %>% filter(dist <= 10)
regPntsClip <- regPntsClip[,1:3]
write.las("D:/grid_points3.las", lhead, regPntsClip)


withSphere <- mcNNindex(as.matrix(summary2), as.matrix(regPnts), k=1)
regPnts$id <- withSphere

regPnts2 <- left_join(regPnts, summary1, by="id") 
regPnts3 <- regPnts2 %>% dplyr::select(X.x, Y.x, Z.x, Intensity)
regPnts3$Intensity <- as.integer(regPnts3$Intensity)
names(regPnts3) <- c("X", "Y", "Z", "Intensity")

write.las("D:/grid_points2.las", lhead, regPnts3)

for(i in 1:199){
  las1 <- readLAS(paste0("set_", as.character(i), "_frames_1_to_1.las"))
  test <- on_sphere_multi_v2(tlsData=las1, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)
  summary2 <- summarizeSpheres(test, 26)
  summary1 <- bind_rows(summary1, summary2)
}


lasGrid <- readLAS("D:/sphere_points2.las")
l1 <- filter_poi(lasGrid, Intensity <101 & Intensity > 0) 
writeLAS(l1, "D:/occlusions2.las")
