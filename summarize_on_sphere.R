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

#Declaring constants
height = 2 #height of the scanner 
thetaDiv=.5 #Angle of theta (side-to-side) covered by each voxel
phiDiv=.5 #Angle of phi (up-and-down) covered by each voxel
r_blank <- raster(ncol= 360/thetaDiv, nrow = 180/phiDiv, 
                  xmn=-180, xmx=180,ymn=-90,ymx=90) #Create an empty raster to use as a template fo the results


cld_fun<-function(r,tlsData){
  #Aim: Calculate coordinates where each ray would intersect the sphere
  # create a dataframe to calculate distances to the scanner
  cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height)) #Extract x,y,z coordinates
  cld$a <- cld$x^2 + cld$y^2 + cld$z^2 #Cacluate distance from scanner to point using the Pythagorean Theorem
  #Math to calculate where ray would intersect a sphere with a given radius
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
  
  #Calculate theta and phi (must convert from radians to degrees)
  s2$theta <- s$theta*(180/pi)
  s2$phi <- s$phi*(180/pi)
  p2$theta <- p$theta*(180/pi)
  p2$phi <- p$phi*(180/pi)
  
  #Convert to a spatial object
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
  # Aim: Transforms a point cloud into a raster showing the points that passed or do not pass a certain radius.
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
      
      gridPre <- gridS - gridP # make temporary raster of difference. This represents number of points remaining once returns from prior bins are subtracted.
      
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
      out2 <- stack(out2, out) # Stack raster grids of n radius 
    }
  }#Close radius loop
  names(out2) <- as.character(radii)# name layers
  return(out2) # set function output. 
} 

#Convert back to a df with x,y,z coordinates and theta, phi, r at center of each voxel
#data column represents the proportion of pulses passing through the voxel that returned form it
#9999 indicates an occluded voxel
summarizeSpheres <- function(rast, bins){
  rast[is.na(test)]<- 9999 #Convert  NA to 9999
  asPnts <- as.data.frame(rasterToPoints(rast)) #Convert first band of raster to a dataframe
  asPnts2 <- asPnts[,c(1:3)] #Extract out only the needed columns
  asPnts2$r <- .5 # Add radius
  names(asPnts2) <- c("t", "p", "data", "r")
  for(i in 2:bins){#Loop through the remaining columns to add to the point data sets
    asPntsx <- asPnts[,c(1:2, i+2)] 
    asPntsx$r <- i - .5 #Define radius at center of voxel
    names(asPntsx) <- c("t", "p", "data", "r")
    asPnts2 <- bind_rows(asPnts2, asPntsx)
  }
  
  #Calculate the areas of each voxel
  asPnts2$area <- ((pi*((asPnts2$r+.5)^3)*(4/3)) - (pi*((asPnts2$r-.5)^3)*(4/3)))/(ncol(rast)*nrow(rast))
  
  #Convert from spherical coordinates back to Cartesian
  asPnts2$X2 <- asPnts2$r*cos(asPnts2$t*0.01745329)*cos(asPnts2$p*0.01745329)
  asPnts2$Y2 <- asPnts2$r*sin(asPnts2$t*0.01745329)*cos(asPnts2$p*0.01745329)
  asPnts2$Z2 <- asPnts2$r*sin(asPnts2$p*0.01745329)
  #Add Cartesian coordinates back to the table
  #Add height of the scanner back to the points so coordinates match the original point cloud
  asPnts2$Z2 <- asPnts2$Z + 2
  #return(sphere_metrics)
  return(asPnts2)
}

#Read in a test point cloud
setwd("D:/")
las1 <- readLAS("101_laz_frames_40_to_40.las")

#Apply function to obtain the raster of proportion of returns at different distances
test <- on_sphere_multi_v2(tlsData=las1, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

#Write out the raster (not necessary, just for visualization/to check output)
#writeRaster(test, "D:/testout2.tif")

#Obtain table of results
summary1 <- summarizeSpheres(test,26)

#Post-Process table
names(summary1) <- c("t", "p", "Intensity", "r", "Area", "X", "Y", "Z") #Rename columns
summary1 <- mutate(summary1, dist = sqrt(X^2+Y^2)) #Calculate distance from sensor in X/Y plane
#Augment intensity data to that data are integer and occlusion is coded to 101
summary1$Intensity <- summary1$Intensity*100 
summary1$Intensity <- as.integer(summary1$Intensity)
summary1 <- summary1 %>% mutate(Intensity2 = ifelse(Intensity == 999900, 101, Intensity))
#Extract only points within 10 m of the sensor in the X/Y plane
summary2 <- summary1 %>% filter((X<10 & X >-10) & (Y<10 & Y >-10) & (Z >= 0 & Z <= 20))
#Subset needed fields
summary3 <- summary2 %>% dplyr::select(X, Y, Z, Intensity2)
names(summary3) <- c("X", "Y", "Z", "Intensity")
summary3$Intensity <- as.integer(summary3$Intensity)
#Can write out las, but not necessary (will need to drop area column and rename Intesity2 to Intensity)
#write.las("FILENAME", lhead, summary3)

#Calculate Features
occluded <- summary3 %>% filter(Intensity == 101)
returned <- summary3 %>% filter(Intensity >0 & Intensity <101)
returned$Intensity <- as.integer(returned$Intensity)
gaps <- summary3 %>% filter(Intensity == 0)
write.las("D:/returned.las", lhead, returned)
write.las("D:/occluded.las", lhead, occluded)
write.las("D:/gaps.las", lhead, gaps)

per_occluded = ((sum(occluded$Area))/(sum(summary3$Area)))*100

per_returned = ((sum(returned$Area))/(sum(summary3$Area)))*100

per_gaps = ((sum(gaps$Area))/(sum(summary3$Area)))*100

