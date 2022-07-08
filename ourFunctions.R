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

#Function to convert PTX to LAS
ptxToLas <- function(ptx, outfile, header){
  dat <- fread(ptx, skip = 10, header = FALSE)
  names(dat) <- c("X", "Y", "Z", "Intensity", "R", "G", "B")
  dat2 <- dat %>% filter(X!=0 & Y!=0 & Z!=0)
  dat2$Intensity <- as.integer(dat2$Intensity)
  write.las(outfile, header, dat2)
}

# Convert PTX files to LAS files
clipIt <- function(las, outfile, radius1, radius2) {
  lasSub <- filter_poi(las, (sqrt(X^2+Y^2) >= radius1) & (sqrt(X^2+Y^2) <= radius2))
  writeLAS(lasSub, outfile)
}

#Function to orient plot so that placard is oriented due north. 
tlsRotate <- function(tls, x, y, z){
  tls2 <- tls
  #Convert placard center coordinates from Cartesian to spherical
  toSph <- as.data.frame(cart2sph(x=x, y=y, z=z))
  #Radian to degrees
  toSph$theta <- toSph$theta*57.2958
  toSph$phi <- toSph$phi*57.2958
  #Get theta correction factor 
  toSph <- toSph %>% mutate(theta2 = case_when(theta < 0 ~ 180 + (180-abs(theta)), TRUE ~ theta))
  rotCorrP <- 90 - toSph$theta*57.2958
  #Convert TLS coordinates from Cartesian to spherical
  toSphTLS <- as.data.frame(cart2sph(x=tls2@data$X, y=tls2@data$Y, z=tls2@data$Z))
  toSphTLS$theta <- toSphTLS$theta*57.2958
  toSphTLS$phi <- toSphTLS$phi*57.2958
  toSphTLS <- toSphTLS %>% mutate(theta2 = case_when(theta < 0 ~ 180 + (180-abs(theta)), TRUE ~ theta))
  #Apply theta correction
  toSphTLS$theta2 <- toSphTLS$theta2 + rotCorrP 
  #Convert back to Cartesian
  toSphTLS <- toSphTLS %>% mutate(theta3 = case_when(theta2 > 180 ~ -180 + abs(180-theta2), TRUE ~ theta2))
  toCarTLS <- as.data.frame(sph2cart(theta=toSphTLS$theta3/57.2958, phi=toSphTLS$phi/57.2958, r=toSphTLS$r))
  #Replace coordinates in TLS data
  tls2@data$X <- toCarTLS$x
  tls2@data$Y <- toCarTLS$y
  tls2@data$Z <- toCarTLS$z
  #Write result
  return(tls2)
}

heightBinMetrics = function(la, circleRad, b1, b2, b3, b4, bMax){
  #Input ground classified and height normalized data
  #Clip to plot radius extent
  laClip <- clip_circle(la, 0, 0, circleRad)
  #Filter out ground classified points
  lg <- filter_poi(laClip, Classification == 2) 
  #Filter out not ground classified points
  lng <- filter_poi(laClip, Classification != 2)
  #Break data into five height bins (defined using b1 through b4 variables)
  l1 <- filter_poi(lng, Z > 0 & Z <= b1) 
  l2 <- filter_poi(lng, Z > b1 & Z <= b2) 
  l3 <- filter_poi(lng, Z > b2 & Z <= b3) 
  l4 <- filter_poi(lng, Z > b3 & Z <= b4) 
  l5 <- filter_poi(lng, Z > b4 & Z <= bMax)
  
  #Get plot scale metrics
  #Count of ground points
  ground_cnt <- length(lg@data$Z) 
  #Count of not ground points
  not_ground_cnt <- length(lng@data$Z) 
  #Percent of returns that are ground
  per_ground <- (ground_cnt/(not_ground_cnt+ground_cnt))*100 
  #Calculate triangular greenness index from RGB data
  ng_tgi <- mean(((670-480)*(as.double(lng@data$R)-as.double(lng@data$G)))-
                   ((670-550)*(as.double(lng@data$R)-as.double(lng@data$B)))/2, na.rm=TRUE)/(2^15)
  #Calculate Visual Atmospheric Resistance Index from RGB data
  ng_vari <- mean((as.double(lng@data$G)-as.double(lng@data$R))/
                    (as.double(lng@data$G)+as.double(lng@data$R)-as.double(lng@data$B)+.0001))
  
  #Get l1 metrics
  l1_cnt <- length(l1@data$Z)
  l1_per <- (l1_cnt/(not_ground_cnt))*100 
  l1_mean <- mean(l1@data$Z) 
  l1_median <- median(l1@data$Z) 
  l1_std <- sd(l1@data$Z) 
  l1_tgi <- mean(((670-480)*(as.double(l1@data$R)-as.double(l1@data$G)))-
                   ((670-550)*(as.double(l1@data$R)-as.double(l1@data$B)))/2, na.rm=TRUE)/(2^15)
  l1_vari <- mean((as.double(l1@data$G)-as.double(l1@data$R))/
                    (as.double(l1@data$G)+as.double(l1@data$R)-as.double(l1@data$B)+0.0001))
  l1_skew <- skewness(l1@data$Z)
  l1_kurt <- kurtosis(l1@data$Z)
  
  #Get l2 metrics
  l2_cnt <- length(l2@data$Z)
  l2_per <- (l2_cnt/(not_ground_cnt))*100 
  l2_mean <- mean(l2@data$Z) 
  l2_median <- median(l2@data$Z) 
  l2_std <- sd(l2@data$Z) 
  l2_tgi <- mean(((670-480)*(as.double(l2@data$R)-as.double(l2@data$G)))-
                   ((670-550)*(as.double(l2@data$R)-as.double(l2@data$B)))/2, na.rm=TRUE)/(2^15)
  l2_vari <- mean((as.double(l2@data$G)-as.double(l2@data$R))/
                    (as.double(l2@data$G)+as.double(l2@data$R)-as.double(l2@data$B)+0.0001))
  l2_skew <- skewness(l2@data$Z)
  l2_kurt <- kurtosis(l2@data$Z)
  
  #Get l3 metrics
  l3_cnt <- length(l3@data$Z) 
  l3_per <- (l3_cnt/(not_ground_cnt))*100
  l3_mean <- mean(l3@data$Z) 
  l3_median <- median(l3@data$Z) 
  l3_std <- sd(l3@data$Z) 
  l3_tgi <- mean(((670-480)*(as.double(l3@data$R)-as.double(l3@data$G)))-
                   ((670-550)*(as.double(l3@data$R)-as.double(l3@data$B)))/2, na.rm=TRUE)/(2^15)
  l3_vari <- mean((as.double(l3@data$G)-as.double(l3@data$R))/
                    (as.double(l3@data$G)+as.double(l3@data$R)-as.double(l3@data$B)+0.0001))
  l3_skew <- skewness(l3@data$Z)
  l3_kurt <- kurtosis(l3@data$Z)
  
  #Get L4 metrics
  l4_cnt <- length(l4@data$Z) 
  l4_per <- (l4_cnt/(not_ground_cnt))*100 
  l4_mean <- mean(l4@data$Z) 
  l4_median <- median(l4@data$Z) 
  l4_std <- sd(l4@data$Z)
  l4_tgi <- mean(((670-480)*(as.double(l4@data$R)-as.double(l4@data$G)))-
                   ((670-550)*(as.double(l4@data$R)-as.double(l4@data$B)))/2, na.rm=TRUE)/(2^15)
  l4_vari <- mean((as.double(l4@data$G)-as.double(l4@data$R))/
                    (as.double(l4@data$G)+as.double(l4@data$R)-as.double(l4@data$B)+0.0001))
  l4_skew <- skewness(l4@data$Z)
  l4_kurt <- kurtosis(l4@data$Z)
  
  #Get L5 metrics
  l5_cnt <- length(l5@data$Z) 
  l5_per <- (l5_cnt/(not_ground_cnt))*100 
  l5_mean <- mean(l5@data$Z) 
  l5_median <- median(l5@data$Z) 
  l5_std <- sd(l5@data$Z) 
  l5_tgi <- mean(((670-480)*(as.double(l5@data$R)-as.double(l5@data$G)))-
                   ((670-550)*(as.double(l5@data$R)-as.double(l5@data$B)))/2, na.rm=TRUE)/(2^15)
  l5_vari <- mean((as.double(l5@data$G)-as.double(l5@data$R))/
                    (as.double(l5@data$G)+as.double(l5@data$R)-as.double(l5@data$B)+0.0001))
  l5_skew <- skewness(l5@data$Z)
  l5_kurt <- kurtosis(l5@data$Z)
  
  #Cloud-level metrics from lidR package cloud_metrics function
  standard_metrics<-cloud_metrics(lng,.stdmetrics_z) %>% do.call('cbind',.) 
  
  #Merge metrics to dataframe
  all_metrics <- as.data.frame(cbind(filename, ground_cnt, not_ground_cnt, per_ground, ng_tgi, ng_vari,
                                     l1_cnt, l1_per, l1_mean, l1_median, l1_std, l1_tgi, l1_vari, l1_skew,l1_kurt,
                                     l2_cnt, l2_per, l2_mean, l2_median, l2_std, l2_tgi, l2_vari, l2_skew,l2_kurt,
                                     l3_cnt, l3_per, l3_mean, l3_median, l3_std, l3_tgi, l3_vari, l3_skew,l3_kurt,
                                     l4_cnt, l4_per, l4_mean, l4_median, l4_std, l4_tgi, l4_vari, l4_skew,l4_kurt,
                                     l5_cnt, l5_per, l5_mean, l5_median, l5_std, l5_tgi, l5_vari, l5_skew,l5_kurt,
                                     standard_metrics))
  
  
  #Remove entropy metric
  drop <- c("zentropy")
  all_metrics = all_metrics[,!(names(all_metrics) %in% drop)]
  
  #Rename all columns
  cNames <- names(all_metrics)
  cNames2 <- paste0("h_", cNames)
  names(all_metrics) <- cNames2
  
  return(all_metrics) 
}


#Get header info (This is from Jeff)
ptx.header<-function(input_file){
  list(
    name = input_file,
    col_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, 
                                  header = FALSE)[1,]),
    row_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, 
                                  header = FALSE)[2,]),
    scan_center = as.matrix(read.csv(input_file, sep = "", skip = 2,nrows = 1, 
                                     header = FALSE)),
    reg = as.matrix(read.csv(input_file, sep = "", skip = 3,nrows = 3, 
                             header = FALSE)),
    trans_matrix = as.matrix(read.csv(input_file, sep = "", skip = 6,nrows = 4, 
                                      header = FALSE))
  )
}

#get proportion of gaps
perGapFunc <- function(dat){
  data_gap <- dat %>% filter(x==0 & y==0 & z==0)
  per_gap <- (nrow(data_gap)/nrow(dat))*100
  return(data.frame(per_gap))
}

#noise filter
nnFilter2 = function(las, d = 0.05, n = 2){
  rnn = knn(las %>% las2xyz, k = n+1)$nn.dists[,-1]
  keep = rep(T, nrow(las@data))
  for(i in 1:ncol(rnn)){
    keep = keep & rnn[,i] < d
  }
  las = filter_poi(las, keep)
  return(las)
}

