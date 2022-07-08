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


#Define function to calculate metrics for height bins
heightBinMetrics = function(la, b1, b2, b3, b4, bMax, prefix){
  #Input ground classified and height normalized data
  #Clip to plot radius extent, commented here because data are already clipped
  #laClip <- clip_circle(la, 0, 0, circleRad)
  #Filter out ground classified points
  lg <- filter_poi(la, Classification == 2) 
  #Filter out not ground classified points
  lng <- filter_poi(la, Classification != 2)
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

  #Get l1 metrics
  l1_cnt <- length(l1@data$Z)
  l1_per <- (l1_cnt/(not_ground_cnt))*100 
  l1_mean <- mean(l1@data$Z, na.rm="TRUE") 
  l1_median <- median(l1@data$Z, na.rm="TRUE") 
  l1_std <- sd(l1@data$Z, na.rm="TRUE") 
  l1_skew <- skewness(l1@data$Z, na.rm="TRUE")
  l1_kurt <- kurtosis(l1@data$Z, na.rm="TRUE")
  
  #Get l2 metrics
  l2_cnt <- length(l2@data$Z)
  l2_per <- (l2_cnt/(not_ground_cnt))*100 
  l2_mean <- mean(l2@data$Z, na.rm="TRUE") 
  l2_median <- median(l2@data$Z, na.rm="TRUE") 
  l2_std <- sd(l2@data$Z, na.rm="TRUE") 
  l2_skew <- skewness(l2@data$Z, na.rm="TRUE")
  l2_kurt <- kurtosis(l2@data$Z, na.rm="TRUE")
  
  #Get l3 metrics
  l3_cnt <- length(l3@data$Z) 
  l3_per <- (l3_cnt/(not_ground_cnt))*100
  l3_mean <- mean(l3@data$Z, na.rm="TRUE") 
  l3_median <- median(l3@data$Z, na.rm="TRUE") 
  l3_std <- sd(l3@data$Z, na.rm="TRUE") 
  l3_skew <- skewness(l3@data$Z, na.rm="TRUE")
  l3_kurt <- kurtosis(l3@data$Z, na.rm="TRUE")
  
  #Get L4 metrics
  l4_cnt <- length(l4@data$Z) 
  l4_per <- (l4_cnt/(not_ground_cnt))*100 
  l4_mean <- mean(l4@data$Z, na.rm="TRUE") 
  l4_median <- median(l4@data$Z, na.rm="TRUE") 
  l4_std <- sd(l4@data$Z, na.rm="TRUE")
  l4_skew <- skewness(l4@data$Z, na.rm="TRUE")
  l4_kurt <- kurtosis(l4@data$Z, na.rm="TRUE")
  
  #Get L5 metrics
  l5_cnt <- length(l5@data$Z) 
  l5_per <- (l5_cnt/(not_ground_cnt))*100 
  l5_mean <- mean(l5@data$Z, na.rm="TRUE") 
  l5_median <- median(l5@data$Z, na.rm="TRUE") 
  l5_std <- sd(l5@data$Z, na.rm="TRUE") 
  l5_skew <- skewness(l5@data$Z, na.rm="TRUE")
  l5_kurt <- kurtosis(l5@data$Z, na.rm="TRUE")
  
  #Metrics for all non-ground returns
  ng_cnt <- length(lng@data$Z) 
  ng_per <- (ng_cnt/(not_ground_cnt))*100 
  ng_mean <- mean(lng@data$Z, na.rm="TRUE")
  ng_median <- median(lng@data$Z, na.rm="TRUE")
  ng_std <- sd(lng@data$Z, na.rm="TRUE")
  ng_skew <- skewness(lng@data$Z, na.rm="TRUE")
  ng_kurt <- kurtosis(lng@data$Z, na.rm="TRUE")
  ng_quart <- unname(quantile(lng@data$Z, na.rm="TRUE", probs=c(.1,.2,.3,.4,.5,.6, .7,.8,.9)))
  q1 <- ng_quart[1]
  q2 <- ng_quart[2]
  q3 <- ng_quart[3]
  q4 <- ng_quart[4] 
  q5 <- ng_quart[5]
  q6 <- ng_quart[6]
  q7 <- ng_quart[7]
  q8 <- ng_quart[8]
  q9 <- ng_quart[9]
  
  #Merge metrics to dataframe
  all_metrics <- as.data.frame(cbind(ground_cnt, not_ground_cnt, per_ground, ng_tgi, ng_vari,
                                     ng_cnt, ng_per, ng_mean, ng_median, ng_std, ng_skew,ng_kurt,
                                     q1, q2, q3, q4, q5, q6, q7, q8, q9,
                                     l1_cnt, l1_per, l1_mean, l1_median, l1_std, l1_skew,l1_kurt,
                                     l2_cnt, l2_per, l2_mean, l2_median, l2_std, l2_skew,l2_kurt,
                                     l3_cnt, l3_per, l3_mean, l3_median, l3_std, l3_skew,l3_kurt,
                                     l4_cnt, l4_per, l4_mean, l4_median, l4_std, l4_skew,l4_kurt,
                                     l5_cnt, l5_per, l5_mean, l5_median, l5_std, l5_skew,l5_kurt))
  

  
  #Rename all columns
  cNames <- names(all_metrics)
  cNames2 <- paste0(prefix, cNames)
  names(all_metrics) <- cNames2
  
  return(all_metrics) 
}




#Set input and output directories
inDir <- "D:/ptx/out/for_comparison/all_subsets/clip/"
outDir <- "D:/ptx/out/for_comparison/all_subsets/"

files <-list.files(pattern=".las", inDir)  

b1 <- 0.5
b2 <- 1
b3 <- 1.5
b4 <- 2
bMax <- 50

for(i in 1:length(files)){
  file <- files[i]
  lasP <- readLAS(paste0(inDir, file))
  lasD <- readLAS(paste0(outDir, "donut/", file))
  lasT <- readLAS(paste0(outDir, "transects/", file))
  lasN <- readLAS(paste0(outDir, "siteN/", file))
  lasS <- readLAS(paste0(outDir, "siteS/", file))
  lasW <- readLAS(paste0(outDir, "siteW/", file))
  lasE <- readLAS(paste0(outDir, "siteE/", file))
  inPTX <- fread(paste0("D:/ptx/in/", substr(file, 1, nchar(file)-4), ".ptx"), 
                 skip = 10, 
                 header = FALSE)
  colnames(inPTX)[1:4] <- c("x","y","z", "Intensity")
  gapPulses <- inPTX %>% filter(x==0 & y==0 & z==0)
  returnPulses <- inPTX %>% filter(x!=0 & y!=0 & z!=0)
  per_gap <- as.data.frame((nrow(gapPulses)/nrow(inPTX))*100)
  names(per_gap) <- "per_gap"
  per_pulse_plot <- as.data.frame((length(lasP@data$X)/nrow(inPTX))*100)
  names(per_pulse_plot) <- "per_pulse_plot"
  per_return_plot <- as.data.frame((length(lasP@data$X)/nrow(returnPulses))*100)
  names(per_return_plot) <- "per_return_plot"
  clipNG <- filter_poi(lasP, Classification != 2)
  per_return_plot_NG <- as.data.frame((length(clipNG@data$X)/nrow(returnPulses))*100)
  names(per_return_plot_NG) <- "NGper_return_plot"
  per_pulse_plot_NG <- as.data.frame((length(clipNG@data$X)/nrow(inPTX))*100)
  names(per_pulse_plot_NG) <- "NGper_pulse_plot"
  allMets <- as.data.frame(file)
  names(allMets) <- "plot"
  allMets <- bind_cols(allMets, per_gap)
  allMets <- bind_cols(allMets, per_pulse_plot)
  allMets <- bind_cols(allMets, per_return_plot)
  allMets <- bind_cols(allMets, per_return_plot_NG)
  allMets <- bind_cols(allMets, per_pulse_plot_NG)
  pMet <- heightBinMetrics(lasP, b1, b2, b3, b4, bMax, "p_")
  dMet <- heightBinMetrics(lasD, b1, b2, b3, b4, bMax, "d_")
  tMet <- heightBinMetrics(lasT, b1, b2, b3, b4, bMax, "t_")
  nMet <- heightBinMetrics(lasN, b1, b2, b3, b4, bMax, "n_")
  sMet <- heightBinMetrics(lasS, b1, b2, b3, b4, bMax, "s_")
  eMet <- heightBinMetrics(lasE, b1, b2, b3, b4, bMax, "e_")
  wMet <- heightBinMetrics(lasW, b1, b2, b3, b4, bMax, "w_")
  allMets <- bind_cols(allMets,pMet)
  allMets <- bind_cols(allMets, dMet)
  allMets <- bind_cols(allMets, tMet)
  allMets <- bind_cols(allMets, nMet)
  allMets <- bind_cols(allMets, sMet)
  allMets <- bind_cols(allMets, eMet)
  allMets <- bind_cols(allMets, wMet)
  write.csv(allMets, paste0(outDir, "metrics/", substr(file, 1, nchar(file)-4), ".csv"))
}