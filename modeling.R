library(dplyr)
library(caret)
library(randomForest)
library(Metrics)
library(yardstick)

#Data Prep===================================================================================================

#Read in ground data table
grnD <- read.csv("C:/Maxwell_Data/Dropbox/ABGSL/fuel_load/summer_plots_study/USFWS_2021_field_data.csv")

grnD <- grnD[,1:23]

nmLst <- c("Plot_name", "SiteN", "PlotN", "Date", "Collectors", "Time", "GPS", 
           "fbN", "lN", "dN",
           "fbS", "lS", "dS",
           "fbW", "lW", "dW",
           "fbE", "lE", "dE",
           "h1", "h10", "h100", "h1000")

names(grnD) <- nmLst

grnD <- grnD %>% filter(Plot_name != 'GAPDR_0109_20210607')

grnD$fbS <- as.integer(grnD$fbS)

grnD <- grnD %>% mutate(fbd = ((fbN+fbS+fbE+fbW)/4))
grnD <- grnD %>% mutate(hr = h1+h10+h100+h1000)

grnD %>% group_by(SiteN) %>% count()

#Read in and merge metrics tables
inDir <- "D:/tls_browns/metrics/metrics/"

files <- list.files(inDir, pattern = '\\.csv$', full.names = TRUE)

all_data <- do.call(rbind, lapply(files, function(x) 
  transform(read.csv(x), File = basename(x))))

all_data$Plot_name <- substr(all_data$plot, 1, nchar(all_data$plot)-4)

tls_grd <- inner_join(grnD, all_data, by="Plot_name")

tls_grd %>% group_by(SiteN) %>% count()

tls_grd[is.na(tls_grd)] <- 0



#Train and Predict Function =========================================================================
modelFunc <- function(site, data, replicates){
  metOut <- data.frame(siteName=character(), 
                       replicate=numeric(), 
                       set=character(), 
                       variable=character(), 
                       rmse=numeric(),
                       rsq=numeric())
  siteName <- site
  dataIn <- data %>% filter(SiteN == site)
  pltN <- 1
  fblCol <- 24
  fblNCol <- 8
  fblSCol <- 11
  fblECol <- 14
  fblWCol <- 17
  hourCol <- 25
  
  pCols <-c(28:86)
  dCols <- c(87:140)
  tCols <- c(141:194)
  nCols <- c(195:248)
  sCols <- c(249:302)
  eCols <- c(303:356)
  wCols <- c(357:410)
  
  p_fbl <- dataIn[,c(fblCol, pCols)]
  d_fbl <- dataIn[,c(fblCol, dCols)]
  t_fbl <- dataIn[,c(fblCol, tCols)]
  p_hr <- dataIn[,c(hourCol, pCols)]
  t_hr <- dataIn[,c(hourCol, tCols)]
  s_fblN <- dataIn[,c(pltN, fblNCol, nCols)]
  s_fblS <- dataIn[,c(pltN, fblSCol, sCols)]
  s_fblE <- dataIn[,c(pltN, fblECol, eCols)]
  s_fblW <- dataIn[,c(pltN, fblWCol, wCols)]
  reName <- names(s_fblN)
  names(s_fblS) <- reName
  names(s_fblE) <- reName
  names(s_fblW) <- reName
  s_fbl <- bind_rows(s_fblN, s_fblS)
  s_fbl <- bind_rows(s_fbl, s_fblE)
  s_fbl <- bind_rows(s_fbl, s_fblW)
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "FBD"
    train <- p_fbl %>% sample_frac(.7)
    test <- setdiff(p_fbl, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,2:ncol(train)], y=train[,1], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,1], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Donut"
    variable1 = "FBD"
    train <- d_fbl %>% sample_frac(.7)
    test <- setdiff(d_fbl, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,2:ncol(train)], y=train[,1], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,1], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "FBD"
    train <- t_fbl %>% sample_frac(.7)
    test <- setdiff(t_fbl, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,2:ncol(train)], y=train[,1], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,1], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "Hourly"
    train <- p_hr %>% sample_frac(.7)
    test <- setdiff(p_hr, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,2:ncol(train)], y=train[,1], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,1], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "Hourly"
    train <- t_hr %>% sample_frac(.7)
    test <- setdiff(t_hr, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,2:ncol(train)], y=train[,1], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,1], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Site"
    variable1 = "FBD"
    siteLst <- as.data.frame(levels(as.factor(s_fbl$Plot_name)))
    names(siteLst) <- "S"
    trainSites <- siteLst %>% sample_frac(.7)
    trainSites2 <- trainSites$S
    train <- s_fbl %>% filter(Plot_name %in% trainSites2)
    test <- setdiff(s_fbl, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- randomForest(x=train[,3:ncol(train)], y=train[,2], ntree=1000)
    p1 <- predict(m1, test)
    rmse1 <- Metrics::rmse(test[,2], p1)
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1)
    metOut <- bind_rows(metOut, result1)
  }
  return(metOut)
}

out_scfmf <- modelFunc("SCFMF", tls_grd, 50)
out_scfmf <- modelFunc("SCFMF", tls_grd, 50)
out_scfmf <- modelFunc("SCFMF", tls_grd, 50)


