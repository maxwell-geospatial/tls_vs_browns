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
  trainctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
  metOut <- data.frame(siteName=character(), 
                       replicate=numeric(), 
                       set=character(), 
                       variable=character(), 
                       rmse=numeric(),
                       rsq=numeric(),
                       mae=numeric(),
                       ccc=numeric(),
                       rpd=numeric(),
                       rpiq=numeric())
  siteName <- site
  dataIn <- data %>% filter(SiteN == site)
  pltN <- 1
  fblCol <- 24
  fblNCol <- 8
  fblSCol <- 11
  fblECol <- 14
  fblWCol <- 17
  hourCol <- 25
  h1 <- 20
  h10 <- 21
  h100 <- 22
  h1000 <- 23
  
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
  
  p_hr1 <- dataIn[,c(h1, pCols)]
  t_hr1 <- dataIn[,c(h1, tCols)]
  p_hr10 <- dataIn[,c(h10, pCols)]
  t_hr10 <- dataIn[,c(h10, tCols)]
  p_hr100 <- dataIn[,c(h100, pCols)]
  t_hr100 <- dataIn[,c(h100, tCols)]
  p_hr1000 <- dataIn[,c(h1000, pCols)]
  t_hr1000 <- dataIn[,c(h1000, tCols)]
  
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
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                      tuneLength = 10,
                      ntree=500,
                      trControl = trainctrl,
                      metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
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
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
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
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
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
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
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
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "h1"
    train <- p_hr1 %>% sample_frac(.7)
    test <- setdiff(p_hr1, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "h1"
    train <- t_hr1 %>% sample_frac(.7)
    test <- setdiff(t_hr1, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "h10"
    train <- p_hr10 %>% sample_frac(.7)
    test <- setdiff(p_hr10, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "h10"
    train <- t_hr10 %>% sample_frac(.7)
    test <- setdiff(t_hr10, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "h100"
    train <- p_hr100 %>% sample_frac(.7)
    test <- setdiff(p_hr100, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "h100"
    train <- t_hr100 %>% sample_frac(.7)
    test <- setdiff(t_hr100, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Plot"
    variable1 = "h1000"
    train <- p_hr1000 %>% sample_frac(.7)
    test <- setdiff(p_hr1000, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  for(r in 1:replicates){
    repNum <- r
    set1 = "Transect"
    variable1 = "h1000"
    train <- t_hr1000 %>% sample_frac(.7)
    test <- setdiff(t_hr1000, train)
    train <- data.frame(train)
    test <- data.frame(test)
    m1 <- train(x=train[,2:ncol(train)], y=train[,1], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,1]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
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
    m1 <- train(x=train[,3:ncol(train)], y=train[,2], method = "ranger", 
                tuneLength = 10,
                ntree=500,
                trControl = trainctrl,
                metric="RMSE")
    p1 <- predict(m1, test)
    rmse1 <- yardstick::rmse_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    rsq1 <- yardstick::rsq_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    ccc1 <- yardstick::ccc_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    mae1 <- yardstick::mae_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    rpd1 <- yardstick::rpd_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    rpiq1 <- yardstick::rpiq_vec(truth=as.numeric(test[,2]), estimate=as.numeric(p1))
    result1 <- data.frame(siteName=site, 
                          replicate=repNum, 
                          set=set1, 
                          variable=variable1, 
                          rmse=rmse1,
                          rsq=rsq1,
                          ccc=ccc1,
                          rpd=rpd1,
                          rpiq=rpiq1)
    metOut <- bind_rows(metOut, result1)
  }
  return(metOut)
}

out_scfmf <- modelFunc("SCFMF", tls_grd, 50)
out_gapdr <- modelFunc("GAPDR", tls_grd, 50)
out_flsmr <- modelFunc("FLSMR", tls_grd, 50)

all_data2 <- bind_rows(out_scfmf, out_gapdr)
all_data2 <- bind_rows(all_data2, out_flsmr)

write.csv(all_data2, "D:/tls_browns_results.csv")
