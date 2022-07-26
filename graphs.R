library(tidyverse)
library(cowplot)
library(Rmpfr)

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

myR <- read.csv("C:/Maxwell_Data/Dropbox/ABGSL/fuel_load/summer_plots_study/tls_browns_results.csv")

myR %>% filter(variable == "FBD") %>% group_by(siteName, set) %>% summarize(med = mean(rsq))
myR %>% filter(variable == "FBD") %>% group_by(siteName, set) %>% summarize(med = mean(rmse))

grndDSub <- grnD %>% group_by(SiteN) %>% sample_n(20)
sub2 <- as.data.frame(grndDSub)
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMax = max(fbS, fbN, fbE, fbW))
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMin = min(fbS, fbN, fbE, fbW))
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMn = mean(fbS, fbN, fbE, fbW))

sub2 <- sub2 %>% select(Plot_name, SiteN, fbS, fbN, fbE, fbW)
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMax = max(fbS, fbN, fbE, fbW))
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMin = min(fbS, fbN, fbE, fbW))
sub2 <- sub2 %>% rowwise() %>% mutate(fbdMn = (fbS+fbN+fbE+fbW)/4)

Rmpfr::pmax(c(sub2$FbS, sub2$FbN, sub2$FbE, sub2$FbW))
pmax(sub2$FbN, sub2$FbN, sub2$FbE, sub2$FbW)

gF <-sub2 %>% filter(SiteN == "FLSMR") %>% ggplot(aes(x=Plot_name, y=fbN))+
  geom_errorbar(aes(ymax=fbdMax, ymin=fbdMin))+
  geom_point()+
  geom_point(aes(y=fbS))+
  geom_point(aes(y=fbE))+
  geom_point(aes(y=fbW))+
  geom_point(aes(y=fbdMn), color="Red")+
  scale_x_discrete(labels = c(as.character(seq(1,20, by=1))))+
  scale_y_continuous(expand=c(0.01, 0.01), breaks = c(seq(0, 140, by = 10)))+
  theme_classic()+
  labs(x="Site", y="FBD (cm)")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

gG <- sub2 %>% filter(SiteN == "GAPDR") %>% ggplot(aes(x=Plot_name, y=fbN))+
  geom_errorbar(aes(ymax=fbdMax, ymin=fbdMin))+
  geom_point()+
  geom_point(aes(y=fbS))+
  geom_point(aes(y=fbE))+
  geom_point(aes(y=fbW))+
  geom_point(aes(y=fbdMn), color="Red")+
  scale_x_discrete(labels = c(as.character(seq(1,20, by=1))))+
  scale_y_continuous(expand=c(0.01, 0.01), breaks = c(seq(0, 140, by = 10)))+
  theme_classic()+
  labs(x="Site", y="FBD (cm)")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

gS <- sub2 %>% filter(SiteN == "SCFMF") %>% ggplot(aes(x=Plot_name, y=fbN))+
  geom_errorbar(aes(ymax=fbdMax, ymin=fbdMin))+
  geom_point()+
  geom_point(aes(y=fbS))+
  geom_point(aes(y=fbE))+
  geom_point(aes(y=fbW))+
  geom_point(aes(y=fbdMn), color="Red")+
  scale_x_discrete(labels = c(as.character(seq(1,20, by=1))))+
  scale_y_continuous(expand=c(0.01, 0.01), breaks = c(seq(0, 140, by = 10)))+
  theme_classic()+
  labs(x="Site", y="FBD (cm)")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))
  
plot_grid(gF, gG, gS, nrow= 3)



inFBD <- ggplot(tls_grd, aes(y=fbd, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="FBD (cm)")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 25, 50, 75, 100, 125, 150), limits= c(0, 150))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))
  

inHr <- ggplot(tls_grd, aes(y=hr, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="Count")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 100, 200, 300, 400, 500), limits= c(0, 500))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

inH1 <- ggplot(tls_grd, aes(y=h1, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="Count")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 100, 200, 300, 400, 500), limits= c(0, 500))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

inH10 <- ggplot(tls_grd, aes(y=h10, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="Count")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 50, 100, 150, 200), limits= c(0, 200))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

inH100 <- ggplot(tls_grd, aes(y=h100, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="Count")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 25, 50, 75), limits= c(0, 75))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

inH1000 <- ggplot(tls_grd, aes(y=h1000, x=SiteN, fill=SiteN))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  labs(x="Site", y="Count")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 5, 10, 15, 20, 25, 30), limits= c(0, 30))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))


plot_grid(inFBD, inHr, inH1, inH10, inH100, inH1000, nrow= 3)

rFBR <- myR %>% filter(variable == "FBD") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (cm)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), limits= c(0, 45))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rHr <- myR %>% filter(variable == "Hourly") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (count)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), limits= c(0, 200))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1 <- myR %>% filter(variable == "h1") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (count)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 25, 50, 75, 100, 125), limits= c(0, 125))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH10 <- myR %>% filter(variable == "h10") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (count)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 10, 20, 30, 40, 50), limits= c(0, 50))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH100 <- myR %>% filter(variable == "h100") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (count)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 4, 8, 12, 16, 20), limits= c(0, 20))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1000 <- myR %>% filter(variable == "h1000") %>% ggplot(aes(y=rmse, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE (count)", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, 2, 4, 6, 8, 10, 12), limits= c(0, 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

plot_grid(rFBR, rHr, rH1, rH10, rH100, rH1000, nrow= 3)

myR %>% filter(variable == "FBD") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()

myR %>% filter(variable == "Hourly") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()

myR %>% filter(variable == "h1") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()

myR %>% filter(variable == "FBD") %>% group_by(siteName) %>% summarize(mn = mean(rmse))

rsqSummary <- myR %>% filter(variable == "FBD") %>% group_by(siteName, variable, set) %>% summarize(mn = mean(rsq))


metSummary <- myR %>% group_by(siteName, variable, set) %>% summarize(rmseMn = mean(rmse), rmseMd = median(rmse), rmseSD = sd(rmse), rsqMn = mean(rsq), rsqMd = median(rsq), rsqSD = sd(rsq))
write.csv(metSummary, "D:/aggregated_summary.csv")


rFBR <- myR %>% filter(variable == "FBD") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rHr <- myR %>% filter(variable == "Hourly") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1 <- myR %>% filter(variable == "h1") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH10 <- myR %>% filter(variable == "h10") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH100 <- myR %>% filter(variable == "h100") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1000 <- myR %>% filter(variable == "h1000") %>% ggplot(aes(y=rsq, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))


plot_grid(rFBR, rHr, rH1, rH10, rH100, rH1000, nrow= 3)


rFBR <- myR %>% filter(variable == "FBD") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rHr <- myR %>% filter(variable == "Hourly") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1 <- myR %>% filter(variable == "h1") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH10 <- myR %>% filter(variable == "h10") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH100 <- myR %>% filter(variable == "h100") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

rH1000 <- myR %>% filter(variable == "h1000") %>% ggplot(aes(y=ccc, x=set, fill=siteName))+
  geom_boxplot()+
  theme_classic()+
  labs(y="R-Squared", x="TLS Subset", fill="Site")+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  scale_y_continuous(expand = c(0, 0),breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), limits= c(0, 1))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.title = element_text(size=14), axis.text=element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))


plot_grid(rFBR, rHr, rH1, rH10, rH100, rH1000, nrow= 3)



