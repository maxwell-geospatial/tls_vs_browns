grdD <- read.csv("C:/Maxwell_Data/Dropbox/ABGSL/fuel_load/summer_plots_study/USFWS_2021_field_data_Edited.csv")

library(ggplot2)
library(cowplot)

ggplot(grdD, aes(x=Site, y=totFD/4, fill=Site))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  labs(y="Fuel Bed Depth (cm)")+
  theme(axis.text.y = element_text(size=12, color="gray40"))+
  theme(axis.text.x = element_text(size=12, color="gray40"))+
  theme(plot.title = element_text(face="bold", size=18, color="gray40"))+
  theme(axis.title = element_text(size=14, color="gray40"))+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))
  

g1 <- ggplot(grdD, aes(x=Site, y=hr1, fill=Site))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  labs(y="Count of 1-Hour Fuels")+
  theme(axis.text.y = element_text(size=12, color="gray40"))+
  theme(axis.text.x = element_text(size=12, color="gray40"))+
  theme(plot.title = element_text(face="bold", size=18, color="gray40"))+
  theme(axis.title = element_text(size=14, color="gray40"))+
  theme(legend.position = "None")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))


g10 <- ggplot(grdD, aes(x=Site, y=hr10, fill=Site))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  labs(y="Count of 10-Hour Fuels")+
  theme(axis.text.y = element_text(size=12, color="gray40"))+
  theme(axis.text.x = element_text(size=12, color="gray40"))+
  theme(plot.title = element_text(face="bold", size=18, color="gray40"))+
  theme(axis.title = element_text(size=14, color="gray40"))+
  theme(legend.position = "None")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

g100 <- ggplot(grdD, aes(x=Site, y=hr100, fill=Site))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  labs(y="Count of 100-hour Fuels")+
  theme(axis.text.y = element_text(size=12, color="gray40"))+
  theme(axis.text.x = element_text(size=12, color="gray40"))+
  theme(plot.title = element_text(face="bold", size=18, color="gray40"))+
  theme(axis.title = element_text(size=14, color="gray40"))+
  theme(legend.position = "None")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))


g1000 <- ggplot(grdD, aes(x=Site, y=hr1000, fill=Site))+
  geom_violin()+
  geom_boxplot(fill="gray", width=.1)+
  theme_classic()+
  scale_fill_manual(values = c("#7F8C1F", "#D2D904", "#A67A44"))+
  labs(y="Count of 1,000-Hour Fuels")+
  theme(axis.text.y = element_text(size=12, color="gray40"))+
  theme(axis.text.x = element_text(size=12, color="gray40"))+
  theme(plot.title = element_text(face="bold", size=18, color="gray40"))+
  theme(axis.title = element_text(size=14, color="gray40"))+
  theme(legend.position = "None")+
  theme(panel.grid.major.y = element_line(colour = "gray40", linetype="dashed"))

plot_grid(g1, g10, g100, g1000, nrow= 2)

