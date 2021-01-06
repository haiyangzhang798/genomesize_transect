rm(list=ls())
setwd("/Users/haiyangzhang/haiyangwork/GitHub/genomesize_transect")

#getwd()
# reading in all the files ------------------------------------------------
# community
comm <- read.csv("Data/2020_COMM_data_median_all.csv")
comm_large <- read.csv("Data/2020_COMM_data_median_2.5_large.csv")
comm_small <- read.csv("Data/2020_COMM_data_median_2.5_small.csv")

# locaiton and environmental data
space <- read.csv("Data/2020_spatial_site.csv", header=T)
env <- read.csv("Data/2020_env_site.csv", header = T)

# load in the funciton for generating the gdm data for making gdm figures
# first we only focus on the turnover component
source("R/func_gdm_data_generation.R")
gdm_final_all_turn <- func_gdm_run(comm = comm, space = space, env = env, dis_select = "turnover")
gdm_final_large_turn <- func_gdm_run(comm = comm_large, space = space, env = env, dis_select = "turnover")
gdm_final_small_turn <- func_gdm_run(comm = comm_small, space = space, env = env, dis_select = "turnover")


##########################################################################################
PSAMPLE <- 200
preddata <- rep(0,times=PSAMPLE)

overlayX <- seq(from=min(gdm_final_all_turn$ecological), to=max(gdm_final_all_turn$ecological), length=PSAMPLE )
overlayY <- 1 - exp( - overlayX )
d1_overlay <- data.frame(overlayX=overlayX, overlayY=overlayY)
tab1 <- as.data.frame(cbind(gdm_final_all_turn$ecological,gdm_final_all_turn$observed))
colnames(tab1) <- c("ecological","observed")

library(ggplot2);
g_total1 <- ggplot(tab1,aes(y=observed, x=ecological))+
  geom_point(pch=20, cex=1, col="grey")+
  scale_x_continuous(expand = c(0,0),limits = c(0.3,3.5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  geom_line(aes(overlayX, overlayY), d1_overlay,lwd=2,col="black")+
  xlab("Environmental disimilarity")+ylab(expression(paste(beta[turnover])))+theme_bw()+
  theme(text=element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),axis.text = element_text(size = 28),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0.37, y=1.05, label= "A", size = 10)

### second figure: first we need to make the data table
IMP_SPECIES_ALL_data <- cbind(gdm_final_large_turn$IMP_SPECIES_ALL,gdm_final_small_turn$IMP_SPECIES_ALL)
colnames(IMP_SPECIES_ALL_data) <- c("Large", "Small")
IMP_SPECIES_ALL_data$Predictor <- row.names(IMP_SPECIES_ALL_data)

IMP_SPECIES_ALL_data_figure <- subset(IMP_SPECIES_ALL_data, Predictor%in%c("Geographic", "soilP", "Alt", "MAP","MAT", "soilN"))
library(plyr)
IMP_SPECIES_ALL_data_figure$Predictor <- revalue(IMP_SPECIES_ALL_data_figure$Predictor, c("Geographic"="Distance","MAP" = "MAP","MAT" = "MAT" ,"soilP"="Soil P", "Alt"="Altitude", "soilN"= "Soil N"))
library(dplyr)
IMP_SPECIES_ALL_data_figure <- IMP_SPECIES_ALL_data_figure %>% 
                               arrange(factor(Predictor, levels =c("Distance","MAT","Altitude","MAP","Soil N", "Soil P")))


library(tidyr)
data_long <- gather(IMP_SPECIES_ALL_data_figure, GS_group, Relative_importance, Large:Small, factor_key=TRUE)
data_long$Predictor <- factor(data_long$Predictor, levels =c("Distance","MAT","Altitude","MAP","Soil N", "Soil P"))
library(ggpubr)

g_var1 <-  ggplot(data=data_long, aes(x=Predictor,  y=Relative_importance))+
  geom_col(aes(color = GS_group, fill = GS_group), position = position_dodge(0.8), width = 0.7) +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"))+ylim(0, 0.25)+ 
  theme_pubr()+theme(legend.position = "none")+xlab("")+ylab("Relative importance")+
  font("xlab", size = 28)+ font("ylab", size = 28)+font("xy.text", size = 28)+
  annotate("text", x=0.5, y=0.25, label= "B", size = 10)



#tiff("Output/transect_fig2abc.tif", width = 7, height = 12, units = "in",res = 150)
library(grid)
gl = lapply(list(g_total1,g_var1), ggplotGrob)     
library(gtable)
g = do.call(rbind, c(gl, size="first"))
g$widths = do.call(unit.pmax, lapply(gl, "[[", "widths"))

#grid.newpage()
#grid.draw(g) 
#dev.off()


##########################################################################################
######################################Fig1d################################################
exSplines_sim <- isplineExtract(gdm_final_all_turn)
exSplines_sim_L <- isplineExtract(gdm_final_large_turn)
exSplines_sim_S <- isplineExtract(gdm_final_small_turn)

tab1_sim <- as.data.frame(exSplines_sim[[1]])
tab2_sim <- as.data.frame(exSplines_sim[[2]])
names(tab1_sim) <- paste(names(tab1_sim),".x", sep="")
names(tab2_sim) <- paste(names(tab2_sim),".y", sep="")
tab_fig1d_sim <- cbind(tab1_sim, tab2_sim)



tab1_L <- as.data.frame(exSplines_sim_L[[1]])
tab2_L <- as.data.frame(exSplines_sim_L[[2]])
names(tab1_L) <- paste(names(tab1_L),".x", sep="")
names(tab2_L) <- paste(names(tab2_L),".y", sep="")
tab_fig1d_sim_L <- cbind(tab1_L, tab2_L)


tab1_S <- as.data.frame(exSplines_sim_S[[1]])
tab2_S <- as.data.frame(exSplines_sim_S[[2]])
names(tab1_S) <- paste(names(tab1_S),".x", sep="")
names(tab2_S) <- paste(names(tab2_S),".y", sep="")
tab_fig1d_sim_S <- cbind(tab1_S, tab2_S)


###########################################plot splines
## 1 Geographic Distance
f1d_distance <- ggplot()+
  #geom_line(aes(Geographic.x, Geographic.y), tab_fig1d_sim,lwd=1.5,col="black")+
  #geom_hline(yintercept=max(tab_fig1d_sim$Geographic.y))+
  geom_line(aes(Geographic.x, Geographic.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(Geographic.x, Geographic.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab("Distance (km)")+ylab("")+theme_bw()+ylim(0, 2)+scale_x_continuous(breaks=c(0,1500, 3000))+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0, y=2, label= "C", size = 10)


#4. Aridity
f1d_MAP <- ggplot()+
  #geom_line(aes(Aridity.x, Aridity.y), tab_fig1d_sim,lwd=1.5,col="black")+
  geom_line(aes(MAP.x, MAP.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(MAP.x, MAP.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab("MAP (mm)")+ylab("")+theme_bw()+ylim(0, 2)+scale_x_continuous(breaks=c(100,250, 400))+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0.3, y=2, label= "F", size = 10)

f1d_MAT <- ggplot()+
  #geom_line(aes(Aridity.x, Aridity.y), tab_fig1d_sim,lwd=1.5,col="black")+
  geom_line(aes(MAT.x, MAT.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(MAT.x, MAT.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab(expression(MAT (degree*C)))+ylab("Partial ecological distance")+theme_bw()+ylim(0, 2)+scale_x_continuous(breaks=c(-4,0, 4, 8))+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=-3.5, y=2, label= "D", size = 10)

## 5 SoilP
f1d_SoilP <- ggplot()+
  #geom_line(aes(exp(soilP.x), soilP.y), tab_fig1d_sim,lwd=1.5,col="black")+
  geom_line(aes(soilP.x, soilP.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(soilP.x, soilP.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab("Soil P (%)")+ylab("")+theme_bw()+ylim(0, 2)+scale_x_continuous(breaks=c(0.02,0.05, 0.08))+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0, y=2, label= "H", size = 10)

## 7 Alt
f1d_Alt <- ggplot()+
  #geom_line(aes(Alt.x^2, Alt.y), tab_fig1d_sim,lwd=1.5,col="black")+
  geom_line(aes(Alt.x, Alt.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(Alt.x, Alt.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab("Altitude (m)")+ylab("")+theme_bw()+ylim(0, 2)+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0, y=2, label= "E", size = 10)




f1d_soilN <- ggplot()+
  #geom_line(aes(soilN.x, soilN.y), tab_fig1d_sim,lwd=1.5,col="black")+
  geom_line(aes(soilN.x, soilN.y), tab_fig1d_sim_L,lwd=1.5,col="#0073C2FF")+
  geom_line(aes(soilN.x, soilN.y), tab_fig1d_sim_S,lwd=1.5,col="#EFC000FF")+
  xlab("Soil N (%)")+ylab("")+theme_bw()+ylim(0, 2)+
  theme(text=element_text(size=28),axis.text = element_text(size=28),panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour = NA),
        axis.line = element_line(colour="black"), 
        panel.grid.major=element_blank())+
  annotate("text", x=0, y=2, label= "G", size = 10)



png("Output/transect_fig2.png", width = 22, height = 12, units = "in",res = 300)
library(gridExtra)
#P=multiplot(f1d_distance,f1d_Alt,f1d_MAT, f1d_PET,f1d_soilC_N, f1d_SoilP,f1d_SoilpH,f1d_soilN_P, f1d_UVB, cols=2)
grid.arrange(g, f1d_distance,f1d_MAT,f1d_MAP, f1d_soilN, f1d_SoilP, ncol = 4, 
             layout_matrix = cbind(c(1,1,1), c(1,1,1), c(2,3,4), c(5, 6)))
dev.off()

