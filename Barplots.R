library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(ggthemes)
library(tibble)
num_DE_clones<-fread("num_DE_clones.csv")
View(num_DE_clones)
num_DE_clones<-melt(num_DE_clones)

  # ggplot(num_DE_clones, aes(num, value))+
  # geom_bar(aes(fill = variable), stat = "identity", position = "dodge", width=0.2, color="black")+
  # #scale_fill_brewer(palette = "Set1")+
  # theme_light()+
  # facet_wrap(~donor)+ylim(0, 1000)+
  # scale_y_continuous(expand=c(0,0), limits = c(0, 1000))+
  # scale_fill_manual(values = c("red2", "dodgerblue3", "seagreen4", "darkorange"))+
  # theme(axis.text.x = element_blank())+
  # theme(axis.title.x= element_blank())+
  # theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
  # theme(strip.text = element_text(colour="black", face="bold"))+
  # theme(panel.grid.major.y = element_blank())
  # 
 


num_DE_clones$num<-rep("aaa", times=ncol(num_DE_clones))
#good palletes det1, dark2, tableu10medium, tableau10
ggplot(num_DE_clones, aes(variable, value, group=donor))+
  geom_bar(aes(fill = donor), stat = "identity", position = "dodge", width=0.6, color="black", size=0.2)+
  scale_fill_manual(values =my_six_colors)+
  theme_light()+
  scale_y_continuous(expand=c(0,0), limits = c(0, 1000))+
  theme(panel.grid.major.x = element_blank())+
#scale_fill_tableau("tableau10medium")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=10, color="black"))+
  ylab(label = "Number of expanded clones comparing to day 0")+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  scale_x_discrete(labels = c("day -7", "day 7", "day 15", "day 45"))+
  theme(legend.title = element_text(size=12, face="bold"))
  
  
load("similarity.rda")
data<-as.data.frame(data)
data$wat<-rep("data", times=ncol(data))
simul_means<-as.data.frame(simul_means)
simul_means$wat<-rep("simul", times=ncol(simul_means))
data2<-rbind(data, simul_means)
data2<-rownames_to_column(data2)
data2$rowname<-gsub(pattern = "_15_F1", replacement = "", x = data2$rowname)
simul_sd<-as.data.frame(simul_sd)
data2<-as.data.frame(melt(data2))
data2$sd<-rep(0, times=nrow(data2))
data2$sd[data2$wat=="simul"&data2$variable=="V1"]<-simul_sd$V1
data2$sd[data2$wat=="simul"&data2$variable=="V2"]<-simul_sd$V2
data2$sd[data2$wat=="simul"&data2$variable=="V3"]<-simul_sd$V3
data3<-data2
data3$value=data3$value/(0.5*simul_nums*(simul_nums-1))[data3$rowname]
data3$sd=data3$sd/(0.5*simul_nums*(simul_nums-1))[data3$rowname]
data3$ymin<-log10(data3$value-data3$sd/10)
data3$ymax<-log10(data3$value+data3$sd/10)

ggplot(data3[data3$variable=="V1", ], aes(rowname, log10(value), group=wat, fill=wat))+
  geom_bar(stat = "identity", position = "dodge", width=0.6, color="black", size=0.2)+
  theme_light()+
  #scale_y_continuous(expand=c(0,0), limits = c(0, 45))+
  theme(panel.grid.major.x = element_blank())+
  scale_fill_manual(values = c("#2B2D42", "#D90429"), labels=c("real data", "simulated data"))+
  #facet_wrap(~variable)
  #geom_errorbar(data=simul_sd, aes(ymin=sd, ymax=sd),  width=0.3, size=1)
  [data3$variable=="V1", ]
  [data3$variable=="V3", ]
[data3$variable=="V3", ]

ggplot(data3, aes(wat, log10(value), color=wat))+
  geom_pointrange(data=data3, aes(wat, log10(value), ymin=ymin, ymax=ymax),position = position_jitter(width = 0.2, height = 0), size=0.75)+
  #geom_point()+
  theme_light()+
  #scale_y_continuous(expand=c(0,0), limits = c(0, 45))+
  theme(panel.grid.major.x = element_blank())+
  facet_wrap(~variable)+
  scale_color_manual(values = c("#2B2D42", "#D90429"), labels=c("real data", "simulated data"))+
  theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
  theme(strip.text = element_text(colour="black", face="bold"))+
  theme(legend.title = element_text(size=12, face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=10, color="black"))+
    theme(axis.text.y = element_text(face="bold", size=10, color="black"))+
  ylab(label = "log10")+
  theme(axis.title.y = element_text(size=12, face="bold"))


