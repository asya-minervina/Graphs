library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

# 1395BA  light blue
# 0D3C55 dark blue
# F16C20 orange
# A2B86C green
# C02E1D	red
# EBC844 yellow

my_six_colors<-c("#1395BA","#0D3C55", "#F16C20", "#A2B86C", "#C02E1D", "#EBC844")

CD4_clones<-fread("CD4_clones.csv")
CD8_clones<-fread("CD8_clones.csv")
CD8_clones<-as.data.frame(melt.data.table(as.data.table(CD8_clones)))
CD4_clones<-as.data.frame(melt.data.table(as.data.table(CD4_clones)))
CD8_clones<-select(CD8_clones, CD8, variable, value)
CD8_clones$CD<-"CD8"
CD4_clones$CD<-"CD4"
colnames(CD8_clones)<-c("donor", "variable", "value", "CD")
colnames(CD4_clones)<-c("donor", "variable", "value", "CD")
CD_clones<-rbind(CD4_clones, CD8_clones)
View(CD_clones)
CD_clones$variable<-str_replace_all(CD_clones$variable, "d", "")
CD_clones$variable<-as.numeric(CD_clones$variable)


ggplot(CD_clones, aes(variable, value, group=donor, color=donor))+
  geom_line(size=1.3, aes(color=donor))+
  facet_wrap(~CD)+
  theme_light()+
  scale_x_continuous(breaks = c(0, 7,15,  45), labels = c("day 0", "day 7", "day 15", "day 45"))+
  theme(panel.grid.minor.x = element_blank())+
  theme(panel.grid.minor.y = element_blank())+
  geom_point(size=3.5, color="white", shape=21, aes(fill=donor))+
  scale_fill_manual(values = my_six_colors)+
  scale_color_manual(values = my_six_colors)+
  theme(panel.grid.major.x = element_line(linetype="dashed"))+
  theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
  theme(strip.text = element_text(colour="black", face="bold"))+
  theme(legend.title = element_text(size=12, face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=10, color="black"))+
  ylab(label = "Fraction of YF17D-specific cells")+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  theme(panel.spacing.x = unit(2, "lines"))
