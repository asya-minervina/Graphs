library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(ggthemes)
library(tibble)
library(gridExtra)
library(grid)
library(stringi)
library(scales)
library(readr)

#make six color pallete #1395BA-light blue; #0D3C55-dark blue; #F16C20-orange; #A2B86C-green; #C02E1D-red; #EBC844-yellow
my_six_colors<-c("#1395BA","#0D3C55", "#F16C20", "#A2B86C", "#C02E1D", "#EBC844")

#function for beautiful labels
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

#Barplots for diff between day0 and other days, place name of the file in the first function

make_barplot <- function(x, plotname, uplim) {
  num_DE_clones<-fread(x)
  num_DE_clones<-melt(num_DE_clones)
  
  
  ggplot(num_DE_clones, aes(variable, value, group=V1))+
    geom_bar(aes(fill = V1), stat = "identity", position = "dodge", width=0.6, color="black", size=0.2)+
    scale_fill_manual(values = my_six_colors, name="Donor")+
    theme_light()+
    scale_y_continuous(expand=c(0,0), limits = c(0, uplim))+
    theme(panel.grid.major.x = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(face="bold", size=10, color="black"))+
    ylab(label = "Number of expanded clones comparing to day 0")+
    theme(axis.title.y = element_text(size=12, face="bold"))+
    scale_x_discrete(labels = c("day -7", "day 7", "day 15", "day 45"))+
    theme(legend.title = element_text(size=12, face="bold"))+
    ggtitle(plotname)+
    theme(plot.title = element_text(hjust = 0.5, size=20))
}

#Lineplots for CD4 CD8 dynamics

make_lineplot <- function(CD4, CD8, plotname) {
  CD4_clones<-read.csv(CD4)
  CD8_clones<-read.csv(CD8)
  CD8_clones<-melt(CD8_clones)
  CD4_clones<-melt(CD4_clones)
  CD8_clones$CD<-"CD8"
  CD4_clones$CD<-"CD4"
  colnames(CD8_clones)<-c("donor", "variable", "value", "CD")
  colnames(CD4_clones)<-c("donor", "variable", "value", "CD")
  CD_clones<-rbind(CD4_clones, CD8_clones)
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
    theme(panel.spacing.x = unit(2, "lines"))+
    ggtitle(plotname)+
    theme(plot.title = element_text(hjust = 0.5, size=20))
}

#Make Pointrange plots about shared tcrs

make_pointrange <- function(filename, plotname) {
  
  shVJ_tbl<-read.csv(filename)
  shVJ_tbl<-shVJ_tbl[shVJ_tbl$Uniq==T, ]
  shVJ_tbl$pos<-"a"
  shVJ_tbl$min<-shVJ_tbl$Freq-shVJ_tbl$SD
  shVJ_tbl$max<-shVJ_tbl$Freq+shVJ_tbl$SD
  shVJ_tbl$dn<-paste0(shVJ_tbl$Var1, "/", shVJ_tbl$Var2)
  shVJ_tbl<-arrange(shVJ_tbl, -Freq)
  shVJ_tbl$num<-c(1:15)

  
  ggplot(shVJ_tbl, aes(x=num, y=Freq, col=twins))+
    geom_point(size=3.2)+
    scale_color_manual(values = c( "#2B2D42", "#D90429"), labels=c("non-twin pairs", "twin pairs"))+
    theme_light()+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    theme(axis.title.x = element_blank())+
    ylab(label = "Frequency")+
    theme(axis.title.y = element_text(size=12, face="bold"))+
    geom_errorbar(data=shVJ_tbl, aes(ymin=min, ymax=max),  width=0.3, size=1)+
    scale_x_continuous( breaks=c(1:15) ,labels = shVJ_tbl$dn)+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank())+
    theme(axis.title.y = element_text(size=12, face="bold"))+
    theme(axis.text.x = element_text(face="bold", size=10, color="black", angle = 90))+
    #theme(legend.background = element_rect(size=0.8, linetype="solid", colour = "black"))+
    theme(legend.text = element_text(size=10, face="bold"))+
    scale_y_continuous(labels=c("0", scientific_10(1e-5), scientific_10(2e-5),scientific_10(3e-5), scientific_10(4e-5)))+
    theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))+
    ggtitle(plotname)+
    theme(plot.title = element_text(hjust = 0.5, size=20))
}

