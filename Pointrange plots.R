library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(gridExtra)
library(grid)
library(stringi)
library(scales)
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

load("shVJ_tbl.rda")

shVJ_tbl$pos<-"a"
shVJ_tbl$min<-shVJ_tbl$Freq-shVJ_tbl$SD
shVJ_tbl$max<-shVJ_tbl$Freq+shVJ_tbl$SD
shVJ_tbl$dn<-paste0(shVJ_tbl$Var1, "/", shVJ_tbl$Var2)
shVJ_tbl<-arrange(shVJ_tbl, -Freq)
shVJ_tbl$num<-c(1:15)
#shVJ_tbl$pos<-1+rnorm(n = nrow(shVJ_tbl), mean = 1, sd = 0.5)
#shVJ_tbl$pos<-runif(n = nrow(shVJ_tbl), min = 0, max=1)
 View(shVJ_tbl)
#set.seed(1)
#pd<-position_jitter(width=0.1, height=0)

 
#Graph with pointrange
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
  scale_x_continuous(breaks=c(1:15), labels = shVJ_tbl$dn)+
  theme(panel.grid.major.x = element_blank())+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(face="bold", size=10, color="black", angle = 90))+
  #theme(legend.background = element_rect(size=0.8, linetype="solid", colour = "black"))+
  theme(legend.text = element_text(size=10, face="bold"))+
  scale_y_continuous(labels=c("0", scientific_10(1e-5), scientific_10(2e-5),scientific_10(3e-5),scientific_10(4e-5)))+
  theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))




  
#Same but barblot
ggplot(shVJ_tbl, aes(x=num, y=Freq))+
  geom_bar(aes(fill = twins), stat = "identity", position = "dodge", width=0.6, color="black", size=0.8)+
  scale_fill_manual(values = c( "#8D99AE", "#D90429"), labels=c("non-twin pairs", "twin pairs"))+
  theme_light()+
  theme(panel.grid.minor.x = element_blank())+
  theme(panel.grid.minor.y = element_blank())+
  theme(axis.title.x = element_blank())+
  ylab(label = "Frequency")+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  geom_errorbar(data=shVJ_tbl, aes(ymin=min, ymax=max),  width=0.3, size=1)+
  scale_x_continuous(breaks=c(1:15), labels = shVJ_tbl$dn)+
  theme(panel.grid.major.x = element_blank())+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(face="bold", size=10, color="black", angle = 90))+
  #theme(legend.background = element_rect(size=0.8, linetype="solid", colour = "black"))+
  theme(legend.text = element_text(size=10, face="bold"))
  


