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

#color palette two #440154FF or #482677FF and #29AF7FFF or #3CBB75FF
#red and grey c( "#2B2D42", "#D90429")

#function for beautiful labels
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

#Barplots for diff between day0 and other days, place name of the file in the first function

make_barplot <- function(x, plotname, uplim=1600) {
  num_DE_clones<-fread(x)
  num_DE_clones<-melt(num_DE_clones)
  
  
  ggplot(num_DE_clones, aes(variable, value, group=V1))+
    geom_bar(aes(fill = V1), stat = "identity", position = "dodge", width=0.8, color="black", size=0.2)+
    scale_fill_manual(values = my_six_colors, name="donor\n ")+
    theme_light()+
    scale_y_continuous(expand=c(0,0), limits = c(0, uplim))+
    theme(panel.grid.major.x = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=18, color="black"))+
    ylab(label = "Number of expanded clones\ncomparing to day 0")+
    theme(axis.title.y = element_text(size=22))+
    scale_x_discrete(labels = c("day -7", "day 7", "day 15", "day 45"))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=20, color="black"))+
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
    geom_line(size=2.2, aes(color=donor))+
    facet_wrap(~CD)+
    theme_light()+
    scale_x_continuous(breaks = c(0, 7,15,  45), labels = c("0", "7", "15", "45"))+
    scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), labels = c(0, 0.01, 0.02, 0.03, 0.04))+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    theme(panel.grid.major.x = element_line(linetype="dashed"))+
    geom_point(size=4, color="white", shape=21, aes(fill=donor))+
    theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
    theme(strip.text = element_text(colour="black", face="bold", size=20))+
    ylab(label = "Fraction of YF17D-specific cells")+
    theme(axis.title.x = element_text(size=22))+
    xlab("Days after vaccination")+
    theme(panel.spacing.x = unit(2, "lines"))+
    scale_color_manual(values = my_six_colors)+
    scale_fill_manual(values = my_six_colors)+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=20, color="black"))+
    theme(axis.title.y = element_text(size=22))+
    theme(axis.text.x = element_text(size=18, color="black"))+
    theme(axis.text.y = element_text(size=16, color="black"))+
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
    scale_color_manual(values = c( "#3D9637", "#F9C80E"), labels=c("non-twin pairs", "twin pairs"))+
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


#Plots for similarity

make_similarity_plot <- function(filename, plotname) {
  load(filename)
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
  data3$variable<-str_replace_all(string = data3$variable, pattern = "V1", replacement = "0 mismatches")
  data3$variable<-str_replace_all(string = data3$variable, pattern = "V2", replacement = "1 mismatch")
  data3$variable<-str_replace_all(string = data3$variable, pattern = "V3", replacement = "2 mismatches")
  
  ggplot(data3, aes(wat, log10(value), fill=wat))+
    geom_pointrange(data=data3, aes(wat, log10(value), ymin=ymin, ymax=ymax),position = position_jitter(width = 0.4, height = 0), size=1, shape=21, color="#052F5F")+
    theme_light()+
    theme(panel.grid.major.x = element_blank())+
    facet_wrap(~variable)+
    scale_fill_manual(values = c("#F9C80E", "#052F5F" ), labels=c("real\ndata", "simulated\ndata"), name=" name")+
    theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
    theme(strip.text = element_text(colour="black", face="bold", size=17))+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=18, color="black"))+
    #theme(axis.text.x = element_blank())+
    theme(legend.position = "none")+
    theme(axis.text.y = element_text(size=20, color="black"))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=20, color="black"))+
    ylab(label = "Title y")+
    theme(axis.title.y = element_text(size=22, face="bold"))+
  ggtitle(plotname)+
    theme(plot.title = element_text(hjust = 0.5, size=20))+
    scale_x_discrete(labels=c("real\ndata", "simulated\ndata"))
}

#color palette two #440154FF or #482677FF and #29AF7FFF or #3CBB75FF    "#3D9637"
#make histogram plot
make_histplot <- function(histdatacsv2) {
  histplot<-read.csv2(histdatacsv2)
  histplot<-select(histplot, 2:4)
  histplot2<-data.frame(mismatches=c(histplot$mismatches, histplot$mismatches), Freq=c(histplot$YF_LLW.Freq, histplot$CMV_NLV.Freq),antigen=c(rep("YF",times=21),rep("CMV",times=21)),stringsAsFactors = F)
  ggplot(histplot2, aes(mismatches, Freq, fill=antigen))+geom_histogram(stat="identity", position = "dodge", color="black")+
    theme_light()+
    theme(panel.grid.major.x = element_blank())+
    theme(panel.grid.minor.x = element_blank())+
    scale_fill_manual(values = c("#052F5F", "#F9C80E"), name="Legend\nname")+
    theme(legend.title = element_text(size=22, face="bold"))+
    #theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=22, color="black"))+
    theme(axis.text.y = element_text(size=22, color="black"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2, 'lines'))+
    ylab(label = "Title y")+
    xlab(label="Title x")+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=24, face="bold"))+
    theme(axis.title.x = element_text(size=24, face="bold"))+
    scale_y_continuous(expand=c(0,0), limits=c(0, 101))+
    scale_x_continuous( breaks = c(0:12), limits=c(0,13))
}

#Heatmap
make_heatmap <- function(data_heat_csv, plotname) {
  heatmap_gr<-read.csv(data_heat_csv)
  sht<-as.matrix(heatmap_gr[, 2:7])
  rownames(sht)<-heatmap_gr[,1]
  row.names(sht)<-gsub("_0_F1","",row.names(sht))
  colnames(sht)<-gsub("_0_F1","",colnames(sht))
  print(sht)
  x <- sht
  dd.col <- as.dendrogram(hclust(dist(x)))
  dd.row <- as.dendrogram(hclust(dist(t(x))))
  dx <- dendro_data(dd.row)
  dy <- dendro_data(dd.col)
  ggdend <- function(df) {
    ggplot() +
      geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend), size=1.2) +
      labs(x = "", y = "") + theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), plot.margin = unit(c(0,7,0,3), "cm"))
  }
  px <- ggdend(dx$segments)
  py <- ggdend(dy$segments) + coord_flip()
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  xx <- sht[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$car <- xx_names[[1]]
  df$car <- with(df, factor(car, levels=car, ordered=TRUE))
  mdf <- reshape2::melt(df, id.vars="car")
  p <- ggplot(mdf, aes(x = car, y = variable)) + geom_tile(aes(fill = value))+theme_minimal()+
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(expand=c(0,0))+
    theme(plot.margin = unit(c(0,1,0,1), "cm"))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    theme(legend.text = element_text(colour="black", size = 24,  hjust = 1))+
    theme(legend.title=element_text(colour="black", size = 24, face = "bold"))+
    theme(legend.text.align=0)+
    theme(legend.key = element_rect(size = 22),legend.key.size = unit(2, 'lines'))
    #scale_fill_viridis( na.value = "grey60",breaks=c(0, 1e-5, 2e-5, 3e-5), labels=c(" 0", scientific_10(1e-5), scientific_10(2e-5),scientific_10(3e-5)),
                        #name="Probability\n ", direction = 1)
  library(gridExtra)
  library(grid)
  pp1<-grid.arrange(px, p, ncol=1, heights=c(2,6), top=plotname)
  pp1
}

#Scatterplot

make_scatter <- function(DR_data_csv) {
  plot<-read.csv(DR_data_csv)
  ggplot(plot, aes(log10(DRhi.Read.proportion), log10(CD8.Read.proportion), fill=changed))+
    geom_point(size=3.8, color="black", shape=21, alpha=0.75)+
    theme_light()+
    theme(panel.grid.minor.y = element_blank())+
    theme(panel.grid.minor.x = element_blank())+
    scale_fill_manual(values = c("#052F5F", "#F9C80E"), labels=c("not changed\nclones\n ", "changed clones"), name="Legend\nname")+
    theme(legend.title = element_text(size=22, face="bold"))+
    #theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=22, color="black"))+
    theme(axis.text.y = element_text(size=22, color="black"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2, 'lines'))+
    ylab(label = "Title y")+
    xlab(label="Title x")+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=24, face="bold"))+
    theme(axis.title.x = element_text(size=24, face="bold"))+
      geom_abline(slope = 1, intercept = 0, size=1.3)
}

#Roc
roc<-read.csv("fig_data/leave_one_out_roc.csv")
roc_tet<-read.csv("fig_data/tet_roc.csv")
ggplot(roc, aes(1-specificity, sensitivity,col=donor))+
  geom_line(size=1.5, aes(color=donor))+
  theme_light()+
  scale_x_continuous(expand=c(0,0), limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
  scale_y_continuous(expand=c(0,0), limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
  theme(panel.grid.minor.x = element_blank())+
  theme(panel.grid.minor.y = element_blank())+
  xlab(label = "1-Specificity")+
  theme(axis.title.x = element_text(size=22))+
  ylab("Sensitivity")+
  scale_color_manual(values = c(my_six_colors))+
  theme(legend.title = element_text(size=24, face="bold"))+
  theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
  theme(legend.text = element_text(size=20, color="black"))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.x = element_text(size=18, color="black"))+
  theme(axis.text.y = element_text(size=16, color="black"))+
  geom_abline(slope = 1, intercept = 0, size=1)
  
  ggtitle(plotname)+
  theme(plot.title = element_text(hjust = 0.5, size=20))
  