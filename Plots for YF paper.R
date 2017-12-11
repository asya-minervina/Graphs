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
library(ggdendro)
library(viridis)

load("C:/Users/Asya/Dropbox/YF_gr/Fig_data/donor_list_rev.rda")
d_rename<-donor_list_rev
#make six color pallete #1395BA-light blue; #0D3C55-dark blue; #F16C20-orange; #99BA46-green; #C02E1D-red; #EBC844-yellow   3C6628-dark green
my_six_colors<-c("#0D3C55","#1395BA", "#C02E1D","#F16C20","#3C6628","#99BA46")

#color palette two #440154FF or #482677FF and #29AF7FFF or #3CBB75FF
#red and grey c( "#2B2D42", "#D90429")

#function for beautiful labels
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

#Barplots for diff between day0 and other days, place name of the file in the first function
#save as 7 9 inches

make_barplot <- function(x, uplim=1600) {
  num_DE_clones<-fread(x)
  num_DE_clones<-melt(num_DE_clones)
  num_DE_clones$V1<-d_rename[num_DE_clones$V1]
  ggplot(num_DE_clones, aes(variable, value, group=V1))+
    geom_bar(aes(fill = V1), stat = "identity", position = "dodge", width=0.8, color="black", size=0.2)+
    scale_fill_manual(values = my_six_colors, name="donor")+
    theme_light()+
    scale_y_continuous(expand=c(0,0), limits = c(0, uplim))+
    theme(panel.grid.major.x = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    ylab(label = "Number of expanded clones\ncomparing to day 0")+
    theme(axis.title.y = element_text(size=26))+
    scale_x_discrete(labels = c("day -7", "day 7", "day 15", "day 45"))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
    # ggtitle(plotname)+
    # theme(plot.title = element_text(hjust = 0.5, size=20))
}

#num_DE_clones<-fread("C:/Users/Asya/Dropbox/YF_gr/Fig_data/edger_res.csv")
#num_DE_clones<-melt(num_DE_clones)
make_barplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/edger_res.csv", uplim = 1200)
make_barplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/Max_thres.csv", uplim = 1700)
#Lineplots for CD4 CD8 dynamics

make_lineplot <- function(CD4, CD8, err4, err8) {
  CD4_clones<-read.csv(CD4)
  CD8_clones<-read.csv(CD8)
  CD4err<-read.csv(err4)
  CD8err<-read.csv(err8)
  CD8err<-melt(CD8err)
  CD4err<-melt(CD4err)
  CD8err$CD<-"CD8"
  CD4err$CD<-"CD4"
  CD_err<-rbind(CD4err, CD8err)
  colnames(CD_err)<-c("donor", "var", "error", "CD")
  CD8_clones<-melt(CD8_clones)
  CD4_clones<-melt(CD4_clones)
  CD8_clones$CD<-"CD8"
  CD4_clones$CD<-"CD4"
  colnames(CD8_clones)<-c("donor", "variable", "value", "CD")
  colnames(CD4_clones)<-c("donor", "variable", "value", "CD")
  CD_clones<-rbind(CD4_clones, CD8_clones)
  CD_clones$variable<-str_replace_all(CD_clones$variable, "d", "")
  CD_clones$variable<-as.numeric(CD_clones$variable)
  CD_clones$donor<-as.character(CD_clones$donor)
  CD_clones$donor<-d_rename[CD_clones$donor]
  CD_clones$error<-CD_err$error
  CD_clones$ymin<-CD_clones$value-CD_clones$error
  CD_clones$ymax<-CD_clones$value+CD_clones$error
  
  ggplot(CD_clones, aes(variable, value, group=donor, color=donor))+
    geom_line(size=2.2, aes(color=donor))+
    geom_errorbar(data=CD_clones, aes(ymin=ymin, ymax=ymax),  width=0.5, size=1)+
    facet_wrap(~CD)+
    theme_light()+
    scale_x_continuous(breaks = c(0, 7,15,  45), labels = c("0", "7", "15", "45"))+
    scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels = c(0, 0.025, 0.05, 0.075))+
    #scale_y_continuous(breaks=c(-4, -3, -2, -1),labels=c(scientific_10(1e-4), scientific_10(1e-3),scientific_10(1e-2),scientific_10(1e-1)))+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    theme(panel.grid.major.x = element_line(linetype="dashed"))+
    geom_point(size=4, color="white", shape=21, aes(fill=donor))+
    theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
    theme(strip.text = element_text(colour="black", face="bold", size=20))+
    ylab(label = "Fraction of YF17D-specific cells")+
    theme(axis.title.x = element_text(size=26))+
    xlab("Days after vaccination")+
    theme(panel.spacing.x = unit(2, "lines"))+
    scale_color_manual(values = my_six_colors)+
    scale_fill_manual(values = my_six_colors)+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=26))+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  # ggtitle(plotname)+
  # theme(plot.title = element_text(hjust = 0.5, size=20))
}

make_log10_lineplot <- function(CD4, CD8, err4, err8) {
  CD4_clones<-read.csv(CD4)
  CD8_clones<-read.csv(CD8)
  CD4err<-read.csv(err4)
  CD8err<-read.csv(err8)
  CD8err<-melt(CD8err)
  CD4err<-melt(CD4err)
  CD8err$CD<-"CD8"
  CD4err$CD<-"CD4"
  CD_err<-rbind(CD4err, CD8err)
  colnames(CD_err)<-c("donor", "var", "error", "CD")
  CD8_clones<-melt(CD8_clones)
  CD4_clones<-melt(CD4_clones)
  CD8_clones$CD<-"CD8"
  CD4_clones$CD<-"CD4"
  colnames(CD8_clones)<-c("donor", "variable", "value", "CD")
  colnames(CD4_clones)<-c("donor", "variable", "value", "CD")
  CD_clones<-rbind(CD4_clones, CD8_clones)
  CD_clones$variable<-str_replace_all(CD_clones$variable, "d", "")
  CD_clones$variable<-as.numeric(CD_clones$variable)
  CD_clones$donor<-as.character(CD_clones$donor)
  CD_clones$donor<-d_rename[CD_clones$donor]
  CD_clones$error<-CD_err$error
  CD_clones$ymin<-CD_clones$value-CD_clones$error
  CD_clones$ymax<-CD_clones$value+CD_clones$error
 
  ggplot(CD_clones, aes(variable, log10(value), group=donor, color=donor))+
    geom_line(size=2.2, aes(color=donor))+
    geom_errorbar(data=CD_clones, aes(ymin=log10(ymin), ymax=log10(ymax)),  width=0.5, size=1)+
    facet_wrap(~CD)+
    theme_light()+
    scale_x_continuous(breaks = c(0, 7,15,  45), labels = c("0", "7", "15", "45"))+
    #scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels = c(0, 0.025, 0.05, 0.075))+
    scale_y_continuous(breaks=c(-4, -3, -2, -1),labels=c(scientific_10(1e-4), scientific_10(1e-3),scientific_10(1e-2),scientific_10(1e-1)))+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    theme(panel.grid.major.x = element_line(linetype="dashed"))+
    geom_point(size=4, color="white", shape=21, aes(fill=donor))+
    theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
    theme(strip.text = element_text(colour="black", face="bold", size=20))+
    ylab(label = "Fraction of YF17D-specific cells\n")+
    theme(axis.title.x = element_text(size=26))+
    xlab("Days after vaccination")+
    theme(panel.spacing.x = unit(2, "lines"))+
    scale_color_manual(values = my_six_colors)+
    scale_fill_manual(values = my_six_colors)+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=26))+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
    # ggtitle(plotname)+
    # theme(plot.title = element_text(hjust = 0.5, size=20))
}
CD4_clones<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4.csv")
CD8_clones<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4.csv")
CD4err<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4_err.csv")
CD8err<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8_err.csv")
make_lineplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_edgeR_CD4.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_edgeR_CD8.csv",
              "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4_err.csv","C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8_err.csv" )
make_log10_lineplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8.csv", 
              "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4_err.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8_err.csv")
make_lineplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8.csv", 
                    "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD4_err.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/phen_max_CD8_err.csv")

#Make Pointrange plots about shared tcrs

make_pointrange <- function(filename) {
  
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
    theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))
    # ggtitle(plotname)+
    # theme(plot.title = element_text(hjust = 0.5, size=20))
}


#Plots for similarity

make_similarity_plot <- function(filename) {
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
  n<-c(174,103,169)/(256*(256-1)/2)
  b<-c("data", "data", "data")
  a<-c("V1", "V2", "V3")
  f<-c("tet", "tet", "tet")
  df<-data.frame(f, b, a, n)
  colnames(df)<-colnames(data2)
  data2<-rbind(data2, df)
  data2$wat2<-data2$wat
  data2$wat2[data2$rowname=="tet"]<-"tet"
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
  data3$value[data3$rowname=="tet"]<-n
  data3$sd[data3$rowname=="tet"]<-0
  data3$ymax[data3$rowname=="tet"]<-log10(n)
  data3$ymin[data3$rowname=="tet"]<-log10(n)
  
  ggplot(data3, aes(wat, log10(value), fill=wat2))+
    geom_pointrange(data=data3, aes(wat, log10(value), ymin=ymin, ymax=ymax),position = position_jitter(width = 0.4, height = 0), size=1, shape=21, color="#052F5F")+
    theme_light()+
    theme(panel.grid.major.x = element_blank())+
    facet_wrap(~variable)+
    scale_fill_manual(values = c("#F9C80E", "#052F5F", "#CC0014" ), labels=c("YF\nspecific", "random"), name=" name")+
    theme(strip.background = element_rect(fill  = "grey95", colour="grey80"))+
    theme(strip.text = element_text(colour="black", face="bold", size=17))+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=22, color="black"))+
    #theme(axis.text.x = element_blank())+
    theme(legend.position = "none")+
    theme(axis.text.y = element_text(size=22, color="black"))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    #ylab(label=expression(atop(paste(log[10]," normalized number"), "of similar CDR3aa sequencies")))+
    ylab(label="Normalized number of\nsimilar CDR3aa sequences\n")+
    theme(axis.title.y = element_text(size=26))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_y_continuous(breaks=c(-6, -5, -4, -3, -2), labels=c(scientific_10(1e-6), scientific_10(1e-5),scientific_10(1e-4),scientific_10(1e-3), scientific_10(1e-2)))+
  # ggtitle(plotname)+
  #   theme(plot.title = element_text(hjust = 0.5, size=20))+
    scale_x_discrete(labels=c("YF\nspecific", "random"))
}

make_similarity_plot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/similarity_edger.rda")
make_similarity_plot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/similarity_max.rda")


#color palette two #440154FF or #482677FF and #29AF7FFF or #3CBB75FF    "#3D9637"


#Heatmap
make_heatmap <- function(data_heat_csv) {
  heatmap_gr<-read.csv(data_heat_csv)
  sht<-as.matrix(heatmap_gr[, 2:7])
  rownames(sht)<-heatmap_gr[,1]
  row.names(sht)<-gsub("_0_F1","",row.names(sht))
  row.names(sht)<-d_rename[row.names(sht)]
  colnames(sht)<-gsub("_0_F1","",colnames(sht))
  colnames(sht)<-d_rename[colnames(sht)]
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
    theme(plot.margin = unit(c(0,1,1,1), "cm"))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    theme(legend.text = element_text(colour="black", size = 24,  hjust = 1))+
    theme(legend.title=element_text(colour="black", size = 24, face = "bold"))+
    theme(legend.text.align=0)+
    theme(legend.key = element_rect(size = 22),legend.key.size = unit(2.5, 'lines'))+
    #for edgeR and Max: scale_fill_viridis( na.value = "grey60",breaks=c(0, 1e-5, 2e-5, 3e-5), labels=c(" 0", scientific_10(1e-5), scientific_10(2e-5),scientific_10(3e-5)),name="Normalized\nsharing\n ", direction = 1)
    #for top: scale_fill_viridis( na.value = "grey60",breaks=c(4e-6, 8e-6, 1.2e-5, 1.6e-5), labels=c(scientific_10(4e-6), scientific_10(8e-6), scientific_10(1.2e-5),scientific_10(1.6e-5)),name="Normalized\nsharing\n ", direction = 1)
    #for total:
    scale_fill_viridis( na.value = "grey60",breaks=c(4.25e-8, 4.75e-8, 5.25e-8), labels=c(scientific_10(4.25e-8), scientific_10(4.75e-8), scientific_10(5.25e-8) ),name="Normalized\nsharing\n ", direction = 1)
  library(gridExtra)
  library(grid)
  pp1<-grid.arrange(px, p, ncol=1, heights=c(2,6))
  pp1
}
#hm<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/edger_sharing.csv")
make_heatmap("C:/Users/Asya/Dropbox/YF_gr/Fig_data/Max_sharing.csv")
make_heatmap("C:/Users/Asya/Dropbox/YF_gr/Fig_data/total_prevac_sharing.csv")
make_heatmap("C:/Users/Asya/Dropbox/YF_gr/Fig_data/top_prevac_sharing.csv")
#Scatterplot

make_scatter <- function(DR_data_csv) {
  plot<-read.csv(DR_data_csv)
  plot<-arrange(plot, tetramer)
  
  ggplot(plot[plot$tetramer==F, ], mapping=aes(log10(DRhi.Read.proportion), log10(CD8.Read.proportion), col=changed))+
    geom_point(shape=16, size=4, alpha=0.9)+ 
  geom_point(data=plot[plot$tetramer==T, ], mapping=aes(log10(DRhi.Read.proportion), log10(CD8.Read.proportion), fill=changed), col="#E20D2D", shape=21, size=4, alpha=0.9, stroke=2)+
    theme_light()+
    theme(panel.grid.minor.y = element_blank())+
    theme(panel.grid.minor.x = element_blank())+
    scale_color_manual(values = c("#052F5F", "#F9C80E"))+
    scale_fill_manual(values = c("#052F5F", "#F9C80E"), labels=c("not changed\nclones\n ", "changed clones"), name="Legend\nname")+
    theme(legend.title = element_text(size=24, face="bold"))+
    #theme(axis.title.x = element_blank())+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size=22, color="black"))+
    theme(axis.text.y = element_text(size=22, color="black"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2, 'lines'))+
    ylab(label=expression(paste("Clone frequency ", "in bulk CD8"^"+")))+
    xlab(label=expression(paste("Clone frequency ","in CD8"^"+","CD38"^"+","HLADR"^"+")))+
    scale_x_continuous(breaks=c(-5, -4, -3,  -2), labels=c(scientific_10(1e-5), scientific_10(1e-4),
                                                                        scientific_10(1e-3), scientific_10(1e-2)))+
    scale_y_continuous(breaks=c(-5, -4, -3, -2),labels=c(scientific_10(1e-5), scientific_10(1e-4),scientific_10(1e-3),scientific_10(1e-2)))+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=23))+
    theme(axis.title.x = element_text(size=23))+
      geom_abline(slope = 1, intercept = 0, size=1.3)+
    theme(plot.margin = unit(c(1,2,1,1), "cm"))
   
}
plot<-read.csv("C:/Users/Asya/Dropbox/YF_gr/Fig_data/DR_small.csv")
make_scatter("C:/Users/Asya/Dropbox/YF_gr/Fig_data/DR_small.csv")

#make histogram plot
make_histplot <- function(histdatacsv2) {
  histplot<-read.csv2(histdatacsv2)
  histplot<-select(histplot, 2:4)
  histplot2<-data.frame(mismatches=c(histplot$mismatches, histplot$mismatches), Freq=c(histplot$YF_LLW.Freq, histplot$CMV_NLV.Freq),antigen=c(rep("YF",times=21),rep("CMV",times=21)),stringsAsFactors = F)
  histplot2$Freq[histplot2$antigen=="YF"]<-histplot2$Freq[histplot2$antigen=="YF"]/sum(histplot2$Freq[histplot2$antigen=="YF"])
  histplot2$Freq[histplot2$antigen=="CMV"]<-histplot2$Freq[histplot2$antigen=="CMV"]/sum(histplot2$Freq[histplot2$antigen=="CMV"])
  ggplot(histplot2, aes(mismatches, Freq, fill=antigen))+geom_histogram(stat="identity", position = "dodge", color="black")+
    theme_light()+
    theme(panel.grid.major.x = element_blank())+
    theme(panel.grid.minor.x = element_blank())+
    scale_fill_manual(values = c("#052F5F", "#F9C80E"), name="Epitope\nsource")+
    theme(legend.title = element_text(size=24, face="bold"))+
    #theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2, 'lines'))+
    ylab(label = "Fraction of clonotypes")+
    xlab(label="Mismatches to the nearest neighbour")+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=26))+
    theme(axis.title.x = element_text(size=26))+
    scale_y_continuous(expand=c(0,0), limits=c(0, 0.31), breaks=c(0,0.1, 0.2, 0.3), labels=c(0,0.1, 0.2, 0.3))+
    scale_x_continuous( breaks = c(0:12), limits=c(0,13))+
    ggtitle("A.")+
    theme(plot.title = element_text(size=30, face="bold"))+
    theme(legend.justification=c(1,1), legend.position=c(1,1))+
    theme(legend.box.background = element_rect(colour = "black", size=2))
}
make_histplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/hist_data.csv")
#Roc for twins
make_roc_tw <- function(roc_data) {
  roc<-read.csv(roc_data)
  roc$donor<-d_rename[roc$donor]

  ggplot(roc, aes(1-specificity, sensitivity,col=donor))+
    geom_line(size=1.5, aes(color=donor))+
    theme_light()+
    scale_x_continuous( limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
    scale_y_continuous(expand=c(0,0), limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    xlab(label = "1-Specificity")+
    theme(axis.title.x = element_text(size=26))+
    ylab("Sensitivity")+
    scale_color_manual(values = c(my_six_colors))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 15),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    theme(axis.title.y = element_text(size=26))+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    geom_abline(slope = 1, intercept = 0, size=1)+
    ggtitle("C.")+
    theme(plot.title = element_text(size=30, face="bold"))+
    theme(legend.justification=c(1,0), legend.position=c(1,0))+
    theme(legend.box.background = element_rect(colour = "black", size=2))
    
}
make_roc_tw("C:/Users/Asya/Dropbox/YF_gr/Fig_data/leave_one_out_roc.csv")
#Roc for tetramer

make_roc_tet <- function(roc_data) {
  roc_tet<-read.csv(roc_data)
  ggplot(roc_tet, aes(1-specificity, sensitivity))+
    geom_line(size=1.5, col="#052F5F")+
    theme_light()+
    scale_x_continuous(limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
    scale_y_continuous(expand=c(0,0), limits = c(0,1), labels=c(0, 0.25, 0.5, 0.75, 1))+
    theme(panel.grid.minor.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    xlab(label = "1-Specificity")+
    theme(axis.title.x = element_text(size=26))+
    ylab("Sensitivity")+
    theme(axis.title.y = element_text(size=26))+
    theme(axis.text.x = element_text(size=24, color="black"))+
    theme(axis.text.y = element_text(size=24, color="black"))+
    geom_abline(slope = 1, intercept = 0, size=1)+
    ggtitle("B.")+
    theme(plot.title = element_text(size=30, face="bold"))
    
}
make_roc_tet("C:/Users/Asya/Dropbox/YF_gr/Fig_data/tet_roc.csv")
  
#arrange hist and rocs save as 7 21
pl1<-make_histplot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/hist_data.csv")
pl2<-make_roc_tet("C:/Users/Asya/Dropbox/YF_gr/Fig_data/tet_roc.csv")
pl3<-make_roc_tw("C:/Users/Asya/Dropbox/YF_gr/Fig_data/leave_one_out_roc.csv")
grid.arrange(pl1, pl2,pl3, ncol=3)

#Plots about sharing

make_sharing_plot <- function(sharing, total) {

  shVJ_tbl<-read.csv(sharing)
  shVJ_tbl<-shVJ_tbl[shVJ_tbl$Uniq==T, ]
  #shVJ_tbl$Var1<-d_rename[shVJ_tbl$Var1]
  #shVJ_tbl$Var2<-d_rename[shVJ_tbl$Var2]
  shVJ_tbl$YFspec<-"YF specific"
  shVJ_tbl$dn<-paste0(shVJ_tbl$Var1, "/", shVJ_tbl$Var2)
  shVJ_tbl$Freq[shVJ_tbl$X==18]<-NA
  tot<-read.csv(total)
  tot<-tot[tot$Uniq==T, ]
  tot$Var1<-gsub("_0_F1","",tot$Var1)
  tot$Var2<-gsub("_0_F1","",tot$Var2)
  #tot$Var1<-d_rename[tot$Var1]
  #tot$Var2<-d_rename[tot$Var2]
  tot$YFspec<-"total"
  tot$dn<-paste0(tot$Var1, "/", tot$Var2)
  shar<-rbind(shVJ_tbl, tot)
  shar$YFspec<-factor(x = shar$YFspec, levels=c("YF specific", "total"), ordered=T)
  
  ggplot(shar, aes(YFspec, log10(Freq), fill=twins))+
    geom_point(position = position_jitter(width = 0.3, height = 0), size=4.1, shape=21, color="#032449")+
    theme_light()+
    theme(panel.grid.major.x = element_blank())+
    theme(panel.grid.minor.y = element_blank())+
    scale_fill_manual(values = c("#052F5F", "#CC0014"), labels=c("unrelated\ndonors", "twin pairs"), name="")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=22, color="black"))+
    #theme(legend.position = "none")+
    theme(axis.text.y = element_text(size=22, color="black"))+
    theme(legend.title = element_text(size=24, face="bold"))+
    theme(legend.key = element_rect(size = 20),legend.key.size = unit(2.5, 'lines'))+
    theme(legend.text = element_text(size=22, color="black"))+
    ylab(label="Normalized number of\nshared TCR3aa sequences\n")+
    theme(axis.title.y = element_text(size=26))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_y_continuous(breaks=c(-7, -6, -5), labels=c(scientific_10(1e-7), scientific_10(1e-6),scientific_10(1e-5)))+
    scale_x_discrete(labels=c("YF\nspecific", "total"))
    
    # ggtitle(plotname)+
    #   theme(plot.title = element_text(hjust = 0.5, size=20))+
    
}
make_sharing_plot("C:/Users/Asya/Dropbox/YF_gr/Fig_data/shar_Max.csv", "C:/Users/Asya/Dropbox/YF_gr/Fig_data/shar_total.csv")
