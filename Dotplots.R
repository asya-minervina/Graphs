library(ggplot2)
library(data.table)
library(dplyr) 

#multiple plots
require(gridExtra)
plot1 <- qplot(1)
plot2 <- qplot(1)
grid.arrange(plot1, plot2, ncol=2)

#Very sexy heated scatterplot
ppp<-densCols(paper_list$res5$shape1,
              log10(paper_list$res5$rbig), 
              colramp = colorRampPalette(c("#440154FF", #"#39568CFF", "#1F968BFF", "#73D055FF","#481567FF",
                                            "#33638DFF", #"#20A387FF", "#95D840FF","#482677FF", "#2D708EFF",
                                            "#29AF7FFF", #"#B8DE29FF","#453781FF", "#287D8EFF", "#3CBB75FF",
                                           "#DCE319FF")))#"#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF")))
p3<-ggplot(paper_list$res5, aes(paper_list$res5$shape1,log10(paper_list$res5$rbig)))+
  geom_point(alpha=0.8, size=2,  position = "jitter", color=ppp)+
  theme_light()+
  #geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="#440154FF", size=3)+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="red",  size=10, shape=21, stroke = 1.5)

#Just points
ggplot(paper_list$res5, aes(paper_list$res5$shape1,log10(paper_list$res5$rbig)))+
  geom_point(alpha=0.2, size=2,  position = "jitter", color="black")+
  theme_light()+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="black", size=3)+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="red",  size=10, shape=21, stroke = 1.5)

#hexplot blues
ggplot(paper_list$res5, aes(paper_list$res5$shape1,log10(paper_list$res5$rbig)))+
  geom_hex(bins=100)+
  theme_light()+
  #geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="black", size=2)+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="red",  size=10, shape=21, stroke = 1.5)

#hexplot magma
p4<-ggplot(paper_list$res5, aes(paper_list$res5$shape1,log10(paper_list$res5$rbig)))+
  geom_hex(bins=100)+
  scale_fill_viridis(option="magma")+
  theme_light()+
  #geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="black", size=2)+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="red",  size=10, shape=21, stroke = 1.5)

#hexplot viridis
ggplot(paper_list$res5, aes(paper_list$res5$shape1,log10(paper_list$res5$rbig)))+
  geom_hex(bins=100)+
  scale_fill_viridis()+
  theme_light()+
  #geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="black", size=2)+
  geom_point(data=paper_list$res5[paper_list$res5$target==T, ],aes(paper_list$res5$shape1[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="red",  size=10, shape=21, stroke = 1.5)


#440154FF
#39568CFF
#1F968BFF
#73D055FF
#481567FF
#33638DFF
#20A387FF
#95D840FF
#482677FF
#2D708EFF
#29AF7FFF
#B8DE29FF
#453781FF
#287D8EFF
#3CBB75FF
#DCE319FF
#404788FF
#238A8DFF
#55C667FF
#FDE725FF



plot_fancy1<-function(tmp2){
  ppp<-densCols(tmp2$donors,
                log10(tmp2$rbig), 
                colramp = colorRampPalette(c("#440154FF", "#39568CFF","#1F968BFF", "#73D055FF",
                                             "#95D840FF","#B8DE29FF", "#DCE319FF", 
                                           "#FDE725FF")), nbin = c(1000))
  
ggplot(tmp2, aes(tmp2$donors,log10(tmp2$rbig)))+
    geom_point(alpha=0.8, size=2, color=ppp)+
    theme_light()+
  scale_y_continuous(expand=c(0,0), limits=c(0, 5.5))+
  scale_x_continuous(expand=c(0,0), limits=c(0, 17))+
  labs(x=expression("Number of donors"),y=expression(paste(log[10]," (",italic("in silico"), "rearrangements)")))+
  theme(axis.title = element_text(size=18, color="black", face="bold"))+  
  #geom_point(data=tmp2[paper_list$res5$target==T, ],aes(paper_list$res5$donors[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="#440154FF", size=3)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="red", size=10, shape=21, stroke = 1.5)+
  geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="red", size=2, shape=21)
}

plot_fancy2<-function(tmp2){
  ppp<-densCols(tmp2$donors,
                log10(tmp2$rbig), 
                colramp = colorRampPalette(c("#440154FF", "#39568CFF","#1F968BFF",                                              "#73D055FF",
                                             "#95D840FF","#B8DE29FF", "#DCE319FF", 
                                             "#FDE725FF")), nbin = c(1000))
  
  
  ggplot(tmp2, aes(log10(tmp2$ML),log10(tmp2$P_post)))+
    geom_point(alpha=0.8, size=2.5, color=ppp)+
    theme_light()+
    labs(x=expression(paste(log[10], " ", P[data])),y=expression(paste(log[10], " ", P[post])))+
    theme(axis.title = element_text(size=18, color="black", face="bold"))+ 
    geom_abline(intercept=0, slope=1)+
    scale_y_continuous( limits=c(-8, -2))+
    scale_x_continuous(limits=c(-8,-2))+
    ggtitle("V7-6, J2-1, CMV")+
    theme(plot.title = element_text(hjust = 0.5))+
    #geom_point(data=tmp2[paper_list$res5$target==T, ],aes(paper_list$res5$donors[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="#440154FF", size=3)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(log10(tmp2$ML[tmp2$target]),log10(tmp2$P_post)[tmp2$target]), color="red", size=10, shape=21, stroke = 1.5)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(log10(tmp2$ML[tmp2$target]),log10(tmp2$P_post)[tmp2$target]), color="red", size=2.5, shape=21)
}

plot_fancy1_hex<-function(tmp2) {
  ggplot(tmp2, aes(tmp2$donors,log10(tmp2$rbig)))+
  geom_hex(bins=100)+
  scale_fill_viridis()+
  theme_light()+
  #geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$shape1[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="black", size=2)+
  geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="red",  size=10, shape=21, stroke = 1.5)
}


 