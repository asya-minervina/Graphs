library(ggplot2)
library(gridExtra)
library(grid)

#Plot fancy1
#ff3333

plot_fancy1<-function(tmp2, x, y){
  ppp<-densCols(tmp2$donors,
                log10(tmp2$rbig), 
                colramp = colorRampPalette(c("#440154FF", "#39568CFF","#1F968BFF", "#73D055FF",
                                             "#95D840FF","#B8DE29FF", "#DCE319FF", 
                                             "#FDE725FF")), nbin = c(200))
  
  ggplot(tmp2, aes(tmp2$donors,log10(tmp2$rbig)))+
    geom_point(alpha=0.8, size=2.5, color=ppp)+
    theme_light()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 5.5))+
    scale_x_continuous(expand=c(0,0), limits=c(0, x))+
    theme(panel.grid.minor=element_blank())+
    ggtitle(y)+
    theme(axis.text= element_text(size=15))+
    theme(plot.title = element_text(size=18, color="black", face="bold", hjust = 0.5))+
    labs(x=expression("Number of donors"),y=expression(paste(log[10]," (",italic("in silico "), "rearrangements)")))+
    theme(axis.title = element_text(size=17, color="black", face="bold"))+  
    geom_point(data=tmp2[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01, ],aes(tmp2$donors[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01],log10(tmp2$rbig)[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01]), color="white", fill="red2", size=3, shape=21, alpha=0.9, stroke=0.3)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="red1", size=12, shape=21, stroke = 1.8)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="#440154FF", fill="red1", size=3.8, shape=21, stroke=1.2)
}

# Make plots
p1<-plot_fancy1(paper_VJ$res5, 165, y="V5-1, J2-6, CMV")
p2<-plot_fancy1(paper_VJ$res7, 45, y="V7-6, J1-4, CMV")
p3<-plot_fancy1(paper_VJ$res12, 155, y="V12-3,-4, J1-2, CMV")
p4<-plot_fancy1(paper_VJ$resGAD, 17, y="V5-1, J1-1, T1D")

#fancy 2 plot
plot_fancy2<-function(tmp2, x1, x2, y1=-8){
  ppp<-densCols(tmp2$donors,
                log10(tmp2$rbig), 
                colramp = colorRampPalette(c("#440154FF", "#39568CFF","#1F968BFF","#73D055FF",
                                             "#95D840FF","#B8DE29FF", "#DCE319FF", 
                                             "#FDE725FF")), nbin = c(200))
  
  
  ggplot(tmp2, aes(log10(tmp2$ML),log10(tmp2$P_post)))+
    geom_point(alpha=0.8, size=2.5, color=ppp)+
    theme_light()+
    labs(x=expression(paste(log[10], " ", P[data])),y=expression(paste(log[10], " ", P[post])))+
    theme(axis.title = element_text(size=17, color="black", face="bold"))+ 
    geom_abline(intercept=0, slope=1)+
    scale_y_continuous( limits=c(y1, -2))+
    scale_x_continuous(limits=c(x1,x2))+
    theme(panel.grid.minor=element_blank())+
    theme(axis.text= element_text(size=15))+       
    geom_point(data=tmp2[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01, ],aes(log10(tmp2$ML[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01]),log10(tmp2$P_post)[tmp2$rbig>0&p.adjust(tmp2$pval_post,method="holm")<0.01]), color="white", fill="red2", size=3, shape=21, alpha=0.9, stroke=0.3)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(log10(tmp2$ML[tmp2$target]),log10(tmp2$P_post)[tmp2$target]), color="red1", size=12, shape=21, stroke = 1.8)+
    geom_point(data=tmp2[tmp2$target==T, ],aes(log10(tmp2$ML[tmp2$target]),log10(tmp2$P_post)[tmp2$target]), color="#440154FF", fill="red1", size=3.5, shape=21, stroke=1.2)
    
}

#make plots2
p5<-plot_fancy2(paper_VJ$res5, -5,-2)
p6<-plot_fancy2(paper_VJ$res7, -4.5, -2.5, y1=-7)
p7<-plot_fancy2(paper_VJ$res12, -5.5, -2.5)
p8<-plot_fancy2(paper_VJ$resGAD, -4.5, -2.5, y1=-8)

#arrange plots
title1<-textGrob("A.", x=0, just="left", gp=gpar(fontface="bold", fontsize=22))
title2<-textGrob("B.", x=0, just="left", gp=gpar(fontface="bold", fontsize=22))
pp1<-grid.arrange(p1, p2, p3, p4, top=title1, ncol=4)
pp2<-grid.arrange(p5,p6,p7,p8, top=title2, ncol=4)
grid.arrange(pp1, pp2, ncol=1)

#http://pdf2tiff.com/
  
