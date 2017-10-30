Vlist<-c("TCRBV10-01","TCRBV10-02","TCRBV10-03","TCRBV11-01","TCRBV11-02","TCRBV11-03","TCRBV12-03","TCRBV12-04","TCRBV12-05","TCRBV13-01","TCRBV14-01","TCRBV15-01","TCRBV16-01","TCRBV18-01","TCRBV19-01","TCRBV02-01","TCRBV20-01","TCRBV21-01","TCRBV23-01","TCRBV24-01","TCRBV25-01","TCRBV27-01","TCRBV28-01","TCRBV29-01","TCRBV03-01","TCRBV30-01","TCRBV04-01","TCRBV04-02","TCRBV04-03","TCRBV05-01","TCRBV05-04","TCRBV05-05","TCRBV05-06","TCRBV05-08","TCRBV06-01","TCRBV06-02","TCRBV06-03","TCRBV06-04","TCRBV06-05","TCRBV06-06","TCRBV06-07","TCRBV07-01","TCRBV07-02","TCRBV07-03","TCRBV07-04","TCRBV07-06","TCRBV07-07","TCRBV07-08","TCRBV07-09","TCRBV09-01")
Jlist<-c("TCRBJ01-01","TCRBJ01-02","TCRBJ01-03","TCRBJ01-04","TCRBJ01-05","TCRBJ01-06","TCRBJ02-01","TCRBJ02-02","TCRBJ02-03","TCRBJ02-04","TCRBJ02-05","TCRBJ02-06","TCRBJ02-07")

preproc_script<-function(){command<-"mkdir TCRBV07-06_TCRBJ01-04\nfor fn in *.tsv\ndo\ngrep \"TCRBV07-06.*TCRBJ01-04\" ${fn} > TCRBV07-06_TCRBJ01-04/${fn}_TCRBV07-06_TCRBJ01-04.tsv;\npaste <(cut -f2 TCRBV07-06_TCRBJ01-04/${fn}_TCRBV07-06_TCRBJ01-04.tsv) <(cut -f8,22 TCRBV07-06_TCRBJ01-04/${fn}_TCRBV07-06_TCRBJ01-04.tsv)|sed 's/TCR/TR/g' > TCRBV07-06_TCRBJ01-04/${fn}_TCRBV07-06_TCRBJ01-04_p.tsv;\nrm TCRBV07-06_TCRBJ01-04/${fn}_TCRBV07-06_TCRBJ01-04.tsv;\ndone\n\n"
code=""
for (V in Vlist)
  for(J in Jlist)
  {
    code<-paste0(code,gsub("TCRBJ01-04",J,gsub("TCRBV07-06",V,command)),collapse="")
  }
write(code,"preproc.sh")
}

rob_full_names<-c("nucleotide", "aminoAcid", "count (templates/reads)", "frequencyCount (%)", 
                  "cdr3Length", "vMaxResolved", "vFamilyName", "vGeneName", "vGeneAllele", 
                  "vFamilyTies", "vGeneNameTies", "vGeneAlleleTies", "dMaxResolved", 
                  "dFamilyName", "dGeneName", "dGeneAllele", "dFamilyTies", "dGeneNameTies", 
                  "dGeneAlleleTies", "jMaxResolved", "jFamilyName", "jGeneName", 
                  "jGeneAllele", "jFamilyTies", "jGeneNameTies", "jGeneAlleleTies", 
                  "vDeletion", "n1Insertion", "d5Deletion", "d3Deletion", "n2Insertion", 
                  "jDeletion", "vIndex", "n1Index", "dIndex", "n2Index", "jIndex", 
                  "estimatedNumberGenomes", "sequenceStatus", "cloneResolved", 
                  "vOrphon", "dOrphon", "jOrphon", "vFunction", "dFunction", "jFunction", 
                  "fractionNucleated", "vAlignLength", "vAlignSubstitutionCount", 
                  "vAlignSubstitutionIndexes", "vAlignSubstitutionGeneThreePrimeIndexes", 
                  "vSeqWithMutations")

read_rob_no_header<-function(folder){
  DTlist<-lapply(list.files(folder,full.names = T),fread,header=F,col.names=rob_full_names)
  names(DTlist)<-list.files(folder,full.names = F)
  lapply(DTlist,setkey)
  DTlist<-lapply(DTlist,unique)
  DTlist<-lapply(DTlist,function(x)x[aminoAcid!="",,])
  DTlist
}
#> Robl12_01_02<-lapply(Robl12,function(x)x[jGeneName=="TCRBJ01-02",unique(aminoAcid),])


read_rob<-function(folder){
  DTlist<-lapply(list.files(folder,full.names = T),fread)
  names(DTlist)<-list.files(folder,full.names = F)
  lapply(DTlist,setkey)
  DTlist<-lapply(DTlist,unique)
  DTlist<-lapply(DTlist,function(x)x[aminoAcid!="",,])
  DTlist
}

merge_dt_list_rob<-function(DTlist){
  #get all in one dt(two colummns - CDR3 and colname)
  #clnames<-c("CDR3.nucleotide.sequence",colname)
  DTlist<-lapply(DTlist,function(x)x[,.(aminoAcid),])
  DTlist<-lapply(DTlist,function(x)x[,p:=1,])
  #change names
  lapply(names(DTlist),function(x)setnames(DTlist[[x]],"p",paste0(x,".","p")))
  DTlistm<-rbindlist(DTlist,fill=T)
  print("Binded")
  DTlistm<-DTlistm[,lapply(.SD, sum,na.rm=T),aminoAcid]
  #print(names(DTlistm)[grepl(colname, names(DTlistm))])
  print("Merged")
  DTlistm[, Sum := Reduce(`+`, .SD), .SDcols=grep("p", names(DTlistm))][]
  print("Sums done")
  DTlistm
  #DTlistm[, Mean :=Sum/(length(grep(colname, names(DTlistm)))),][]
  #print("Means done")
}

#> Robl12_01_02<-lapply(Robl12,function(x)x[jGeneName=="TCRBJ01-02",.(aminoAcid),][!duplicated(aminoAcid),,][aminoAcid!="",,])
# shrep<-as.data.frame(merge_dt_list_rob(lapply(Robl12_01_02,FUN = function(x)x))[Sum>9,,]);save(shrep,file=("specialTCRBV12_TCRBJ01-02.rda"))


make_tonn_rda<-function(DTlist,thres=9,index=""){
k<-0
res<-list()
for (V in Vlist)
  for(J in Jlist)
  {
    k<-k+1
    print(k)
    print(format(Sys.time(), "%a %b %d %X %Y"))
    print(V)
    print(J)
    shrep<-as.data.frame(merge_dt_list_rob(lapply(DTlist,FUN = function(x)x[vGeneName==V&jGeneName==J,,]))[Sum>thres,,])
    names(shrep)[1]<-"CDR3.amino.acid.sequence"
    print(c(nrow(shrep),max(shrep$Sum)))
    res[[k]]<-c(nrow(shrep),max(shrep$Sum))
    save(shrep,file=paste0(index,V,"_",J,".rda",collapse = ""))
  }
res
}

resolve_rob<-function(robdt){# this function adds TCRBV20-1 and TCRBV12-03 into the dataset
  robdt[vMaxResolved=="TCRBV20",vGeneName:="TCRBV20-01",]
  robdt[vMaxResolved=="TCRBV12",vGeneName:="TCRBV12-03",]
  robdt[vMaxResolved=="TCRBV06",vGeneName:="TCRBV06-02",]
  robdt
 # TCRBV06-02
}

annotate_narco<-function(resRob){
  load("narco_names.rda")
  resRob$narco=rowSums(resRob[,narco_patients]!=0)
  resRob$control=rowSums(resRob[,narco_healthy]!=0)
  resRob
}

plot_narco<-function(resRob,title=""){
  plot(resRob$control+resRob$narco,log10(resRob$rbig+1),col=((resRob$control==0)+1),main=title,ylab="log10 (# assemblies+1)",xlab="number of donors with this TCR");
  abline(v=2.5)
  abline(h=2)
}

annotate_narco_2<-function(resRob,sample_inds,donor_inds,donor_status)#general annotation function. Picks ill_inds, healthy_inds and donor_inds - these are patient names. 
{
  #tst<-list();for (i in 1:27){print(i);tst[[i]]<-rowSums(resRob[,(which(T1D_cohorts$bin==i)+2),drop=F])}
  tst<-list();
  for (i in 1:length(unique(donor_inds)))
  {
    tst[[i]]<-rowSums(resRob[,sample_inds][,which(donor_inds==unique(donor_inds)[i]),drop=F]!=0)!=0
  #print("tst")
  }
  resRob$narco<-rowSums(do.call(cbind,tst[donor_status=="narco"]))
  resRob$control<-rowSums(do.call(cbind,tst[donor_status=="hd"]))
  resRob
}

#V="TCRBV05-01"
#J="TCRBJ02-07"
#shrep<-as.data.frame(merge_dt_list_rob(lapply(Robl,FUN = function(x)x[vGeneName==V&jGeneName==J,,]))[Sum>9,,])
#names(shrep)[1]<-"CDR3.amino.acid.sequence"
#print(c(nrow(shrep),max(shrep$Sum)))
#res[[k]]<-c(nrow(shrep),max(shrep$Sum))
#save(shrep,file=paste0(V,"_",J,"no_unique.rda",collapse = ""))
getp3<-function(i,p)dbeta(p,shape1=sum(resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==1),shape2=sum(N9_1_2_7vec[CMV_ind1-2][resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==0])-sum(resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==1))
getp2<-function(i,p)p^sum(resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==1)*prod((1-p)^N9_1_2_7vec[CMV_ind1-2][resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==0])
getp<-function(i,p)prod(1-exp(-N9_1_2_7vec[CMV_ind1-2][resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==1]*p))*prod(exp(-N9_1_2_7vec[CMV_ind1-2][resRob[i,-c(1,2,c(789:ncol(resRob)))][CMV_ind1-2]==0]*p))
plot_i<-function(i){
  plot(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp,i=i)[-1]*diff(seq(1e-6,1e-4,by = 1e-6))))
  points(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp2,i=i)[-100]*diff(seq(1e-6,1e-4,by = 1e-6))),col="red")
   points(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp3,i=i)[-100]*diff(seq(1e-6,1e-4,by = 1e-6))),col="green")}

plot_i3<-function(i){
  plot(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp3,i=i)[-1]*diff(seq(1e-6,1e-4,by = 1e-6))))
  #points(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp2,i=1)[-100]*diff(seq(1e-6,1e-4,by = 1e-6))),col="red")
  #points(seq(1e-6,1e-4,by = 1e-6)[-100],prop.table(sapply(seq(1e-6,1e-4,by = 1e-6),getp3,i=1)[-100]*diff(seq(1e-6,1e-4,by = 1e-6))),col="green")}
}

add_shape<-function(resRob,indvec,sizevec){
  print(nrow(resRob))
  resRob<-as.data.frame(resRob)
  #if (nrow(resRob==0)){return(NULL);break}
  shape1=numeric(nrow(resRob))
  shape2=numeric(nrow(resRob))
  arrInt<-resRob[,indvec]
  #present= arrInt==1
 # absent= arrInt==0
  print(sizevec[indvec-2])
  #sizemat=as.numeric(t(matrix(data = sizevec[indvec-2],nrow=length(sizevec[indvec-2]),ncol=nrow(arrInt))))
 
  #for (i in 1:nrow(resRob))
  {
  shape1=rowSums(arrInt)
  #shape2=(as.numeric(sum(as.numeric(sizemat[absent|present])))-as.numeric(shape1)) # aint that bullshit??? why absent?
  shape2=as.numeric(sum(sizevec[indvec-2]))-shape1 # aint this better?
  }
  resRob$shape1=shape1
  resRob$shape2=shape2
  resRob
}

#process
#sizevec_f<-lapply(list.files(pattern = "resTCRBV"),function(f){sapply(Robl,FUN = function(x)x[vGeneName==gsub("res","",unlist(strsplit(f,split="_"))[1])&jGeneName==gsub(".rda","",unlist(strsplit(f,split="_"))[2]),.N,])})
#rob_resl<-list();for (f in list.files(pattern = "resTCRBV")){load(f);rob_resl[[f]]<-resRob;rm(resRob)}
#rob_resls<-lapply(list.files(pattern = "resTCRBV"),function(f){print(f);add_shape(resRob=rob_resl[[f]],indvec = CMV_ind1,sizevec = sizevec_f[[f]])})

plot_rob<-function(resRob){#shape1/(shape1+shape2)
  smoothScatter(log10(qbeta(0.975,resRob$shape1,resRob$shape2)),log10(resRob$rbig))
  points(log10(qbeta(0.975,resRob[rob_int$AA,]$shape1,resRob[rob_int$AA,]$shape2)),log10(resRob[rob_int$AA,]$rbig),pch=19,col="red")
}

plot_robl<-function(resRobl,ind){#shape1/(shape1+shape2)
  plot(log10(qbeta(0.975,resRobl[[ind]]$shape1,resRobl[[ind]]$shape2)),log10(resRobl[[ind]]$rbig+1),main=ind,xlab="P-estimate",ylab="Log10(rearangements)")
  print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  print(rob_int_f)
  points(log10(qbeta(0.975,resRobl[[ind]][rob_int_f$AA,]$shape1,resRobl[[ind]][rob_int_f$AA,]$shape2)),log10(resRobl[[ind]][rob_int_f$AA,]$rbig+1),pch=19,col="red")
}

plot_robl2<-function(resRobl,ind){#shape1/(shape1+shape2)
  smoothScatter(log10(resRobl[[ind]]$shape1/(resRobl[[ind]]$shape1+resRobl[[ind]]$shape2)),log10((resRobl[[ind]]$rbig+1)/(100e6/3)),main=ind,xlab="P-estimate",ylab="Log10(rearangements)",xlim=c(-6,-2))
  #graphics::segments(y0 = ,y1=)
  print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  print(rob_int_f)
  points(log10(resRobl[[ind]][rob_int_f$AA,]$shape1/(resRobl[[ind]][rob_int_f$AA,]$shape1+resRobl[[ind]][rob_int_f$AA,]$shape2)),log10((resRobl[[ind]][rob_int_f$AA,]$rbig+1)/(100e6/3)),pch=19,col="red")
}

plot_robl3<-function(resRobl,ind){#shape1/(shape1+shape2)
  plot((resRobl[[ind]]$shape1),log10((resRobl[[ind]]$rbig+1)/(100e6/3)),main=ind,xlab="P-estimate",ylab="Log10(rearangements)",xlim=c(0,210),ylim=c(-7.5,-3.5))
  #graphics::segments(y0 = ,y1=)
  print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  print(rob_int_f)
  points((resRobl[[ind]][rob_int_f$AA,]$shape1),log10((resRobl[[ind]][rob_int_f$AA,]$rbig+1)/(100e6/3)),pch=19,col="red")
}


values_robl2<-function(resRobl,ind){
  if (!is.null(resRobl[[ind]])){
  rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  resRobl[[ind]]$q=10^resq[[ind]]
  resRobl[[ind]][rob_int_f$AA,,]}
}

plot_rob_s<-function(){
  
  for (i in 1:97)if(!is.null(rob_resls2[[i]])){
    #print(nrow(rob_resls2[[i]]))
   # print(i)
    pdf(file=paste(names(rob_resls2)[i], "2.pdf", sep = ""))
    plot_robl3(rob_resls2,names(rob_resls2)[i])
    dev.off()
  }
}

plot_narco_s<-function(res_narco){
  
  for (i in 1:length(res_narco))if(!is.null(res_narco[[i]])){
    #print(nrow(rob_resls2[[i]]))
    # print(i)
    pdf(file=paste(names(res_narco)[i], "_narco.pdf", sep = ""))
    plot_narco(annotate_narco(res_narco[[i]]),names(res_narco)[i])
    dev.off()
  }
}


get_q<-function(x){
xtmp<-x[x$shape1>0,]
coef(lm(log10(xtmp$shape1/xtmp$shape2)~offset(log10((xtmp$rbig+1)/(100e6/3))),singular.ok = F))
}

get_q_offset<-function(x,total=2e9){
  xtmp<-x[x$shape1>0,]
  coef(lm(log10(xtmp$shape1/xtmp$shape2)~offset(log10((xtmp$rbig+1)/(total/3))),singular.ok = F))
}

get_q_offset_ML<-function(x,total=2e9){
  xtmp<-x#[x$shape1>0,]
  coef(lm(log10(xtmp$ML)~offset(log10((xtmp$rbig+1)/(total/3))),singular.ok = F))
}


plot_rob_err<-function(resRobl,ind){
  plot(log10(resRobl[[ind]]$shape1/(resRobl[[ind]]$shape1+resRobl[[ind]]$shape2)),log10((resRobl[[ind]]$rbig+1)/(100e6/3)),main=ind,xlab="P-estimate",ylab="Log10(rearangements)",xlim=c(-6,-2))
  errdata1<-qbeta(0.025,shape1 = resRobl[[ind]]$shape1,resRobl[[ind]]$shape2)
  errdata2<-qbeta(0.975,shape1 = resRobl[[ind]]$shape1,resRobl[[ind]]$shape2)
  graphics::segments(x0=log10(errdata1),x1=log10(errdata2),y0=log10((resRobl[[ind]]$rbig+1)/(100e6/3)),y1=log10((resRobl[[ind]]$rbig+1)/(100e6/3)))
  errexp1<-qbeta(0.025,shape1 = resRobl[[ind]]$rbig,100e6/3)
  errexp2<-qbeta(0.975,shape1 = resRobl[[ind]]$rbig,100e6/3)
  graphics::segments(y0=log10(errexp1),y1=log10(errexp2),x0=log10(resRobl[[ind]]$shape1/(resRobl[[ind]]$shape1+resRobl[[ind]]$shape2)),x1=log10(resRobl[[ind]]$shape1/(resRobl[[ind]]$shape1+resRobl[[ind]]$shape2)))
  
  #graphics::segments(y0 = ,y1=)
  print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  print(rob_int_f)
  points(log10(resRobl[[ind]][rob_int_f$AA,]$shape1/(resRobl[[ind]][rob_int_f$AA,]$shape1+resRobl[[ind]][rob_int_f$AA,]$shape2)),log10((resRobl[[ind]][rob_int_f$AA,]$rbig+1)/(100e6/3)),pch=19,col="red")
    
}

plot_rob_err_paper<-function(resRob,total=2e9,q_offset=get_q_offset(resRob),...){
  plot(log10(resRob$shape1/(resRob$shape1+resRob$shape2)),log10((resRob$rbig)/(total/3))+q_offset,xlab="P-post estimate, data",ylab="P-post, model",...)
  errdata1<-qbeta(0.025,shape1 = resRob$shape1,resRob$shape2)
  errdata2<-qbeta(0.975,shape1 = resRob$shape1,resRob$shape2)
  graphics::segments(x0=log10(errdata1),x1=log10(errdata2),y0=log10((resRob$rbig)/(total/3))+q_offset,y1=log10((resRob$rbig)/(total/3))+q_offset,col=rgb(0,0,0,0.3))
  errexp1<-qbeta(0.025,shape1 = resRob$rbig,total/3)
  errexp2<-qbeta(0.975,shape1 = resRob$rbig,total/3)
  graphics::segments(y0=log10(errexp1)+q_offset,y1=log10(errexp2)+q_offset,x0=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),x1=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),col=rgb(0,0,0,0.3))
  
  #graphics::segments(y0 = ,y1=)
  #print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  #print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  #rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
 # print(rob_int_f)
  points(log10(resRob[resRob$target,]$shape1/(resRob[resRob$target,]$shape1+resRob[resRob$target,]$shape2)),log10((resRob[resRob$target,]$rbig)/(total/3))+q_offset,pch=19,col="red")
  abline(c(0,1))
}
plot_rob_noerr_paper<-function(resRob,total=2e9,q_offset=get_q_offset(resRob),...){
  plot(log10(resRob$shape1/(resRob$shape1+resRob$shape2)),log10((resRob$rbig)/(total/3))+q_offset,xlab="P-post estimate, data",ylab="P-post, model",...)
  errdata1<-qbeta(0.025,shape1 = resRob$shape1,resRob$shape2)
  errdata2<-qbeta(0.975,shape1 = resRob$shape1,resRob$shape2)
  #graphics::segments(x0=log10(errdata1),x1=log10(errdata2),y0=log10((resRob$rbig)/(total/3))+q_offset,y1=log10((resRob$rbig)/(total/3))+q_offset,col=rgb(0,0,0,0.3))
  errexp1<-qbeta(0.025,shape1 = resRob$rbig,total/3)
  errexp2<-qbeta(0.975,shape1 = resRob$rbig,total/3)
 # graphics::segments(y0=log10(errexp1)+q_offset,y1=log10(errexp2)+q_offset,x0=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),x1=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),col=rgb(0,0,0,0.3))
  
  #graphics::segments(y0 = ,y1=)
  #print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  #print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  #rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  # print(rob_int_f)
  points(log10(resRob[resRob$target,]$shape1/(resRob[resRob$target,]$shape1+resRob[resRob$target,]$shape2)),log10((resRob[resRob$target,]$rbig)/(total/3))+q_offset,pch=19,col="red")
  abline(c(0,1))
}

plot_rob_simple_paper<-function(resRob,total=2e9,q_offset=get_q_offset(resRob),...){
  plot(resRob$shape1,log10((resRob$rbig)/total),xlab="Number of donors",ylab="fraction of in silico rearrangements",...)
  #errdata1<-qbeta(0.025,shape1 = resRob$shape1,resRob$shape2)
  #errdata2<-qbeta(0.975,shape1 = resRob$shape1,resRob$shape2)
  #graphics::segments(x0=log10(errdata1),x1=log10(errdata2),y0=log10((resRob$rbig)/(total/3))+q_offset,y1=log10((resRob$rbig)/(total/3))+q_offset,col=rgb(0,0,0,0.3))
  errexp1<-qbeta(0.025,shape1 = resRob$rbig,total/3)
  errexp2<-qbeta(0.975,shape1 = resRob$rbig,total/3)
  # graphics::segments(y0=log10(errexp1)+q_offset,y1=log10(errexp2)+q_offset,x0=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),x1=log10(resRob$shape1/(resRob$shape1+resRob$shape2)),col=rgb(0,0,0,0.3))
  
  #graphics::segments(y0 = ,y1=)
  #print(gsub("res","",unlist(strsplit(ind,split="_"))[1]))
  #print(gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])))
  #rob_int_f<-rob_int[V==gsub("res","",unlist(strsplit(ind,split="_"))[1])&J==gsub(".rda","",gsub("res","",unlist(strsplit(ind,split="_"))[2])),,]
  # print(rob_int_f)
  points(resRob[resRob$target,]$shape1,log10((resRob[resRob$target,]$rbig)/total),pch=19,col="red")
  
  #points(log10(resRob[resRob$target,]$shape1/(resRob[resRob$target,]$shape1+resRob[resRob$target,]$shape2)),log10((resRob[resRob$target,]$rbig)/(total/3))+q_offset,pch=19,col="red")
  #abline(c(0,1))
}

#add_p_value_conv(paper_list$res5,sizevec=sizevec_f$`resTCRBV05-01_TCRBJ02-06.rda`,indvec=CMVinds)
add_p_value_conv<-function(resRob,total=2e9,q_offset=get_q_offset(resRob),indvec,sizevec){
  p_val<-numeric(nrow(resRob))
  pred_conv<-numeric(nrow(resRob))
  effect_size<-numeric(nrow(resRob))
  for (i in 1:nrow(resRob))
  {
    conv_dist<-(return_convolution(1-exp(-sizevec[indvec-2]*(10^q_offset)*resRob$rbig[i]/(total/3))))
    pred_conv[i]<-which.max(conv_dist)-1
    p_val[i]<-sum(conv_dist[-c(1:(resRob$shape1[i]+1))])
    effect_size[i]<-log10(resRob$shape1[i]/(resRob$shape1[i]+resRob$shape2[i]))-log10((10^q_offset)*resRob$rbig[i]/(total/3))
  }
  resRob$p_val<-p_val
  resRob$pred_conv<-pred_conv
  resRob$effect_size<-effect_size
  resRob
}
add_p_value_conv2<-function(resRob,total=2e9,q_offset=get_q_offset_ML(resRob),indvec,sizevec){
  p_val<-numeric(nrow(resRob))
  pred_conv<-numeric(nrow(resRob))
  effect_size<-numeric(nrow(resRob))
  for (i in 1:nrow(resRob))
  {
    conv_dist<-(return_convolution(1-exp(-sizevec[indvec-2]*(10^q_offset)*resRob$rbig[i]/(total/3))))
    pred_conv[i]<-which.max(conv_dist)-1
    p_val[i]<-sum(conv_dist[-c(1:(resRob$shape1[i]+1))])
    effect_size[i]<-log10(resRob$ML[i])-log10((10^q_offset)*resRob$rbig[i]/(total/3))
  }
  resRob$p_val<-p_val
  resRob$pred_conv<-pred_conv
  resRob$effect_size<-effect_size
  resRob
}

plot_plot_rob_noerr_paper_sign<-function(resRob,total=2e9,q_offset=get_q_offset(resRob),...)
{
  plot_rob_noerr_paper(resRob,cex=0.5,col=rgb(0,0,0,0.3),...=...)  
  points(log10(resRob[p.adjust(resRob$p_val,method = "holm")<0.01,]$shape1/(resRob[p.adjust(resRob$p_val,method = "holm")<0.01,]$shape1+resRob[p.adjust(resRob$p_val,method = "holm")<0.01,]$shape2)),log10((resRob[p.adjust(resRob$p_val,method = "holm")<0.01,]$rbig)/(total/3))+q_offset,col="blue")
}

get_conv_est_dist<-function(rbig){(return_convolution(1-exp(-sizevec_f$`resTCRBV05-01_TCRBJ02-06.rda`[CMVinds-2]*7.07*rbig/(2e9/3))))}
get_conv_est_5<-function(rbig){which.max(return_convolution(1-exp(-sizevec_f$`resTCRBV05-01_TCRBJ02-06.rda`[CMVinds-2]*7.07*rbig/(2e9/3))))-1}
get_p_est_5<-function(rbig,ndonors){sum((return_convolution(1-exp(-sizevec_f$`resTCRBV05-01_TCRBJ02-06.rda`[CMVinds-2]*7.07*rbig/(2e9/3))))[-(1:(ndonors+1))])}

get_P_Data_posterior<-function(x_vec,n_vec,ps=10^-seq(8,0,length.out = 100)){
  #print(ps)
  likelihood=sapply(ps,function(p)prod(c(1-exp(-n_vec[x_vec]*p)),exp(-n_vec[!x_vec]*p)))
  posteriorP=prop.table(diff(ps)*likelihood[-1])
  cs<-cumsum(posteriorP)
  list(probs=posteriorP,cumsum=cs,CI=c(ps[which(cs>0.025)[1]],ps[which(cs>0.975)[1]]),maxP=ps[which.max(likelihood)])
}

add_p_data_intervals<-function(resRob,total=2e9,indvec,sizevec){
  left<-numeric(nrow(resRob))
  right<-numeric(nrow(resRob))
  ML<-numeric(nrow(resRob))
  donors<-numeric(nrow(resRob))
  
  #effect_size<-numeric(nrow(resRob))
  for (i in 1:nrow(resRob))
  {
  #  print (100*i/nrow(resRob))
    x_vec<-(resRob[i,grepl(".p",names(resRob),fixed = T)]!=0)[indvec-2]
    n_vec<-sizevec[indvec-2]
    posterior=get_P_Data_posterior(x_vec = x_vec,n_vec = n_vec)
    left[i]=posterior$CI[1]
    right[i]=posterior$CI[2]
    ML[i]=posterior$maxP
    donors[i]=sum(x_vec)
    #conv_dist<-(return_convolution(1-exp(-sizevec[indvec-2]*(10^q_offset)*resRob$rbig[i]/(total/3))))
    #pred_conv[i]<-which.max(conv_dist)-1
    #p_val[i]<-sum(conv_dist[-c(1:(resRob$shape1[i]+1))])
   # $effect_size[i]<-log10(resRob$shape1[i]/(resRob$shape1[i]+resRob$shape2[i]))-log10((10^q_offset)*resRob$rbig[i]/(total/3))
  }
  resRob$left<-left
  resRob$right<-right
  resRob$ML<-ML
  resRob$donors<-donors
  #resRob$effect_size<-effect_size
  resRob
}

add_p_value_pdata<-function(resRob,q_offset=get_q_offset_ML(resRob),total=2e9,sizevec,indvec){
  resRob$P_post=(resRob$rbig*3/total)*10^q_offset
  ps=10^-seq(8,0,length.out = 100)
  pval<-numeric(nrow(resRob))
  effect_size<-numeric(nrow(resRob))
  
  for (i in 1:nrow(resRob))
  {
    #print (100*i/nrow(resRob))
    x_vec<-(resRob[i,grepl(".p",names(resRob),fixed = T)]!=0)[indvec-2]
    n_vec<-sizevec[indvec-2]
    
    posterior=get_P_Data_posterior(x_vec = x_vec,n_vec = n_vec,ps = ps)
    #left[i]=posterior$CI[1]
    #right[i]=posterior$CI[2]
    #ML[i]=posterior$maxP
    pval[i]=sum(posterior$cumsum[which(ps>resRob$P_post[i])[1] ])
    effect_size[i]<-log10(resRob$ML[i])-log10(resRob$P_post[i])
  }
  resRob$pval_post=pval
  resRob$effect_size<-effect_size
  resRob
}

all_analysis<-function(resRob,total,sizevec,indvec){
  tmp<-add_p_data_intervals(resRob,total =total,sizevec = sizevec, indvec=indvec )
  add_p_value_pdata(tmp,total =total,sizevec = sizevec, indvec=indvec)
}

plot_fancy1<-function(tmp2){
  
  ppp<-densCols(tmp2$donors,
                 log10(tmp2$rbig), 
                 colramp = colorRampPalette(c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF","#481567FF",
                                                "#33638DFF", "#20A387FF", "#95D840FF","#482677FF", "#2D708EFF",
                                                "#29AF7FFF", "#B8DE29FF","#453781FF", "#287D8EFF", "#3CBB75FF",
                                                "#DCE319FF")),nbin = 1000)#"#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF")))
   
     ggplot(tmp2, aes(tmp2$donors,log10(tmp2$rbig)))+
     geom_point(alpha=0.8, size=2, position = "jitter", color=ppp)+
     theme_light()+
     #geom_point(data=tmp2[paper_list$res5$target==T, ],aes(paper_list$res5$donors[paper_list$res5$target],log10(paper_list$res5$rbig)[paper_list$res5$target]), color="#440154FF", size=3)+
     geom_point(data=tmp2[tmp2$target==T, ],aes(tmp2$donors[tmp2$target],log10(tmp2$rbig)[tmp2$target]), color="red", size=10, shape=21, stroke = 1.5)  
}

plot_fancy1_hex<-function(resRob){
  
  ggplot(resRob, aes(resRob$donors,log10(resRob$rbig)))+
    geom_hex()+
    theme_light()+
    xlab(label = "Number of donors")+
    ylab(label = "log10 number of in silico recombinations" )+
    #geom_point(data=resRob[resRob$target==T, ],aes(resRob$donors[resRob$target],log10(resRob$rbig)[resRob$target]), color="#440154FF", size=3)+
    geom_point(data=resRob[resRob$target==T, ],aes(resRob$donors[resRob$target],log10(resRob$rbig)[resRob$target]), color="red", size=10, shape=21, stroke = 1.5)+
    scale_fill_continuous(viridis)
  
}

plot_fancy2<-function(resRob){
  ppp<-densCols(resRob$ML,
                log10(resRob$P_post), 
                colramp = colorRampPalette(c("#440154FF", #"#39568CFF", "#1F968BFF", "#73D055FF","#481567FF",
                                             "#33638DFF", #"#20A387FF", "#95D840FF","#482677FF", "#2D708EFF",
                                             "#29AF7FFF", #"#B8DE29FF","#453781FF", "#287D8EFF", "#3CBB75FF",
                                             "#DCE319FF")))#"#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF")))
  
  ggplot(resRob, aes(resRob$ML,log10(resRob$P_post)))+
    geom_point(alpha=0.8, size=2, position = "jitter", color=ppp)+
    theme_light()+
    xlab(label = "Pdata")+
    ylab(label = "P_post" )+
    #geom_point(data=resRob[resRob$target==T, ],aes(resRob$ML[resRob$target],log10(resRob$P_post)[resRob$target]), color="#440154FF", size=3)+
    geom_point(data=resRob[resRob$target==T, ],aes(resRob$ML[resRob$target],log10(resRob$P_post)[resRob$target]), color="red", size=10, shape=21, stroke = 1.5)
  
}


#process rob
#rob_resl<-list();for (f in list.files(pattern = "resTCRBV")){load(f);rob_resl[[f]]<-resRob;rm(resRob)} deserialize
# get size vecs
#sizevec_f<-lapply(list.files(pattern = "resTCRBV"),function(f){sapply(Robl,FUN = function(x)x[vGeneName==gsub("res","",unlist(strsplit(f,split="_"))[1])&jGeneName==gsub(".rda","",unlist(strsplit(f,split="_"))[2]),.N,])})
#names(sizevec_f)<-list.files(pattern = "resTCRBV")
#CMVinds<-which(sapply(names(rob_resl[[1]]),function(x)sum(sapply(CMVind1,grepl,x = x))!=0))
#rob_resls2<-lapply(list.files(pattern = "resTCRBV"),function(f){print(f);if(nrow(rob_resl[[f]])!=0)add_shape(resRob=rob_resl[[f]],indvec = CMV_inds,sizevec = sizevec_f[[f]])})
#Robl<-read_rob("Rob_cmv_short/")
#resq<-numeric(97);for (i in 1:97){if (!is.null(nrow(rob_resls2[[i]]))){resq[i]<-get_q(rob_resls2[[i]])}}
# get shapes for CMV ind cohort 1