library(GenomicRanges)
library(parallel)
library(data.table)
library(dplyr)
library(readr)

####################################
# 1.make_prematrix_gb 
# 2.make_prematrix_cgi
# 3.make_prematrix_promoter
# 4.make_matrix_infinium
# 5.make_matrix_infinium_promoter
####################################

make_prematrix_gb<-function(LISTS,gf){
  
  
  #generation of GRanges oject for gene annotation data  
  colnames(gf)=c("chr","start","end","strand","symbol")
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  
  d=c(1:length(LISTS))
  Cores=round(detectCores()/2,0)+1
  
  preMatrix<-mclapply(d, function(d) {
    gc(); gc(); gc();
    df0=fread(LISTS[d])
    colnames(df0)=c("chr","start","end","strand","mC","C")
    gr=with(df0, GRanges(chr, IRanges(start, end), mC=mC, C=C, strand=strand))
    
    fO=findOverlaps(gr,annotations_gr,ignore.strand=TRUE)
    selectedregions=cbind(df0[attributes(fO)$from,])
    x=annotations_gr[attributes(fO)$to,]
    xx=as.data.frame(x)
    
    y=df0[attributes(fO)$from,]
    
    z=cbind(y,xx)
    colnames(z)[8]="Start"
    colnames(z)[9]="End"
    z=data.frame(z,loci=paste(z$seqnames,z$Start,z$End,z$symbol,sep="__"),meth=z$mC/z$C*100,blockID=d)
    z_group <- z %>% group_by(loci)
    z_group_methR=z_group %>% summarise(avg_methR = round(mean(meth),1),blocK=median(blockID))
  },mc.cores=Cores)
  
  datalist=NULL;names_list=NULL;datanames=NULL
  for (i in 1:length(preMatrix)){
    blockID=as.numeric(unique(preMatrix[[i]][,3]))
    dataname=LISTS[blockID]
    datanames=c(datanames,dataname)
    
    X_names=as.matrix(preMatrix[[i]]$loci)   ##for chr split
    names_list=c(names_list,X_names)
    X=preMatrix[[i]]$avg_methR
    names(X)=X_names
    datalist=c(datalist,list(na.omit(X)))
    names(datalist)[i]=dataname
  }
  
  x=table(names_list)
  allentry=names(x[x==length(datalist)])
  
  datamat=NULL
  
  for (i in 1:length(datalist)){
    datamat=cbind(datamat,datalist[[i]][allentry])
  }
  colnames(datamat)=datanames
  
  
  targetloci=matrix(unlist(strsplit(rownames(datamat), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],datamat)
  
  #datamat=data.frame(loci=rownames(datamat),datamat)
  return(datamat)
  
}



make_prematrix_1stIntron<-function(LISTS,gf){
  
  #generation of GRanges oject for gene annotation data  
  colnames(gf)=c("chr","start","end","strand","symbol")
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  
  d=c(1:length(LISTS))
  Cores=round(detectCores()/2,0)+1
  
  preMatrix<-mclapply(d, function(d) {
    gc(); gc(); gc();
    df0=fread(LISTS[d])
    colnames(df0)=c("chr","start","end","strand","mC","C")
    gr=with(df0, GRanges(chr, IRanges(start, end), mC=mC, C=C, strand=strand))
    
    fO=findOverlaps(gr,annotations_gr,ignore.strand=TRUE)
    selectedregions=cbind(df0[attributes(fO)$from,])
    x=annotations_gr[attributes(fO)$to,]
    xx=as.data.frame(x)
    
    y=df0[attributes(fO)$from,]
    
    z=cbind(y,xx)
    colnames(z)[8]="Start"
    colnames(z)[9]="End"
    z=data.frame(z,loci=paste(z$seqnames,z$Start,z$End,z$symbol,sep="__"),meth=z$mC/z$C*100,blockID=d)
    z_group <- z %>% group_by(loci)
    z_group_methR=z_group %>% summarise(avg_methR = round(mean(meth),1),blocK=median(blockID))
  },mc.cores=Cores)
  
  datalist=NULL;names_list=NULL;datanames=NULL
  for (i in 1:length(preMatrix)){
    blockID=as.numeric(unique(preMatrix[[i]][,3]))
    dataname=LISTS[blockID]
    datanames=c(datanames,dataname)
    
    X_names=as.matrix(preMatrix[[i]]$loci)   ##for chr split
    names_list=c(names_list,X_names)
    X=preMatrix[[i]]$avg_methR
    names(X)=X_names
    datalist=c(datalist,list(na.omit(X)))
    names(datalist)[i]=dataname
  }
  
  x=table(names_list)
  allentry=names(x[x==length(datalist)])
  
  datamat=NULL
  
  for (i in 1:length(datalist)){
    datamat=cbind(datamat,datalist[[i]][allentry])
  }
  colnames(datamat)=datanames
  
  
  targetloci=matrix(unlist(strsplit(rownames(datamat), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],datamat)
  
  #datamat=data.frame(loci=rownames(datamat),datamat)
  return(datamat)
  
}


make_prematrix_cgi<-function(LISTS,gf){
  
  #generation of GRanges oject for CGI annotation data  
  colnames(gf)=c("chr","start","end")
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end)))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  d=c(1:length(LISTS))
  Cores=round(detectCores()/2,0)+1
  
  preMatrix<-mclapply(d, function(d) {
    gc(); gc(); gc();
    df0=fread(LISTS[d])
    colnames(df0)=c("chr","start","end","strand","mC","C")
    gr=with(df0, GRanges(chr, IRanges(start, end), mC=mC, C=C, strand=strand))
    
    fO=findOverlaps(gr,annotations_gr,ignore.strand=TRUE)
    selectedregions=cbind(df0[attributes(fO)$from,])
    x=annotations_gr[attributes(fO)$to,]
    xx=as.data.frame(x)
    
    y=df0[attributes(fO)$from,]
    
    z=cbind(y,xx)
    colnames(z)[8]="Start"
    colnames(z)[9]="End"
    z=data.frame(z,loci=paste(z$seqnames,z$Start,z$End,"CGI",sep="__"),meth=z$mC/z$C*100,blockID=d)
    z_group <- z %>% group_by(loci)
    z_group_methR=z_group %>% summarise(avg_methR = round(mean(meth),1),blocK=median(blockID))
  },mc.cores=Cores)
  
  
  datalist=NULL;names_list=NULL;datanames=NULL
  for (i in 1:length(preMatrix)){
    blockID=as.numeric(unique(preMatrix[[i]][,3]))
    dataname=LISTS[blockID]
    datanames=c(datanames,dataname)
    
    X_names=as.matrix(preMatrix[[i]]$loci)   ##for chr split
    names_list=c(names_list,X_names)
    X=preMatrix[[i]]$avg_methR
    names(X)=X_names
    datalist=c(datalist,list(na.omit(X)))
    names(datalist)[i]=dataname
  }
  
  x=table(names_list)
  allentry=names(x[x==length(datalist)])
  
  datamat=NULL
  
  for (i in 1:length(datalist)){
    datamat=cbind(datamat,datalist[[i]][allentry])
  }
  colnames(datamat)=datanames
  
  targetloci=matrix(unlist(strsplit(rownames(datamat), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],datamat)
  
  return(datamat)
  
}




make_prematrix_promoter<-function(LISTS,gf,starT,enD){
  
  colnames(gf)=c("chr","start","end","strand","symbol")
  annotations_gb_gr=with(gf, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
  annotations_gb_gr=subset(annotations_gb_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  annotations_gr=promoters(annotations_gb_gr,upstream=abs(starT),downstream=abs(enD))
  
  d=c(1:length(LISTS))
  Cores=round(detectCores()/2,0)+1
  
  preMatrix<-mclapply(d, function(d) {
    gc(); gc(); gc();
    df0=fread(LISTS[d])
    colnames(df0)=c("chr","start","end","strand","mC","C")
    gr=with(df0, GRanges(chr, IRanges(start, end), mC=mC, C=C, strand=strand))
    
    fO=findOverlaps(gr,annotations_gr,ignore.strand=TRUE)
    selectedregions=cbind(df0[attributes(fO)$from,])
    x=annotations_gr[attributes(fO)$to,]
    xx=as.data.frame(x)
    
    y=df0[attributes(fO)$from,]
    
    z=cbind(y,xx)
    colnames(z)[8]="Start"
    colnames(z)[9]="End"
    z=data.frame(z,loci=paste(z$seqnames,z$Start,z$End,z$symbol,sep="__"),meth=z$mC/z$C*100,blockID=d)
    z_group <- z %>% group_by(loci)
    z_group_methR=z_group %>% summarise(avg_methR = round(mean(meth),1),blocK=median(blockID))
  },mc.cores=Cores)
  
  
  datalist=NULL;names_list=NULL;datanames=NULL
  for (i in 1:length(preMatrix)){
    blockID=as.numeric(unique(preMatrix[[i]][,3]))
    dataname=LISTS[blockID]
    datanames=c(datanames,dataname)
    
    X_names=as.matrix(preMatrix[[i]]$loci)   ##for chr split
    names_list=c(names_list,X_names)
    X=preMatrix[[i]]$avg_methR
    names(X)=X_names
    datalist=c(datalist,list(na.omit(X)))
    names(datalist)[i]=dataname
  }
  
  x=table(names_list)
  allentry=names(x[x==length(datalist)])
  
  datamat=NULL
  
  for (i in 1:length(datalist)){
    datamat=cbind(datamat,datalist[[i]][allentry])
  }
  colnames(datamat)=datanames
  
  
  targetloci=matrix(unlist(strsplit(rownames(datamat), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],datamat)
  
  return(datamat)
  
  
}



make_matrix_infinium_cgi <- function(LISTS,gf,Infinium,FEATURE){

  #loading genomic feature coordinates
  colnames(gf)=c("chr","start","end")
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end)))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  #loading methylome data (matrix)
  dm0=fread(LISTS[1])
  dm1_gr=Infinium[names(Infinium) %in% dm0$ID]
  
  fO=findOverlaps(dm1_gr,annotations_gr,ignore.strand=TRUE) ###gr=query, annotations_gr=subject
  dm1_annotation=dm1_gr[attributes(fO)$from,]               ### selection of data
  annotation_coordinates=annotations_gr[attributes(fO)$to,] ### selection of annotation
  
  dm1_annotation=data.frame(ranges(dm1_annotation),ID_REF=names(dm1_annotation))
  
  z=cbind(as.data.frame(dm1_annotation),as.data.frame(annotation_coordinates))
  

    
    #z2=unique(z[,c(1,5,6,7)])
    z2=unique(z[,c(5,6,7,8)])
    xxx2=merge(z2,dm0,by.x="ID_REF",by.y="ID",all = F)
    xxx2=data.frame(loci=paste(xxx2$seqnames,xxx2$start,xxx2$end,"CGI",sep="__"),xxx2)
    xxx3=xxx2[,c(1,6:ncol(xxx2))]
    
    xxx3_naomit=na.omit(xxx3)
    xxx4<-xxx3_naomit %>%
      dplyr::group_by(loci) %>%
      dplyr::summarise_at(dplyr::funs(mean),.vars=colnames(xxx3_naomit)[c(2:ncol(xxx3_naomit))])
    
    
    targetloci=matrix(unlist(strsplit(as.matrix(xxx4[,1]), "__")),ncol=4,byrow=T)
    datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],round(xxx4[,c(2:ncol(xxx4))],2))
    return(na.omit(datamat))
    
  
}



make_matrix_infinium_cgishores <- function(LISTS,gf,Infinium,FEATURE){
  
  #loading genomic feature coordinates
  colnames(gf)=c("chr","start","end")
  
  gf_shores1=gf
  gf_shores1$start=gf$start - 2000
  gf_shores1$end  =gf$start - 1
  
  gf_shores2=gf
  gf_shores2$start=gf$end +1
  gf_shores2$end  =gf$end +2000
  
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end)))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  annotations_gr_shores1=with(gf_shores1, GRanges(chr, IRanges(start, end)))
  annotations_gr_shores1=subset(annotations_gr_shores1,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  annotations_gr_shores2=with(gf_shores2, GRanges(chr, IRanges(start, end)))
  annotations_gr_shores2=subset(annotations_gr_shores2,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  
  #loading methylome data (matrix)
  dm0=fread(LISTS[1])
  dm1_gr=Infinium[names(Infinium) %in% dm0$ID]
  
  fO=findOverlaps(dm1_gr,annotations_gr,ignore.strand=TRUE) ###gr=query, annotations_gr=subject
  fO_shore1=findOverlaps(dm1_gr,annotations_gr_shores1,ignore.strand=TRUE)
  fO_shore2=findOverlaps(dm1_gr,annotations_gr_shores2,ignore.strand=TRUE)
  
  dm1_annotation=dm1_gr[attributes(fO)$from,]               ### selection of data
  y1=dm1_gr[attributes(fO_shore1)$from,]
  y2=dm1_gr[attributes(fO_shore2)$from,]
  
  annotation_coordinates =annotations_gr[attributes(fO)$to,]         ## selection of annotation
  annotation_coordinates1=annotations_gr[attributes(fO_shore1)$to,]  #######annotations_gr 
  annotation_coordinates2=annotations_gr[attributes(fO_shore2)$to,]  #######annotations_gr
  
  dm1_annotation=data.frame(ranges(dm1_annotation),ID_REF=names(dm1_annotation))
  y1=data.frame(ranges(y1),ID_REF=names(y1))
  y2=data.frame(ranges(y2),ID_REF=names(y2))
  
  z =cbind(as.data.frame(dm1_annotation),as.data.frame(annotation_coordinates))
  z1=cbind(as.data.frame(y1),as.data.frame(annotation_coordinates1))
  z2=cbind(as.data.frame(y2),as.data.frame(annotation_coordinates2))
  z1_z2=rbind(z1,z2)
  
  
  z1_z2_uniq=unique(z1_z2[,c(5,6,7,8)])
  xxx1_xxx2=merge(z1_z2_uniq,dm0,by.x="ID_REF",by.y="ID",all = F)
  xxx1_xxx2=data.frame(loci=paste(xxx1_xxx2$seqnames,xxx1_xxx2$start,xxx1_xxx2$end,"CGI",sep="__"),xxx1_xxx2)
  yyy1_yyy2=xxx1_xxx2[,c(1,6:ncol(xxx1_xxx2))]
  yyy1_yyy2_naomit=na.omit(yyy1_yyy2)
  zzz1_zzz2<-yyy1_yyy2_naomit %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise_at(dplyr::funs(mean),.vars=colnames(yyy1_yyy2_naomit)[c(2:ncol(yyy1_yyy2_naomit))])
  
  targetloci=matrix(unlist(strsplit(as.matrix(zzz1_zzz2[,1]), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],round(zzz1_zzz2[,c(2:ncol(zzz1_zzz2))],2))
  
  #z2=unique(z[,c(5,6,7,8)])
  #xxx2=merge(z2,dm0,by.x="ID_REF",by.y="ID",all = F)
  #xxx2=data.frame(loci=paste(xxx2$seqnames,xxx2$start,xxx2$end,"CGI",sep="__"),xxx2)
  #xxx3=xxx2[,c(1,6:ncol(xxx2))]
  
  #xxx3_naomit=na.omit(xxx3)
  #xxx4<-xxx3_naomit %>%
  #  dplyr::group_by(loci) %>%
  #  dplyr::summarise_at(dplyr::funs(mean),.vars=colnames(xxx3_naomit)[c(2:ncol(xxx3_naomit))])
  
  #targetloci=matrix(unlist(strsplit(as.matrix(xxx4[,1]), "__")),ncol=4,byrow=T)
  #datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],round(xxx4[,c(2:ncol(xxx4))],2))
  
  
   return(na.omit(datamat))
  
  
}



#gene body or 1st intron
make_matrix_infinium <- function(LISTS,gf,Infinium,FEATURE){
  
  #loading genomic feature coordinates
  colnames(gf)=c("chr","start","end","strand","symbol")
  annotations_gr=with(gf, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
  annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  
  #loading methylome data (matrix)
  dm0=fread(LISTS[1])
  dm1_gr=Infinium[names(Infinium) %in% dm0$ID]
  
  fO=findOverlaps(dm1_gr,annotations_gr,ignore.strand=TRUE) ###gr=query, annotations_gr=subject
  dm1_annotation=dm1_gr[attributes(fO)$from,]               ### selection of data
  annotation_coordinates=annotations_gr[attributes(fO)$to,] ### selection of annotation
  
  dm1_annotation=data.frame(ranges(dm1_annotation),ID_REF=names(dm1_annotation))
  
  z=cbind(as.data.frame(dm1_annotation),as.data.frame(annotation_coordinates))
  
    
    #z2=unique(z[,c(1,5,6,7,10)])
    z2=unique(z[,c(5,6,7,8,11)])
    xxx2=merge(z2,dm0,by.x="ID_REF",by.y="ID",all = F)
    xxx2=data.frame(loci=paste(xxx2$seqnames,xxx2$start,xxx2$end,xxx2$symbol,sep="__"),xxx2)
    xxx3=xxx2[,c(1,7:ncol(xxx2))]
    
    xxx3_naomit=na.omit(xxx3)
    xxx4<-xxx3_naomit %>%
      dplyr::group_by(loci) %>%
      dplyr::summarise_at(dplyr::funs(mean),.vars=colnames(xxx3_naomit)[c(2:ncol(xxx3_naomit))])
    
    targetloci=matrix(unlist(strsplit(as.matrix(xxx4[,1]), "__")),ncol=4,byrow=T)
    datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],round(xxx4[,c(2:ncol(xxx4))],2))
    return(na.omit(datamat))
  
}

make_matrix_infinium_promoter <- function(LISTS,gf,Infinium,starT,enD){
  
  #loading genomic feature coordinates
  colnames(gf)=c("chr","start","end","strand","symbol")
  annotations_gb_gr=with(gf, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
  annotations_gb_gr=subset(annotations_gb_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  annotations_gr=promoters(annotations_gb_gr,upstream=abs(starT),downstream=abs(enD))
  
  #loading methylome data (matrix)
  dm0=fread(LISTS[1])
  dm1_gr=Infinium[names(Infinium) %in% dm0$ID]
  
  fO=findOverlaps(dm1_gr,annotations_gr,ignore.strand=TRUE) ###gr=query, annotations_gr=subject
  dm1_annotation=dm1_gr[attributes(fO)$from,]                  ### selection of data
  annotation_coordinates=annotations_gr[attributes(fO)$to,] ### selection of annotation
  
  dm1_annotation=data.frame(ranges(dm1_annotation),ID_REF=names(dm1_annotation))
  z=cbind(as.data.frame(dm1_annotation),as.data.frame(annotation_coordinates))
  
  z2=unique(z[,c(5,6,7,8,11)])
  
  xxx2=merge(z2,dm0,by.x="ID_REF",by.y="ID",all = F)
  xxx2=data.frame(loci=paste(xxx2$seqnames,xxx2$start,xxx2$end,xxx2$symbol,sep="__"),xxx2)
  xxx3=xxx2[,c(1,7:ncol(xxx2))]
  
  xxx3_naomit=na.omit(xxx3)
  xxx4<-xxx3_naomit %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise_at(dplyr::funs(mean),.vars=colnames(xxx3_naomit)[c(2:ncol(xxx3_naomit))])
  
  targetloci=matrix(unlist(strsplit(as.matrix(xxx4[,1]), "__")),ncol=4,byrow=T)
  datamat=data.frame(chr=targetloci[,1],Start=targetloci[,2],End=targetloci[,3],Symbol=targetloci[,4],round(xxx4[,c(2:ncol(xxx4))],2))
  return(na.omit(datamat))
  
}


