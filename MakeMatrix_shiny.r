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

make_prematrix_gb<-function(LISTS,Genome){
  
  #uploading gene annotation data
  if (Genome %in% "hg19") gb=fread("refGene_hg19_gb_sorted")
  if (Genome %in% "hg38") gb=fread("refGene_hg38_gb_sorted")
  if (Genome %in% "mm10") gb=fread("refGene_mm10_gb_sorted")
  if (Genome %in% "mm9")  gb=fread("refGene_mm9_gb_sorted")
  
  #generation of GRanges oject for gene annotation data  
  colnames(gb)=c("chr","start","end","strand","symbol")
  annotations_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
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



make_prematrix_1stIntron<-function(LISTS,Genome){
  
  #uploading gene annotation data
  if (Genome %in% "hg19") gb=fread("refGene_hg19_1stIntron")
  if (Genome %in% "hg38") gb=fread("refGene_hg38_1stIntron")
  if (Genome %in% "mm10") gb=fread("refGene_mm10_1stIntron")
  if (Genome %in% "mm9")  gb=fread("refGene_mm9_1stIntron")
  
  #generation of GRanges oject for gene annotation data  
  colnames(gb)=c("chr","start","end","strand","symbol")
  annotations_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
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


make_prematrix_cgi<-function(LISTS,Genome){
  
  
  #uploading CGI annotation data
  if (Genome %in% "hg19") cgi=fread("hg19_CGI")
  if (Genome %in% "hg38") cgi=fread("hg38_CGI")
  if (Genome %in% "mm10") cgi=fread("mm10_CGI")
  if (Genome %in% "mm9")  cgi=fread("mm9_CGI")
  
  #generation of GRanges oject for CGI annotation data  
  colnames(cgi)=c("chr","start","end")
  annotations_gr=with(cgi, GRanges(chr, IRanges(start, end)))
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




make_prematrix_promoter<-function(LISTS,Genome,starT,enD){
  
  #uploading gene annotation data
  if (Genome %in% "hg19") gb=fread("refGene_hg19_gb_sorted")
  if (Genome %in% "hg38") gb=fread("refGene_hg38_gb_sorted")
  if (Genome %in% "mm10") gb=fread("refGene_mm10_gb_sorted")
  if (Genome %in% "mm9")  gb=fread("refGene_mm9_gb_sorted")
  
  colnames(gb)=c("chr","start","end","strand","symbol")
  annotations_gb_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
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


make_matrix_infinium <- function(LISTS,Genome,PLATFORM,FEATURE){
  
  #loading genome coordinates of probe IDs
  if (Genome %in% "hg19"){
    if (PLATFORM ==2 ) Infinium=readRDS("EPIC.hg19.manifest.addressA.rds") ##EPIC
    if (PLATFORM ==3 ) Infinium=readRDS("hm450.hg19.manifest.addressA.rds") ##450K
  }
  
  if (Genome %in% "hg38"){
    if (PLATFORM ==2 ) Infinium=readRDS("EPIC.hg38.manifest.addressA.rds") ##EPIC
    if (PLATFORM ==3 ) Infinium=readRDS("hm450.hg38.manifest.addressA.rds") ##450K
  }
  
  
  #loading genomic feature coordinates
  if (FEATURE == 1){
    #loading CGI annotation data
    if (Genome %in% "hg19") cgi=fread("hg19_CGI")
    if (Genome %in% "hg38") cgi=fread("hg38_CGI")
    
    colnames(cgi)=c("chr","start","end")
    annotations_gr=with(cgi, GRanges(chr, IRanges(start, end)))
    annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
    
  }
  
  if (FEATURE == 2){
    #loading gb annotation data
    if (Genome %in% "hg19") gb=fread("refGene_hg19_gb_sorted")
    if (Genome %in% "hg38") gb=fread("refGene_hg38_gb_sorted")
    
    colnames(gb)=c("chr","start","end","strand","symbol")
    annotations_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
    annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  }
  
  if (FEATURE == 3){
    #loading 1st intron annotation data
    if (Genome %in% "hg19") gb=fread("refGene_hg19_1stIntron")
    if (Genome %in% "hg38") gb=fread("refGene_hg38_1stIntron")
    
    colnames(gb)=c("chr","start","end","strand","symbol")
    annotations_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
    annotations_gr=subset(annotations_gr,!seqnames %in% "chrX" & !seqnames %in% "chrY")   ###removed SexChr
  }
  
  
  #loading methylome data (matrix)
  dm0=fread(LISTS[1])
  dm1_gr=Infinium[names(Infinium) %in% dm0$ID]
  
  fO=findOverlaps(dm1_gr,annotations_gr,ignore.strand=TRUE) ###gr=query, annotations_gr=subject
  dm1_annotation=dm1_gr[attributes(fO)$from,]                  ### selection of data
  annotation_coordinates=annotations_gr[attributes(fO)$to,] ### selection of annotation
  
  dm1_annotation=data.frame(ranges(dm1_annotation),ID_REF=names(dm1_annotation))
  
  z=cbind(as.data.frame(dm1_annotation),as.data.frame(annotation_coordinates))
  
  if (FEATURE == 1){
    
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
  
  if (FEATURE > 1){
    
    #z2=unique(z[,c(1,5,6,7,10)])
    z2=unique(z[,c(5,6,1,2,11)])
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
  
}


make_matrix_infinium_promoter <- function(LISTS,Genome,PLATFORM,starT,enD){
  
  #loading genome coordinates of probe IDs
  if (Genome %in% "hg19"){
    if (PLATFORM ==2 ) Infinium=readRDS("EPIC.hg19.manifest.addressA.rds") ##EPIC
    if (PLATFORM ==3 ) Infinium=readRDS("hm450.hg19.manifest.addressA.rds") ##450K
  }
  
  if (Genome %in% "hg38"){
    if (PLATFORM ==2 ) Infinium=readRDS("EPIC.hg38.manifest.addressA.rds") ##EPIC
    if (PLATFORM ==3 ) Infinium=readRDS("hm450.hg38.manifest.addressA.rds") ##450K
  }
  
  
  if (Genome %in% "hg19") gb=fread("refGene_hg19_gb_sorted")
  if (Genome %in% "hg38") gb=fread("refGene_hg38_gb_sorted")
  
  
  colnames(gb)=c("chr","start","end","strand","symbol")
  annotations_gb_gr=with(gb, GRanges(chr, IRanges(start, end), strand=strand, symbol=symbol))
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
  
  #z2=unique(z[,c(1,5,6,7,10)])
  z2=unique(z[,c(5,6,1,2,11)])
  
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

