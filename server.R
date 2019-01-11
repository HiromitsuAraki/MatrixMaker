# Load packages ----
library(shiny)
library(shinyFiles)
library(GenomicRanges)
library(parallel)
library(data.table)
library(dplyr)
library(shinydashboard)

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))}
  inputs
}

# Source MakeMatrix ----
source("MakeMatrix_shiny.R")

###Server
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=3000*1024^2) 
  
  roots=c('/'='/')
  shinyDirChoose(input,'directory',roots=c('/'='/'))

  observeEvent(input$submit,
               {
                 
                 withProgress(message = 'In progress', {

                   #annotation load
                   if (as.character(input$checkGenome)=="hg38"){
                     if (input$checkFeature==1)  gf=fread("hg38_CGI")
                     if (input$checkFeature==2)  gf=fread("refGene_hg38_gb_sorted")
                     if (input$checkFeature>=3)  gf=fread("refGene_hg38_1stIntron")
                     if (input$checkPlatform==2) Infinium=readRDS("EPIC.hg38.manifest.addressA.rds")
                     if (input$checkPlatform==3) Infinium=readRDS("hm450.hg38.manifest.addressA.rds")
                   }
                   
                   if (as.character(input$checkGenome)=="hg19"){
                     if (input$checkFeature==1)  gf=fread("hg19_CGI")
                     if (input$checkFeature==2)  gf=fread("refGene_hg19_gb_sorted")
                     if (input$checkFeature>=3)  gf=fread("refGene_hg19_1stIntron")
                     if (input$checkPlatform==2) Infinium=readRDS("EPIC.hg19.manifest.addressA.rds")
                     if (input$checkPlatform==3) Infinium=readRDS("hm450.hg19.manifest.addressA.rds")
                   }
      
                   if (as.character(input$checkGenome)=="mm10"){
                     if (input$checkFeature==1) gf=fread("mm10_CGI")
                     if (input$checkFeature==2) gf=fread("refGene_mm10_gb_sorted")
                     if (input$checkFeature>=3) gf=fread("refGene_mm10_1stIntron")
                   }
                  
                   if (as.character(input$checkGenome)=="mm9"){
                     if (input$checkFeature==1) gf=fread("mm9_CGI")
                     if (input$checkFeature==2) gf=fread("refGene_mm9_gb_sorted")
                     if (input$checkFeature>=3) gf=fread("refGene_mm9_1stIntron")
                   }

                   
                   #methylome data load
                   roots=c('/'='/')
                   mydirectory = parseDirPath(roots, input$directory)
                   setwd(mydirectory)
                   myfiles = list.files(mydirectory)
                   
                   if (input$checkPlatform==1){
                     if (input$checkFeature==1){
                       #MethMat=make_prematrix_cgi(myfiles,as.character(input$checkGenome),gf)
                       MethMat=make_prematrix_cgi(myfiles,gf)
                     }
                     
                     if (input$checkFeature==2){
                       #MethMat=make_prematrix_gb(myfiles,as.character(input$checkGenome),gf)
                       MethMat=make_prematrix_gb(myfiles,gf)   
                     }

                     if (input$checkFeature==3){
                       #MethMat=make_prematrix_1stIntron(myfiles,as.character(input$checkGenome),gf)
                       MethMat=make_prematrix_1stIntron(myfiles,gf)
                     }
                     
                     if (input$checkFeature==4){
                       #MethMat=make_prematrix_promoter(myfiles,as.character(input$checkGenome),gf,input$PromoterStart,input$PromoterEnd)
                       MethMat=make_prematrix_promoter(myfiles,gf,input$PromoterStart,input$PromoterEnd)
                     }
                     
                   }
                   
                   if (input$checkPlatform>=2){
                     if (input$checkFeature<=3){
                       #MethMat=make_matrix_infinium(myfiles,as.character(input$checkGenome),gf,input$checkPlatform,input$checkFeature)
                       MethMat=make_matrix_infinium(myfiles,gf,Infinium,input$checkFeature)
                     }
                     
                     if (input$checkFeature==4){
                       #MethMat=make_matrix_infinium_promoter(myfiles,as.character(input$checkGenome),gf,input$checkPlatform, input$PromoterStart,input$PromoterEnd)
                       MethMat=make_matrix_infinium_promoter(myfiles, gf, Infinium,input$PromoterStart,input$PromoterEnd)
                     }
                   }
                   
                   
                 })
                 
                 output$MethMat <- downloadHandler(
                   filename = function() {
                     paste("data-", Sys.time(), ".txt", sep="")
                   },
                   content = function(file) {
                     write.table(MethMat,file,sep="\t",row.names=F,quote=F,col.names=T)
                   }
                 )
                 
               }
  )
  
})


