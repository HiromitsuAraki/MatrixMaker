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
                 
                 #output$selected_var <- renderText({
                 #   updateRadioButtons(session, "checkFeature",  value = input$checkFeature)
                 #   updateNumericInput(session, "PromoterStart", value = input$PromoterStart)
                 #   updateNumericInput(session, "PromoterEnd",   value = input$PromoterEnd)
                 #   paste(paste("Feature", input$checkFeature),paste("Start", input$PromoterStart),paste("End", input$PromoterEnd))
                 # })
                 
                 ##ここから
                 
                 
                 roots=c('/'='/')
                 mydirectory = parseDirPath(roots, input$directory)
                 setwd(mydirectory)
                 myfiles = list.files(mydirectory)
                 #Myfiles=paste(mydirectory,myfiles,sep="/")
                 
                 #output$MethMat=renderDataTable(make_prematrix(myfiles))
                 
                 withProgress(message = 'In progress', {
                   #MethMat=make_prematrix(myfiles)
                   if (input$checkPlatform==1){
                     if (input$checkFeature==1){
                       MethMat=make_prematrix_cgi(myfiles,as.character(input$checkGenome))
                     }
                     
                     if (input$checkFeature==2){
                       MethMat=make_prematrix_gb(myfiles,as.character(input$checkGenome))                     
                     }

                     if (input$checkFeature==3){
                       MethMat=make_prematrix_1stIntron(myfiles,as.character(input$checkGenome))                     
                     }
                     
                     if (input$checkFeature==4){
                       #updateNumericInput(session, "PromoterStart", value = input$PromoterStart)
                       #updateNumericInput(session, "PromoterEnd",   value = input$PromoterEnd)
                       MethMat=make_prematrix_promoter(myfiles,as.character(input$checkGenome),input$PromoterStart,input$PromoterEnd)                     
                     }
                     
                   }
                   
                   if (input$checkPlatform>=2){
                     if (input$checkFeature<=3){
                       MethMat=make_matrix_infinium(myfiles,as.character(input$checkGenome),input$checkPlatform,input$checkFeature)                     
                     }
                     
                     if (input$checkFeature==4){
                       #updateNumericInput(session, "PromoterStart", value = input$PromoterStart)
                       #updateNumericInput(session, "PromoterEnd",   value = input$PromoterEnd)
                       MethMat=make_matrix_infinium_promoter(myfiles,as.character(input$checkGenome),input$checkPlatform, input$PromoterStart,input$PromoterEnd)                     
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


