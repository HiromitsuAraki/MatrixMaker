# Load packages ----
library(shiny)
library(shinyFiles)
library(GenomicRanges)
library(parallel)
library(data.table)
library(dplyr)
library(shinydashboard)


# Source MakeMatrix----
source("MakeMatrix_shiny.R")


###UI
shinyUI(dashboardPage(skin="green",
                      dashboardHeader(title = "MatrixMaker"),
                      dashboardSidebar(
                        sidebarMenu(
                          menuItem("Parameters",  tabName = "parameters", icon = icon("th")),
                          menuItem("README",      tabName = "readme",     icon = icon("th"))
                        )
                      ),
                      
                      dashboardBody(
                        tabItems(
                          # First tab content
                          tabItem(tabName = "parameters",
                                  fluidRow(
                                    column(3,
                                           h3("Select Directory"),
                                           shinyDirButton('directory','Select Directory:','Select Directory:')),
                                    
                                    column(3, 
                                           radioButtons("checkPlatform", 
                                                        label = h3("Select Platform"), 
                                                        choices = list("Bisulfite seq"  = 1, 
                                                                       "Infinium (EPIC)"= 2,
                                                                       "Infinium (450K)"= 3),
                                                        selected = "")),
                                    
                                    column(3, 
                                           radioButtons("checkGenome", 
                                                        label = h3("Select Genome"), 
                                                        choices = list(
                                                          "hg38" = "hg38", 
                                                          "hg19" = "hg19", 
                                                          "mm10" = "mm10",
                                                          "mm9"  = "mm9"),
                                                        selected = "")
                                    ),                     
                                    
                                    column(3,
                                           radioButtons("checkFeature", 
                                                        label = h3("Select Feature"), 
                                                        choices = list("CGI"          = 1, 
                                                                       "CGI shore"    = 2, 
                                                                       "Gene body"    = 3, 
                                                                       "1st intron"   = 4, 
                                                                       "Promoter*"    = 5),
                                                        
                                                        
                                                        selected = ""),
                                           tags$span(
                                             h5("*Promoter region (TSS=0)")
                                           ),
                                           fluidRow(
                                             column(4, numericInput("PromoterStart", label = h6("Start (-10000~-1)"), value="",width='100%',min = (-10000), max =(-1))),
                                             column(4, numericInput("PromoterEnd",   label = h6("End (0~10000)"),     value="",width='100%',min = 0, max =10000))
                                           ),
                                           br(),
                                           br(),
                                           
                                           h3("Submit"),
                                           actionButton('submit','Submit'),
                                           
                                           br(),
                                           h3("Download"),
                                           downloadButton("MethMat", "Download")
                                    )
                                  )
                          ),
                          
                          
                          tabItem(tabName = "readme",
                                  h2("README")
                          )
                        )
                      )
)
)

  
