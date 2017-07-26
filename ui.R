library(shiny)
library(shinydashboard)
library(dygraphs)
library(d3heatmap)
library(colourpicker)


scriptHeaders <- function(inputId=NULL) {
  tagList(
    singleton(tags$head(tags$script(src = "js/gtm.js")))
  )
}

textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

shinyUI(fluidPage(theme = "bootstrap.css",
                  
             headerPanel(
               h1("shinyGISPA", 
                  style = "font-family: 'Times New Roman'; font-style:bold;
                  font-size: 80px; line-height: 1.1;
                  color: red")),
             
             tags$h4("Gene Integrated Set Profile Analysis with Shiny", 
                     style = "font-family: 'Times New Roman'; font-style:italic;
                     font-size: 50px;
                     color: red"),
             
             

             fluidRow(
               column(3,
                      wellPanel(
                        h4("Analysis Type"),
                        tags$hr(),
                        radioButtons("analysisType", "", choices = c("Single" = "1d","Two" = "2d","Three" = "3d"), inline=TRUE)
                      ),
                      fluidRow(
                        column(6,
                          wellPanel(
                            conditionalPanel(condition = "input.analysisType == '1d'",
                                   h4("Data Input"),
                                   tags$hr(),
                                   selectInput(inputId="file", label="File Input:", choices=c("User data"="LD","Example data"="ED"), selectize = TRUE, width=NULL),
                                   conditionalPanel("input.file=='LD'", fileInput(inputId="dataset", label="", width=NULL, multiple=TRUE)),
                                   selectInput("profile","GeneSet Profile", c('up','down'),width="130px",selectize = TRUE),
                                   #selectInput(inputId="dataType","Data Type", c('EXP','VAR','CNV','MET','OTHER'), width="130px",selectize = TRUE),
                                   #conditionalPanel(condition = "input.dataType == 'OTHER'",
                                   textInput("dataType", "Data Type", value="", width="149px")
                                   #)
                                                    
                            ),
                            conditionalPanel(condition = "input.analysisType == '2d' || input.analysisType == '3d'",
                                   h4("Data Input 1"),
                                   tags$hr(),
                                   selectInput(inputId="file1", label="File Input 1:", choices=c("User data"="LD","Example data"="ED"), selectize = TRUE, width=NULL),
                                   conditionalPanel("input.file1=='LD'", fileInput(inputId="dataset1", label="", width=NULL, multiple=TRUE)),
                                   selectInput("profile1","GeneSet Profile 1", c('up','down'),width="130px",selectize = TRUE),
                                   #selectInput(inputId="dataType1","Data Type 1", c('EXP','VAR','CNV','MET','OTHER'), width="130px",selectize = TRUE),
                                   #conditionalPanel(condition = "input.dataType1 == 'OTHER'",
                                   textInput("dataType1", "Data Type 1", value="", width="149px")
                                   #)
                            )
                            
                          )),
                        column(6,
                          wellPanel(
                            conditionalPanel(condition = "input.analysisType == '1d'",
                                   h4("Define Samples"),
                                   tags$hr(),         
                                   numericInput(inputId="ref_1d", "Reference", c('3','4','5'),width="130px", value=3),
                                   numericInput(inputId="cmp1_1d", "Sample 1", c('3','4','5'), width="130px", value=4),
                                   numericInput(inputId="cmp2_1d", "Sample 2", c('3','4','5'),width="130px", value=5)
                            ),
                           conditionalPanel(condition = "input.analysisType == '2d' || input.analysisType == '3d'",
                                   h4("Data Input 2"),
                                   tags$hr(),
                                   selectInput(inputId="file2", label="File Input 2:", choices=c("User data"="LD","Example data"="ED"), selectize = TRUE, width=NULL),
                                   conditionalPanel("input.file2=='LD'", fileInput(inputId="dataset2", label="", width=NULL, multiple=TRUE)),
                                   selectInput("profile2","GeneSet Profile 2", c('up','down'),width="130px",selectize = TRUE),
                                   #selectInput(inputId="dataType2","Data Type 2", c('EXP','VAR','CNV','MET','OTHER'), width="130px",selectize = TRUE),
                                   #conditionalPanel(condition = "input.dataType2 == 'OTHER'",
                                   textInput("dataType2", "Data Type 2", value="", width="149px")
                                   #)
                           )
                        ))),
                      fluidRow(
                        column(12,
                        conditionalPanel(condition = "input.analysisType == '2d'",
                          wellPanel(
                                 h4("Define Samples"),
                                 tags$hr(),
                            fluidRow(
                                column(4,
                                        numericInput(inputId="ref_2d", "Reference", c('3','4','5'),width="100px", value=3)),
                                column(4,
                                        numericInput(inputId="cmp1_2d", "Sample 1", c('3','4','5'), width="100px", value=4)),
                                column(4,
                                        numericInput(inputId="cmp2_2d", "Sample 2", c('3','4','5'),width="100px", value=5))
                           ))
                      ))),
                     fluidRow(
                       column(6,
                        conditionalPanel(condition = "input.analysisType == '3d'",
                          wellPanel(
                                   h4("Data Input 3"),
                                   tags$hr(),
                                   selectInput(inputId="file3", label="File Input 3:", choices=c("User data"="LD","Example data"="ED"), selectize = TRUE, width=NULL),
                                   conditionalPanel("input.file3=='LD'", fileInput(inputId="dataset3", label="", width=NULL, multiple=TRUE)),
                                   selectInput("profile3","GeneSet Profile 3", c('up','down'),width="130px",selectize = TRUE),
                                   #selectInput(inputId="dataType3","Data Type 3", c('EXP','VAR','CNV','MET','OTHER'),width="130px",selectize = TRUE),
                                   #conditionalPanel(condition = "input.dataType3 == 'OTHER'",
                                   textInput("dataType3", "Data Type 3", value="", width="149px")
                                   #)
                        ))),
                       column(6,
                        conditionalPanel(condition = "input.analysisType == '3d'",
                          wellPanel(
                                   h4("Define Samples"),
                                   tags$hr(),         
                                   numericInput(inputId="ref_3d", "Reference", c('3','4','5'),width="130px", value=3),
                                   numericInput(inputId="cmp1_3d", "Sample 1", c('3','4','5'), width="130px", value=4),
                                   numericInput(inputId="cmp2_3d", "Sample 2", c('3','4','5'),width="130px", value=5)
                  )))),
                      fluidRow(
                        column(6,
                          wellPanel(
                            h4("Changepoint Input"),
                            tags$hr(),
                            selectInput("changes","Changes Using:", c('mean','var','meanvar'),width="150px",selected="var", selectize = TRUE),
                            selectInput("method","Method:", c('AMOC','PELT','BinSeg','SeqNeigh'),width="150px",selected="BinSeg", selectize = TRUE),
                            numericInput("max","Max Allowed:", min=1, max=60, step=1, width="149px", value=60)
                          )),
                        column(6,
                            wellPanel(
                                h4("Display"),
                                tags$hr(),
                                numericInput("cpt_to_display","Upto", min=1, max=60, step=1, width="149px", value=1)
                            ),
                            wellPanel(
                              #h4("Download"),
                              #tags$hr(),
                              downloadButton('downloadGO', 'Save')
                            )
                        )),
                      fluidRow(
                        column(12,
                          conditionalPanel(condition =  "input.Onedim == 4 || input.Twodim == 4 || input.Threedim == 4",
                             wellPanel(
                               sliderInput(inputId = "opt.cex",
                                           label = "Gene Label Size",                            
                                           min = 0, max = 2, step = 0.25, value = 1,width="400px"),
                               sliderInput(inputId = "opt.cexaxis",
                                           label = "Axis Text Size",                            
                                           min = 0, max = 2, step = 0.25, value = 1,width="400px"),
                               sliderInput(inputId = "opt.gap",
                                           label = "Plot Gaps",                            
                                           min = -20, max = 20, step = 0.1, value = 1,width="400px")
                             ))
                      ))
               ),
               
               
               column(8,
                conditionalPanel(condition = "input.analysisType == '1d'",     
                  tabsetPanel(id="Onedim",
                      #tabPanel("Manual",
                      #         tags$iframe(src = "shinyGISPA_manual.pdf", width = 1200, height = 1000)
                      #),
                      tabPanel ("Input Data",value=1,
                                fluidRow(
                                  column(8,
                                      plotOutput(("boxPlot"), height="600px", width="500px"),
                                      br(),
                                      textOutput("geneNumtext"),
                                      textOutput("probeNumtext"),
                                      textOutput("sampNumtext"),
                                      br(),
                                      textOutput("refSamp"),
                                      textOutput("sampOne"),
                                      textOutput("sampTwo"),
                                      br(),
                                      br(),
                                      dataTableOutput("inputData")
                                  ))
                      ),
                      tabPanel("Results Table",value=2,
                               fluidRow(
                                 column(12,
                                        dataTableOutput("gispa_data")
                                 ))
                      ),
                      tabPanel("Dignostic Plots",value=3,
                               fluidRow(
                                 column(6,
                                        plotOutput(("gispa_plot"), height="600px", width="500px")
                                 ),
                                 column(6,
                                        plotOutput(("slope_plot"), height="600px", width="500px")
                                 )),
                               fluidRow(
                                 column(6,
                                        plotOutput(("profile_boxplot"), height="600px", width="500px")
                                 ))
                      ),
                      tabPanel("GeneSet Profile",value=4,
                               fluidRow(
                                 column(8,
                                        br(),
                                        textOutput("profile_geneNumtext"),
                                        br()
                                 ))
                      ),
                      tabPanel("How to Cite",value=5,
                               #tags$div(class="header", checked=NA,
                               h4("shinyGISPA is a product of the Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University (https://bbisr.winship.emory.edu)."),
                               br(),
                               h4("This work is funded by the Leukemia and Lymphoma Society Translational Research Program Award (Jeanne Kowalski); Georgia Research Alliance Scientist Award (Jeanne Kowalski); a Team Science Seed Funding from the Winship Cancer Institute of Emory University (Lawrence H. Boise, Sagar Lonial, Michael R. Rossi); Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University and NIH/NCI [Award number P30CA138292, in part]. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH."),
                               br(),
                               h4("Bioconductor R package for GISPA is available at the https://www.bioconductor.org/packages/GISPA"),
                               br(),
                               h4("Please cite the GISPA method as: Kowalski J, Dwivedi B, Newman S, Switchenko JM, Pauly R, Gutman DA, Arora J, Gandhi K, Ainslie K, Doho G, Qin Z, Moreno CS, Rossi MR, Vertino PM, Lonial S, Bernal-Mizrachi L, Boise LH. Gene integrated set profile analysis: a context-based approach for inferring biological endpoints. Nucleic Acids Res. 2016 Apr 20;44(7):e69. doi: 10.1093/nar/gkv1503. Epub 2016 Jan 29. PubMed PMID: 26826710; PubMed Central PMCID: PMC4838358."),
                               br(),
                               tags$iframe(src = "shinyGISPA_Manual.pdf", width = 1200, height = 1000)
                      )
                      
                )),

               conditionalPanel(condition = "input.analysisType == '2d'",     
                  tabsetPanel(id="Twodim",
                      #tabPanel("Manual",
                      #         tags$iframe(src = "shinyGISPA_manual.pdf", width = 1200, height = 1000)
                      #),
                      tabPanel ("Input Data", value=1,
                            fluidRow(
                              column(6,
                                     plotOutput(("boxPlot1"), height="600px", width="500px"),
                                     br(),
                                     textOutput("geneNumtext1"),
                                     textOutput("probeNumtext1"),
                                     textOutput("sampNumtext1"),
                                     br(),
                                     textOutput("refSamp1"),
                                     textOutput("sampOne1"),
                                     textOutput("sampTwo1"),
                                     br(),
                                     br(),
                                     dataTableOutput("inputData1")
                              ),
                              column(6,
                                     plotOutput(("boxPlot2"), height="600px", width="500px"),
                                     br(),
                                     textOutput("geneNumtext2"),
                                     textOutput("probeNumtext2"),
                                     textOutput("sampNumtext2"),
                                     br(),
                                     textOutput("refSamp2"),
                                     textOutput("sampOne2"),
                                     textOutput("sampTwo2"),
                                     br(),
                                     br(),
                                     dataTableOutput("inputData2")
                              ))
                    ),
                    tabPanel("Results Table", value=2,
                             fluidRow(
                               column(12,
                                      dataTableOutput("gispa_data2")
                               ))
                    ),
                    tabPanel("Diagnostic Plots", value=3,
                             fluidRow(
                               column(6,
                                      plotOutput(("gispa_plot2"), height="600px", width="600px")
                               ),
                               column(6,
                                      plotOutput(("slope_plot2"), height="600px", width="600px")
                               )),
                             fluidRow(
                               column(6,
                                      plotOutput(("profile_boxplot1"), height="600px", width="600px")
                               ),
                               column(6,
                                      plotOutput(("profile_boxplot2"), height="600px", width="600px")
                               ))
                    ),
                    tabPanel("GeneSet Profile", value=4,
                                      br(),
                                      textOutput("profile_geneNumtext2"),
                                      br()
                    ),
                    tabPanel("How to Cite",value=5,
                             #tags$div(class="header", checked=NA,
                             h4("GISPA is a product of the Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University (https://bbisr.winship.emory.edu)."),
                             br(),
                             h4("This work is funded by the Leukemia and Lymphoma Society Translational Research Program Award (Jeanne Kowalski); Georgia Research Alliance Scientist Award (Jeanne Kowalski); a Team Science Seed Funding from the Winship Cancer Institute of Emory University (Lawrence H. Boise, Sagar Lonial, Michael R. Rossi); Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University and NIH/NCI [Award number P30CA138292, in part]. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH."),
                             br(),
                             h4("Bioconductor R package for GISPA is available at the https://www.bioconductor.org/packages/GISPA"),
                             br(),
                             h4("Please cite the GISPA method as: Kowalski J, Dwivedi B, Newman S, Switchenko JM, Pauly R, Gutman DA, Arora J, Gandhi K, Ainslie K, Doho G, Qin Z, Moreno CS, Rossi MR, Vertino PM, Lonial S, Bernal-Mizrachi L, Boise LH. Gene integrated set profile analysis: a context-based approach for inferring biological endpoints. Nucleic Acids Res. 2016 Apr 20;44(7):e69. doi: 10.1093/nar/gkv1503. Epub 2016 Jan 29. PubMed PMID: 26826710; PubMed Central PMCID: PMC4838358."),
                             br(),
                             tags$iframe(src = "shinyGISPA_Manual.pdf", width = 1200, height = 1000)
                    )
              )),
              
              conditionalPanel(condition = "input.analysisType == '3d'",     
                 tabsetPanel(id="Threedim",
                     #tabPanel("Manual",
                     #         tags$iframe(src = "shinyGISPA_manual.pdf", width = 1200, height = 1000)
                     #),
                     tabPanel ("Input Data", value=1,
                      conditionalPanel(condition = "input.DataDescription == 10",  
                            fluidRow(
                               column(4,
                                        textOutput("geneNumtext3d_1"),
                                        textOutput("probeNumtext3d_1"),
                                        textOutput("sampNumtext3d_1"),
                                        br(),
                                        textOutput("refSamp3d_1"),
                                        textOutput("sampOne3d_1"),
                                        textOutput("sampTwo3d_1"),
                                        br()
                                      
                               ))),
                      conditionalPanel(condition = "input.DataDescription == 20", 
                               fluidRow(
                                 column(4,
                                        textOutput("geneNumtext3d_2"),
                                        textOutput("probeNumtext3d_2"),
                                        textOutput("sampNumtext3d_2"),
                                        br(),
                                        textOutput("refSamp3d_2"),
                                        textOutput("sampOne3d_2"),
                                        textOutput("sampTwo3d_2"),
                                        br()
                                        
                                 ))),
                      conditionalPanel(condition = "input.DataDescription == 30",   
                               fluidRow(
                                 column(4,
                                        textOutput("geneNumtext3d_3"),
                                        textOutput("probeNumtext3d_3"),
                                        textOutput("sampNumtext3d_3"),
                                        br(),
                                        textOutput("refSamp3d_3"),
                                        textOutput("sampOne3d_3"),
                                        textOutput("sampTwo3d_3"),
                                        br()
                                        
                                 )))),
                   tabPanel("Results Table", value=2,
                            fluidRow(
                              column(12,
                                     dataTableOutput("gispa_data3")
                              ))
                   ),
                   tabPanel("Diagnostic Plot", value=3,
                            fluidRow(
                              column(6,
                                     plotOutput(("gispa_plot3"), height="600px", width="600px")
                              ),
                              column(6,
                                     plotOutput(("slope_plot3"), height="600px", width="700px")
                              )),
                            fluidRow(
                              column(4,
                                     plotOutput(("profile_boxplot3d_1"), height="600px", width="500px")
                              ),
                              column(4,
                                     plotOutput(("profile_boxplot3d_2"), height="600px", width="500px")
                              ),
                              column(4,
                                     plotOutput(("profile_boxplot3d_3"), height="600px", width="500px")
                              ))
                   ),
                   tabPanel("GeneSet Profile", value=4,
                      #conditionalPanel(condition = "input.setProp == 10",  
                            fluidRow(
                              column(8,
                                     br(),
                                     textOutput("profile_geneNumtext3"),
                                     br()
                              ))#)
                   ),
                   tabPanel("How to Cite",value=5,
                            #tags$div(class="header", checked=NA,
                            h4("GISPA is a product of the Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University (https://bbisr.winship.emory.edu)."),
                            br(),
                            h4("This work is funded by the Leukemia and Lymphoma Society Translational Research Program Award (Jeanne Kowalski); Georgia Research Alliance Scientist Award (Jeanne Kowalski); a Team Science Seed Funding from the Winship Cancer Institute of Emory University (Lawrence H. Boise, Sagar Lonial, Michael R. Rossi); Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University and NIH/NCI [Award number P30CA138292, in part]. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH."),
                            br(),
                            h4("Bioconductor R package for GISPA is available at the https://www.bioconductor.org/packages/GISPA"),
                            br(),
                            h4("Please cite the GISPA method as: Kowalski J, Dwivedi B, Newman S, Switchenko JM, Pauly R, Gutman DA, Arora J, Gandhi K, Ainslie K, Doho G, Qin Z, Moreno CS, Rossi MR, Vertino PM, Lonial S, Bernal-Mizrachi L, Boise LH. Gene integrated set profile analysis: a context-based approach for inferring biological endpoints. Nucleic Acids Res. 2016 Apr 20;44(7):e69. doi: 10.1093/nar/gkv1503. Epub 2016 Jan 29. PubMed PMID: 26826710; PubMed Central PMCID: PMC4838358."),
                            br(),
                            tags$iframe(src = "shinyGISPA_Manual.pdf", width = 1200, height = 1000)
                   )
                            
              )),
              conditionalPanel(condition = "input.Threedim == 1 & input.analysisType == '3d'",
                 tabsetPanel(id="DataDescription",
                   tabPanel ("Data 1",value=10,
                             fluidRow(
                               column(4,
                                      plotOutput(("boxPlot3d_1"), height="600px", width="500px")
                               )),
                             fluidRow(
                               column(6,
                                      dataTableOutput("inputData3d_1")
                               ))
                    ),
                   tabPanel ("Data 2",value=20,
                             fluidRow(
                               column(4,
                                      plotOutput(("boxPlot3d_2"), height="600px", width="500px")
                               )),
                             fluidRow(
                               column(6,
                                      dataTableOutput("inputData3d_2")
                               ))
                   ),
                   tabPanel ("Data 3",value=30,
                             fluidRow(
                               column(4,
                                      plotOutput(("boxPlot3d_3"), height="600px", width="500px")
                               )),
                             fluidRow(
                               column(6,
                                      dataTableOutput("inputData3d_3")
                               ))
                   )
              )),
              conditionalPanel(condition = "input.Twodim == 4 & input.analysisType == '2d'",
                 tabsetPanel(id="setProp",
                   tabPanel ("Between Sample Differences",value=42,
                             fluidRow(
                               column(6,
                                      plotOutput(("stacked_barplot2"),height="1000px", width="1000px")
                               ))
                   ),
                   tabPanel ("Between Feature Differences",value=41,
                             fluidRow(
                               column(6,
                                      plotOutput(("stacked_propbarplot2"), height="1000px", width="1000px")
                               ))
                   )
              )),
              conditionalPanel(condition = "input.Threedim == 4 & input.analysisType == '3d'",
                 tabsetPanel(id="geneSetProfiles",
                   tabPanel ("Between Sample Differences",
                             fluidRow(
                               column(6,
                                      plotOutput(("stacked_barplot3"),height="1000px", width="1000px")
                               ))
                   ),
                   tabPanel ("Between Feature Differences",
                             fluidRow(
                               column(6,
                                      plotOutput(("stacked_propbarplot3"), height="1000px", width="1000px")
                               ))
                   )
                 )),
              conditionalPanel(condition = "input.Onedim == 4 & input.analysisType == '1d'",
                 tabsetPanel(id="geneSetProfiles",
                   tabPanel ("Between Sample Differences",
                             fluidRow(
                               column(6,
                                      plotOutput(("stacked_barplot"),height="1000px", width="1000px")
                               ))
                   )
                 ))),
            
              column(1,
                       wellPanel(          
                         h4("Sample Colors"),
                         tags$hr(),
                         colourInput(inputId="Ref", label= "Reference", palette = "limited", value = "red"),
                         colourInput(inputId="Cmp1", label= "Sample 1", palette = "limited", value = "blue"),
                         colourInput(inputId="Cmp2", label= "Sample 2", palette = "limited", value = "green")
                       ),
                       wellPanel(          
                         h4("Data Type Colors"),
                         tags$hr(),
                         colourInput(inputId="dt1", label= "Data Type 1", palette = "limited", value = "#F4A460"),
                         colourInput(inputId="dt2", label= "Data Type 2", palette = "limited", value = "#1E90FF"),
                         colourInput(inputId="dt3", label= "Data Type 3", palette = "limited", value = "#9400D3")
                       )
              )
        )
))#shinyServer

