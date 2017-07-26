options(shiny.maxRequestSize=10000*1024^2)
options(bitmapType='cairo')


library(changepoint)
library(colourpicker)
library(data.table)
library(plyr)
library(genefilter)
library(ggplot2)
library(graphics)
library(HH)
library(lattice)
library(latticeExtra)
library(scatterplot3d)
library(shiny)
library(splitstackshape)
library(stats)
library(Biobase)
library(GSEABase)

source("functions/computePS.R")
source("functions/cptModel.R")
source("functions/cptPlot.R")
source("functions/cptSlopeplot.R")
source("functions/GISPA.R")
source("functions/getBoxPlot.R")
source("functions/stackedBarplot.R")
source("functions/propBarplot.R")

shinyServer(function(input, output, session) {
  
  # read in initial data input.csv i.e. a file containing the data
  # this file contains column corrosponding to samples
  # this file contains rows corrosponding to genes
  # the data directory of the app folder contains those .csv files
  
  ##########Read uploaded data!
  ##Input data for 1D
  getData = reactive({
    
    if(input$analysisType == '1d'){
      
      if (input$file == 'ED'){
        dataFile <- "data/MM_CL_exp.txt"
        dat <- read.csv(dataFile, sep="\t", header=T)
      } else if (input$file == 'LD'){
          if (is.null(input$dataset)) return(NULL)
          else{
            dataFile <- input$dataset$datapath
            dat <- read.csv(dataFile, sep="\t", header=T)
          }
        }else{
        return(NULL)
        }
      
      # phenoData:Sample information
      sampleid = data.frame(colnames(dat[,c(input$ref_1d, input$cmp1_1d, input$cmp2_1d)]))
      sampleid$treatment = c(input$ref_1d, input$cmp1_1d, input$cmp2_1d)
      colnames(sampleid) <- c("id", "treatment")
      rownames(sampleid) <- sampleid$id
      pdata <- AnnotatedDataFrame(sampleid)
      # featureData: gene information
      dat_dedup <- dat[!duplicated(dat[2]),]
      dat_dedup <- dat_dedup[,c(1,2,input$ref_1d, input$cmp1_1d, input$cmp2_1d)]
      genelist <- dat_dedup[,c(2,1)]
      rownames(genelist) <- dat_dedup[,2]
      fdata <- AnnotatedDataFrame(genelist)
      # expression daya
      m <- as.matrix(dat_dedup[,c(3,4,5)])
      rownames(m) <- dat_dedup[,2]
      ## create ExpressionSet object:
      eset <- new("ExpressionSet", exprs = m, phenoData = pdata, featureData = fdata)
      samplegrp = colnames(dat[,c(input$ref_1d,input$cmp1_1d,input$cmp2_1d)])
      uniqGenes = dat[!duplicated(dat[1]),]
      uniqProbes = dat[!duplicated(dat[2]),]
      data = list(ori=dat_dedup, obj=eset, geneid=nrow(uniqGenes), probeid=nrow(uniqProbes), samplegrp=samplegrp,  sampleid=nrow(sampleid))
      return(data)
    }
  })
  
  ##Input data for 2D
  getData2 = reactive({
    
    if(input$analysisType == '2d'){
      
      if (input$file1 == 'ED' && input$file2 == 'ED'){
        dataFile1 <- "data/MM_CL_exp.txt"
        dat1 <- read.csv(dataFile1, sep="\t", header=T)
        dataFile2 <- "data/MM_CL_var.txt"
        dat2 <- read.csv(dataFile2, sep="\t", header=T)
      } else if (input$file1 == 'LD'){
          if (is.null(input$dataset1)) return(NULL)
          else{
              dataFile1 <- input$dataset1$datapath
              dat1 <- read.csv(dataFile1, sep="\t", header=T)
          }
          if (is.null(input$dataset2)) return(NULL)
          else{
              dataFile2 <- input$dataset2$datapath
              dat2 <- read.csv(dataFile2, sep="\t", header=T)
          }
      } else{
        return(NULL)
      }
        
      # phenoData:Sample information
      sampleid1 = data.frame(colnames(dat1[,c(input$ref_2d, input$cmp1_2d, input$cmp2_2d)]))
      sampleid1$treatment = c(input$ref_2d, input$cmp1_2d, input$cmp2_2d)
      colnames(sampleid1) <- c("id", "treatment")
      rownames(sampleid1) <- sampleid1$id
      pdata1 <- AnnotatedDataFrame(sampleid1)
      # featureData: gene information
      dat_dedup1 <- dat1[!duplicated(dat1[2]),]
      dat_dedup1 <- dat_dedup1[,c(1,2,input$ref_2d, input$cmp1_2d, input$cmp2_2d)]
      genelist1 <- dat_dedup1[,c(2,1)]
      rownames(genelist1) <- dat_dedup1[,2]
      fdata1 <- AnnotatedDataFrame(genelist1)
      # expression daya
      m1 <- as.matrix(dat_dedup1[,c(3,4,5)])
      rownames(m1) <- dat_dedup1[,2]
      ## create ExpressionSet object:
      eset1 <- new("ExpressionSet", exprs = m1, phenoData = pdata1, featureData = fdata1)
      samplegrp1 = colnames(dat1[,c(input$ref_2d,input$cmp1_2d,input$cmp2_2d)])
      uniqGenes1 = dat1[!duplicated(dat1[1]),]
      uniqProbes1 = dat1[!duplicated(dat1[2]),]
 
      # phenoData:Sample information
      sampleid2 = data.frame(colnames(dat2[,c(input$ref_2d, input$cmp1_2d, input$cmp2_2d)]))
      sampleid2$treatment = c(input$ref_2d, input$cmp1_2d, input$cmp2_2d)
      colnames(sampleid2) <- c("id", "treatment")
      rownames(sampleid2) <- sampleid2$id
      pdata2 <- AnnotatedDataFrame(sampleid2)
      # featureData: gene information
      dat_dedup2 <- dat2[!duplicated(dat2[2]),]
      dat_dedup2 <- dat_dedup2[,c(1,2,input$ref_2d, input$cmp1_2d, input$cmp2_2d)]
      genelist2 <- dat_dedup2[,c(2,1)]
      rownames(genelist2) <- dat_dedup2[,2]
      fdata2 <- AnnotatedDataFrame(genelist2)
      # expression daya
      m2 <- as.matrix(dat_dedup2[,c(3,4,5)])
      rownames(m2) <- dat_dedup2[,2]
      ## create ExpressionSet object:
      eset2 <- new("ExpressionSet", exprs = m2, phenoData = pdata2, featureData = fdata2)
      samplegrp2 = colnames(dat2[,c(input$ref_2d,input$cmp1_2d,input$cmp2_2d)])
      uniqGenes2 = dat2[!duplicated(dat2[1]),]
      uniqProbes2 = dat2[!duplicated(dat2[2]),]
      
      data = list(ori1=dat_dedup1, obj1=eset1, geneid1=nrow(uniqGenes1), probeid1=nrow(uniqProbes1), samplegrp1=samplegrp1,  sampleid1=nrow(sampleid1),
                  ori2=dat_dedup2, obj2=eset2, geneid2=nrow(uniqGenes2), probeid2=nrow(uniqProbes2), samplegrp2=samplegrp2,  sampleid2=nrow(sampleid2))
    }
    
  })  
  
  ##Input data for 3D
  getData3 = reactive({
    
    if(input$analysisType == '3d'){
      
      if (input$file1 == 'ED' && input$file2 == 'ED' && input$file3 == 'ED'){
        dataFile1 <- "data/MM_CL_exp.txt"
        dat1 <- read.csv(dataFile1, sep="\t", header=T)
        dataFile2 <- "data/MM_CL_var.txt"
        dat2 <- read.csv(dataFile2, sep="\t", header=T)
        dataFile3 <- "data/MM_450K_met_gw.txt"
        dat3 <- read.csv(dataFile3, sep="\t", header=T)
      } else if (input$file1 == 'LD' && input$file2 == 'LD' && input$file3 == 'LD'){
          if (is.null(input$dataset1)) return(NULL)
          else{
            dataFile1 <- input$dataset1$datapath
            dat1 <- read.csv(dataFile1, sep="\t", header=T)
          }
          if (is.null(input$dataset2)) return(NULL)
          else{
            dataFile2 <- input$dataset2$datapath
            dat2 <- read.csv(dataFile2, sep="\t", header=T)
          }
          if (is.null(input$dataset3)) return(NULL)
          else{
            dataFile3 <- input$dataset3$datapath
            dat3 <- read.csv(dataFile3, sep="\t", header=T)
          }
      }else{
        # if no data file is uploaded, return NULL
        return(NULL)
      }
      
      # phenoData:Sample information
      sampleid1 = data.frame(colnames(dat1[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)]))
      sampleid1$treatment = c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)
      colnames(sampleid1) <- c("id", "treatment")
      rownames(sampleid1) <- sampleid1$id
      pdata1 <- AnnotatedDataFrame(sampleid1)
      # featureData: gene information
      dat_dedup1 <- dat1[!duplicated(dat1[2]),]
      dat_dedup1 <- dat_dedup1[,c(1,2,input$ref_3d, input$cmp1_3d, input$cmp2_3d)]
      genelist1 <- dat_dedup1[,c(2,1)]
      rownames(genelist1) <- dat_dedup1[,2]
      fdata1 <- AnnotatedDataFrame(genelist1)
      # expression daya
      m1 <- as.matrix(dat_dedup1[,c(3,4,5)])
      rownames(m1) <- dat_dedup1[,2]
      ## create ExpressionSet object:
      eset1 <- new("ExpressionSet", exprs = m1, phenoData = pdata1, featureData = fdata1)
      samplegrp1 = colnames(dat1[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)])
      uniqGenes1 = dat1[!duplicated(dat1[1]),]
      uniqProbes1 = dat1[!duplicated(dat1[2]),]
      
      # phenoData:Sample information
      sampleid2 = data.frame(colnames(dat2[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)]))
      sampleid2$treatment = c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)
      colnames(sampleid2) <- c("id", "treatment")
      rownames(sampleid2) <- sampleid2$id
      pdata2 <- AnnotatedDataFrame(sampleid2)
      # featureData: gene information
      dat_dedup2 <- dat2[!duplicated(dat2[2]),]
      dat_dedup2 <- dat_dedup2[,c(1,2,input$ref_3d, input$cmp1_3d, input$cmp2_3d)]
      genelist2 <- dat_dedup2[,c(2,1)]
      rownames(genelist2) <- dat_dedup2[,2]
      fdata2 <- AnnotatedDataFrame(genelist2)
      # expression daya
      m2 <- as.matrix(dat_dedup2[,c(3,4,5)])
      rownames(m2) <- dat_dedup2[,2]
      ## create ExpressionSet object:
      eset2 <- new("ExpressionSet", exprs = m2, phenoData = pdata2, featureData = fdata2)
      samplegrp2 = colnames(dat2[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)])
      uniqGenes2 = dat2[!duplicated(dat2[1]),]
      uniqProbes2 = dat2[!duplicated(dat2[2]),]
      
      # phenoData:Sample information
      sampleid3 = data.frame(colnames(dat3[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)]))
      sampleid3$treatment = c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)
      colnames(sampleid3) <- c("id", "treatment")
      rownames(sampleid3) <- sampleid3$id
      pdata3 <- AnnotatedDataFrame(sampleid3)
      # featureData: gene information
      dat_dedup3 <- dat3[!duplicated(dat3[2]),]
      dat_dedup3 <- dat_dedup3[,c(1,2,input$ref_3d, input$cmp1_3d, input$cmp2_3d)]
      genelist3 <- dat_dedup3[,c(2,1)]
      rownames(genelist3) <- dat_dedup3[,2]
      fdata3 <- AnnotatedDataFrame(genelist3)
      # expression daya
      m3 <- as.matrix(dat_dedup3[,c(3,4,5)])
      rownames(m3) <- dat_dedup3[,2]
      ## create ExpressionSet object:
      eset3 <- new("ExpressionSet", exprs = m3, phenoData = pdata3, featureData = fdata3)
      samplegrp3 = colnames(dat3[,c(input$ref_3d, input$cmp1_3d, input$cmp2_3d)])
      uniqGenes3 = dat3[!duplicated(dat3[1]),]
      uniqProbes3 = dat3[!duplicated(dat3[2]),]
      
      data = list(ori1=dat_dedup1, obj1=eset1, geneid1=nrow(uniqGenes1), probeid1=nrow(uniqProbes1), samplegrp1=samplegrp1,  sampleid1=nrow(sampleid1),
                  ori2=dat_dedup2, obj2=eset2, geneid2=nrow(uniqGenes2), probeid2=nrow(uniqProbes2), samplegrp2=samplegrp2,  sampleid2=nrow(sampleid2),
                  ori3=dat_dedup3, obj3=eset3, geneid3=nrow(uniqGenes3), probeid3=nrow(uniqProbes3), samplegrp3=samplegrp3,  sampleid3=nrow(sampleid3))
      
       }
    
  })  
  
  
  ########### Results Table Output
  getGISPA_results = reactive({
    
    if(input$analysisType == '1d'){
      if (is.null(getData()$ori)) return(NULL)
      else{
        withProgress(message = 'Running GISPA', value = 0.1, {
          Sys.sleep(0.25)
          
          # Create 0-row data frame which will be used to store data
          dat <- data.frame(x = numeric(0), y = numeric(0))
          
          # withProgress calls can be nested, in which case the nested text appears
          # below, and a second bar is shown.
          withProgress(message = 'Generating data', detail = "part 0", value = 0, {
            for (i in 1:10) {
              # Each time through the loop, add another row of data. This a stand-in
              # for a long-running computation.
              dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
              
              # Increment the progress bar, and update the detail text.
              incProgress(0.1, detail = paste("part", i))
              
              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
            }
            getresults <- GISPA(feature=1, f.sets=c(getData()$obj),f.profiles=input$profile,
                                  ref.samp.idx=1, comp.samp.idx=c(2,3),
                                  cpt.data=input$changes, cpt.method=input$method, cpt.max=input$max)
            getresults$cpt_out <- getresults$cpt_out[,-c(1,7:9,11)]
            colnames(getresults$cpt_out)[1] <- colnames(getData()$ori)[2]
            colnames(getresults$cpt_out)[2] <- colnames(getData()$ori)[1]
            colnames(getresults$cpt_out)[3] <- colnames(getData()$ori)[input$ref_1d]
            colnames(getresults$cpt_out)[4] <- colnames(getData()$ori)[input$cmp1_1d]
            colnames(getresults$cpt_out)[5] <- colnames(getData()$ori)[input$cmp2_1d]
            dat <- list(cpt_out=getresults$cpt_out, cpt_plot=getresults$cpt_plot, all_cutoffs=getresults$all_cutoffs)
          })
        })
      }
    }
    else if (input$analysisType == '2d'){
      if (is.null(getData2()$ori1) || is.null(getData2()$ori2)) return(NULL)
      else{
        withProgress(message = 'Running GISPA', value = 0.1, {
          Sys.sleep(0.25)
          
          # Create 0-row data frame which will be used to store data
          dat <- data.frame(x = numeric(0), y = numeric(0))
          
          # withProgress calls can be nested, in which case the nested text appears
          # below, and a second bar is shown.
          withProgress(message = 'Generating data', detail = "part 0", value = 0, {
            for (i in 1:10) {
              # Each time through the loop, add another row of data. This a stand-in
              # for a long-running computation.
              dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
              
              # Increment the progress bar, and update the detail text.
              incProgress(0.1, detail = paste("part", i))
              
              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
            }
            getresults <- GISPA(feature=2, f.sets=c(getData2()$obj1, getData2()$obj2),
                                f.profiles=c(input$profile1,input$profile2),
                                ref.samp.idx=1, comp.samp.idx=c(2,3),
                                cpt.data=input$changes, cpt.method=input$method, cpt.max=input$max)
            
            getresults$cpt_out <- getresults$cpt_out[,-c(1,7,12:15,17)]
            colnames(getresults$cpt_out)[1] <- colnames(getData2()$ori1)[1]
            colnames(getresults$cpt_out)[2] <- colnames(getData2()$ori1)[2]
            colnames(getresults$cpt_out)[3] <- colnames(getData2()$ori1)[input$ref_2d]
            colnames(getresults$cpt_out)[4] <- colnames(getData2()$ori1)[input$cmp1_2d]
            colnames(getresults$cpt_out)[5] <- colnames(getData2()$ori1)[input$cmp2_2d]
            colnames(getresults$cpt_out)[6] <- colnames(getData2()$ori2)[2]
            colnames(getresults$cpt_out)[7] <- colnames(getData2()$ori2)[input$ref_2d]
            colnames(getresults$cpt_out)[8] <- colnames(getData2()$ori2)[input$cmp1_2d]
            colnames(getresults$cpt_out)[9] <- colnames(getData2()$ori2)[input$cmp2_2d]
            
            dat <- list(cpt_out=getresults$cpt_out, cpt_plot=getresults$cpt_plot, all_cutoffs=getresults$all_cutoffs)
          })
        })
      }
    }
    else if (input$analysisType == '3d'){
      if (is.null(getData3()$ori1) || is.null(getData3()$ori2) || is.null(getData3()$ori3)) return(NULL)
      else{
        withProgress(message = 'Running GISPA', value = 0.1, {
          Sys.sleep(0.25)
          
          # Create 0-row data frame which will be used to store data
          dat <- data.frame(x = numeric(0), y = numeric(0))
          
          # withProgress calls can be nested, in which case the nested text appears
          # below, and a second bar is shown.
          withProgress(message = 'Generating data', detail = "part 0", value = 0, {
            for (i in 1:10) {
              # Each time through the loop, add another row of data. This a stand-in
              # for a long-running computation.
              dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
              
              # Increment the progress bar, and update the detail text.
              incProgress(0.1, detail = paste("part", i))
              
              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
            }
            getresults <- GISPA(feature=3, f.sets=c(getData3()$obj1, getData3()$obj2, getData3()$obj3),
                                f.profiles=c(input$profile1,input$profile2, input$profile3),
                                ref.samp.idx=1, comp.samp.idx=c(2,3),
                                cpt.data=input$changes, cpt.method=input$method, cpt.max=input$max)
            
            getresults$cpt_out <- getresults$cpt_out[,-c(1,7,12,17:21,23)]
            colnames(getresults$cpt_out)[1] <- colnames(getData3()$ori1)[1]
            colnames(getresults$cpt_out)[2] <- colnames(getData3()$ori1)[2]
            colnames(getresults$cpt_out)[3] <- colnames(getData3()$ori1)[input$ref_3d]
            colnames(getresults$cpt_out)[4] <- colnames(getData3()$ori1)[input$cmp1_3d]
            colnames(getresults$cpt_out)[5] <- colnames(getData3()$ori1)[input$cmp2_3d]
            colnames(getresults$cpt_out)[6] <- colnames(getData3()$ori2)[2]
            colnames(getresults$cpt_out)[7] <- colnames(getData3()$ori2)[input$ref_3d]
            colnames(getresults$cpt_out)[8] <- colnames(getData3()$ori2)[input$cmp1_3d]
            colnames(getresults$cpt_out)[9] <- colnames(getData3()$ori2)[input$cmp2_3d]
            colnames(getresults$cpt_out)[10] <- colnames(getData3()$ori3)[2]
            colnames(getresults$cpt_out)[11] <- colnames(getData3()$ori3)[input$ref_3d]
            colnames(getresults$cpt_out)[12] <- colnames(getData3()$ori3)[input$cmp1_3d]
            colnames(getresults$cpt_out)[13] <- colnames(getData3()$ori3)[input$cmp2_3d]
            
            dat <- list(cpt_out=getresults$cpt_out, cpt_plot=getresults$cpt_plot, all_cutoffs=getresults$all_cutoffs)
          })
        })
      }
    }
    
    return (dat)
  })
  
  
  
  ########### Data Summary Panel Output for 1D
  output$boxPlot <- renderPlot({
    if (is.null(getData()$ori)) return(NULL)
    else{
      boxplot(getData()$ori[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1,input$Cmp2), ylab=input$dataType, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData <- renderDataTable({
    if (is.null(getData()$ori)) return(NULL)
    else{
      getData()$ori
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData()$geneid, sep = "")
    }
  })
  
  output$probeNumtext <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData()$probeid, sep = "")
    }
  })
  
  output$sampNumtext <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("Number of Samples: ", getData()$sampleid, sep = "")
    }
  })
  
  output$refSamp <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("'Reference' is ", getData()$samplegrp[1], sep = "")
    }
  })
  
  output$sampOne <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("'Sample 1' is ", getData()$samplegrp[2], sep = "")
    }
  })
  
  output$sampTwo <- renderText({ 
    if (is.null(getData()$ori)) return(NULL)
    else{
      paste("'Sample 2' is ", getData()$samplegrp[3], sep = "")
    }
  })
  
  ########### Data Summary Panel Output for 2D
  output$boxPlot1 <- renderPlot({
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      #default 3,4,5 are the data points of first data type
      boxplot(getData2()$ori1[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1,input$Cmp2), ylab=input$dataType1, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData1 <- renderDataTable({
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      getData2()$ori1
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData2()$geneid1, sep = "")
    }
  })
  
  output$probeNumtext1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData2()$probeid1, sep = "")
    }
  })
  
  output$sampNumtext1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      #default 9,10,11 are the data points of second data type
      paste("Number of Samples: ", getData2()$sampleid1, sep = "")
    }
  })
  
  output$refSamp1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      paste("'Reference' is ", getData2()$samplegrp1[1], sep = "")
    }
  })
  
  output$sampOne1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      paste("'Sample 1' is ", getData2()$samplegrp1[2], sep = "")
    }
  })
  
  output$sampTwo1 <- renderText({ 
    if (is.null(getData2()$ori1)) return(NULL)
    else{
      paste("'Sample 2' is ", getData2()$samplegrp1[3], sep = "")
    }
  })
  
  output$boxPlot2 <- renderPlot({
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      boxplot(getData2()$ori2[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1,input$Cmp2), ylab=input$dataType2, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData2 <- renderDataTable({
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      getData2()$ori2
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData2()$geneid2, sep = "")
    }
  })
  
  output$probeNumtext2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData2()$probeid2, sep = "")
    }
  })
  
  output$sampNumtext2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("Number of Samples: ", getData2()$sampleid2, sep = "")
    }
  })
  
  output$refSamp2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("'Reference' is ", getData2()$samplegrp2[1], sep = "")
    }
  })
  
  output$sampOne2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("'Sample 1' is ", getData2()$samplegrp2[2], sep = "")
    }
  })
  
  output$sampTwo2 <- renderText({ 
    if (is.null(getData2()$ori2)) return(NULL)
    else{
      paste("'Sample 2' is ", getData2()$samplegrp2[3], sep = "")
    }
  })
  
  ########### Data Summary Panel Output for 3D
  output$boxPlot3d_1 <- renderPlot({
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      #default 3,4,5 are the data points of first data type
      boxplot(getData3()$ori1[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1, input$Cmp2), ylab=input$dataType1, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData3d_1 <- renderDataTable({
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      getData3()$ori1
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData3()$geneid1, sep = "")
    }
  })
  
  output$probeNumtext3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData3()$probeid1, sep = "")
    }
  })
  
  output$sampNumtext3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      #default 9,10,11 are the data points of second data type
      paste("Number of Samples: ", getData3()$sampleid1, sep = "")
    }
  })
  
  output$refSamp3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      paste("'Reference' is ", getData3()$samplegrp1[1], sep = "")
    }
  })
  
  output$sampOne3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      paste("'Sample 1' is ", getData3()$samplegrp1[2], sep = "")
    }
  })
  
  output$sampTwo3d_1 <- renderText({ 
    if (is.null(getData3()$ori1)) return(NULL)
    else{
      paste("'Sample 2' is ", getData3()$samplegrp1[3], sep = "")
    }
  })
  
  #~~~~~
  output$boxPlot3d_2 <- renderPlot({
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      boxplot(getData3()$ori2[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1,input$Cmp2), ylab=input$dataType2, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData3d_2 <- renderDataTable({
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      getData3()$ori2
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData3()$geneid2, sep = "")
    }
  })
  
  output$probeNumtext3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData3()$probeid2, sep = "")
    }
  })
  
  output$sampNumtext3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("Number of Samples: ", getData3()$sampleid2, sep = "")
    }
  })
  
  output$refSamp3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("'Reference' is ", getData3()$samplegrp2[1], sep = "")
    }
  })
  
  output$sampOne3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("'Sample 1' is ", getData3()$samplegrp2[2], sep = "")
    }
  })
  
  output$sampTwo3d_2 <- renderText({ 
    if (is.null(getData3()$ori2)) return(NULL)
    else{
      paste("'Sample 2' is ", getData3()$samplegrp2[3], sep = "")
    }
  })
  
  
  output$boxPlot3d_3 <- renderPlot({
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      boxplot(getData3()$ori3[,c(3,4,5)], main = "", notch = FALSE, col = c(input$Ref,input$Cmp1,input$Cmp2), ylab=input$dataType3, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
    }       
  })
  
  output$inputData3d_3 <- renderDataTable({
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      getData3()$ori3
    }
  },
  options=list(scrollX = TRUE))
  
  output$geneNumtext3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("Number of Unique Genes: ", getData3()$geneid3, sep = "")
    }
  })
  
  output$probeNumtext3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("Number of Unique Probes: ", getData3()$probeid3, sep = "")
    }
  })
  
  output$sampNumtext3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("Number of Samples: ", getData3()$sampleid3, sep = "")
    }
  })
  
  output$refSamp3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("'Reference' is ", getData3()$samplegrp3[1], sep = "")
    }
  })
  
  output$sampOne3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("'Sample 1' is ", getData3()$samplegrp3[2], sep = "")
    }
  })
  
  output$sampTwo3d_3 <- renderText({ 
    if (is.null(getData3()$ori3)) return(NULL)
    else{
      paste("'Sample 2' is ", getData3()$samplegrp3[3], sep = "")
    }
  })
  
  ## Obtain GISPA results
  ## 1D
  output$gispa_data <- renderDataTable({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getGISPA_results()$cpt_out
    }
  },
  options=list(scrollX = TRUE))
  
  output$gispa_plot <- renderPlot({
    if (is.null(getGISPA_results()$cpt_plot)) return(NULL)
    else{
      cptPlot(getGISPA_results()$cpt_plot, 
              getGISPA_results()$all_cutoffs,
              cpt=as.numeric(input$cpt_to_display))
    }
  })
  output$slope_plot <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      cptSlopeplot(gispa.output=getGISPA_results()$cpt_out, 
                   feature=1,
                   type=input$dataType,
                   cpt=as.numeric(input$cpt_to_display))
    }
  })
  
  ## 2D
  output$gispa_data2 <- renderDataTable({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getGISPA_results()$cpt_out
    }
  },
  options=list(scrollX = TRUE))
  
  output$gispa_plot2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_plot)) return(NULL)
    else{
      cptPlot(getGISPA_results()$cpt_plot, 
              getGISPA_results()$all_cutoffs,
              cpt=as.numeric(input$cpt_to_display))
    }
  })
  output$slope_plot2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      cptSlopeplot(gispa.output=getGISPA_results()$cpt_out, 
                   feature=2,
                   type=c(input$dataType1, input$dataType2),
                   cpt=as.numeric(input$cpt_to_display))
    }
  })
  
  ## 3D
  output$gispa_data3 <- renderDataTable({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getGISPA_results()$cpt_out
    }
  },
  options=list(scrollX = TRUE))
  
  output$gispa_plot3 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_plot)) return(NULL)
    else{
      cptPlot(getGISPA_results()$cpt_plot, 
              getGISPA_results()$all_cutoffs,
              cpt=as.numeric(input$cpt_to_display))
    }
  })
  output$slope_plot3 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      cptSlopeplot(gispa.output=getGISPA_results()$cpt_out, 
                   feature=3,
                   type=c(input$dataType1, input$dataType2, input$dataType3),
                   cpt=as.numeric(input$cpt_to_display))
    }
  })
  
  ## Gene set profile plots
  ## 1D
  output$profile_boxplot <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,as.numeric(input$cpt_to_display),
                 3,4,5,
                 input$dataType,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  output$profile_geneNumtext <- renderText({ 
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      DF_select <- subset(getGISPA_results()$cpt_out,getGISPA_results()$cpt_out$changepoints==as.numeric(input$cpt_to_display))
      paste("Number of genes in the selected changepoint: ", nrow(DF_select), sep = "")
    }
  })
  output$stacked_barplot <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      stackedBarplot(gispa.output=getGISPA_results()$cpt_out,
                     feature=1,
                     cpt=as.numeric(input$cpt_to_display),
                     type=input$dataType,
                     input.cex=input$opt.cex,
                     input.cex.lab=input$opt.cexaxis,
                     input.gap=input$opt.gap,
                     strip.col=c("yellow","bisque"),
                     samp.col=c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  ## 2D
  output$profile_boxplot1 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,
                 as.numeric(input$cpt_to_display),
                 3,4,5,
                 input$dataType1,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  output$profile_boxplot2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,
                 as.numeric(input$cpt_to_display),
                 7,8,9,
                 input$dataType2,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  output$profile_geneNumtext2 <- renderText({ 
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      DF_select <- subset(getGISPA_results()$cpt_out,getGISPA_results()$cpt_out$changepoints==as.numeric(input$cpt_to_display))
      paste("Number of genes in the selected changepoint: ", nrow(DF_select), sep = "")
    }
  })
  output$stacked_barplot2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      stackedBarplot(gispa.output=getGISPA_results()$cpt_out,
                     feature=2,
                     cpt=as.numeric(input$cpt_to_display),
                     type=c(input$dataType1,input$dataType2),
                     input.cex=input$opt.cex,
                     input.cex.lab=input$opt.cexaxis,
                     input.gap=input$opt.gap,
                     strip.col=c("yellow","bisque"),
                     samp.col=c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  output$stacked_propbarplot2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_plot)) return(NULL)
    else{
      propBarplot(getGISPA_results()$cpt_plot,
                  feature=2,
                  as.numeric(input$cpt_to_display),
                  input$opt.cex,
                  input$opt.cexaxis,
                  c(input$dt1,input$dt2),
                  strip.col="yellow")
    }
  })
  
  ## 3D
  output$profile_boxplot3d_1 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,
                 as.numeric(input$cpt_to_display),
                 input$ref_3d,input$cmp1_3d,input$cmp2_3d,
                 input$dataType1,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  output$profile_boxplot3d_2 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,
                 as.numeric(input$cpt_to_display),
                 7,8,9,
                 input$dataType2,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  output$profile_boxplot3d_3 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      getBoxPlot(getGISPA_results()$cpt_out,
                 as.numeric(input$cpt_to_display),
                 11,12,13,
                 input$dataType3,
                 c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  output$profile_geneNumtext3 <- renderText({ 
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      DF_select <- subset(getGISPA_results()$cpt_out,getGISPA_results()$cpt_out$changepoints==as.numeric(input$cpt_to_display))
      paste("Number of genes in the selected changepoint: ", nrow(DF_select), sep = "")
    }
  })
  
  output$stacked_barplot3 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_out)) return(NULL)
    else{
      stackedBarplot(gispa.output=getGISPA_results()$cpt_out,
                     feature=3,
                     cpt=as.numeric(input$cpt_to_display),
                     type=c(input$dataType1,input$dataType2,input$dataType3),
                     input.cex=input$opt.cex,
                     input.cex.lab=input$opt.cexaxis,
                     input.gap=input$opt.gap,
                     strip.col=c("yellow","bisque"),
                     samp.col=c(input$Ref,input$Cmp1,input$Cmp2))
    }
  })
  
  output$stacked_propbarplot3 <- renderPlot({
    if (is.null(getGISPA_results()$cpt_plot)) return(NULL)
    else{
      propBarplot(getGISPA_results()$cpt_plot,
                  feature=3,
                  as.numeric(input$cpt_to_display),
                  input$opt.cex,
                  input$opt.cexaxis,
                  c(input$dt1,input$dt2,input$dt3),
                  strip.col="yellow")
    }
  })
  
  output$downloadGO <- downloadHandler(
    filename = function() { paste("GISPA_results", '.csv', sep='') },
    content = function(file) {
      write.csv(getGISPA_results()$cpt_out, file,row.names=FALSE)
    }
  )
  
  
})



