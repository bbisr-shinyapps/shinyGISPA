#'@name stackedBarplot
#'@aliases stackedBarplot
#'@title A plotting function for each sample proportion
#'@description Given a gene, this function will plot the proportion of each sample over all the samples for each gene or row within each data type
#'@usage stackedBarplot(gispa.output,feature,cpt,type,input.cex,input.cex.lab,input.gap,samp.col,strip.col)
#'@param gispa.output : A data frame containing genes as rows followed by between gene feature profile statistics for each sample.
#'@param feature : Analysis type i.e., one ('1d'), two ('2d') or three ('3d') dimensional feature analysis.
#'@param cpt : Changepoint gene set to be plotted. 
#'@param type : Type of data, e.g., EXP  (default) for expression, VAR of variants, CNV for copy number change. 
#'@param input.cex : character (or symbol expansion): x- and y-axis labels
#'@param input.cex.lab : character (or symbol expansion) for the horizontal and vertial strip labels 
#'@param input.gap : gap or distance between each data type plot
#'@param samp.col : a vector of colors of the bar for the three sample groups, reference followed by other two samples
#'@param strip.col : color to be used for the vertial strip 
#'@details This function expects the output from GISPA function of GISPA package, and highlights the gene set in the selected changepoints and their proportion in each of the three sample groups.
#'@return Barplot illustrating each sample proportion for each gene in the selected change point
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@importFrom HH likert
#'@importFrom lattice panel.abline
#'@importFrom latticeExtra layer
#'@importFrom stats update
#'@examples
#'id <- 20 ## number of rows
#'s <- 4 ## number of columns
#'dm <- matrix(runif(id*s,min=0,max=100), nrow=id, ncol=s, 
#'                  dimnames=list(paste("gene", 1:id, sep=""), 
#'                  paste("sample", 1:s, sep="")))
#'changepoints <- sort(sample(1:2, id, replace=TRUE))
#'dm <- cbind(dm,changepoints)
#'stackedBarplot(gispa.output=dm, feature=1, cpt=1, type="EXP",
#'               input.cex=1.5, input.cex.lab=1.5, input.gap=0.5,
#'               samp.col=c("red", "green", "blue"), strip.col="yellow")
#'@export

stackedBarplot = function(gispa.output, feature=1, cpt=1, type="EXP", 
                          input.cex=0.5, input.cex.lab=0.5, input.gap=0.1, 
                          samp.col=c("grey0", "grey60", "grey90"), 
                          strip.col = c("yellow", "bisque")){
  
      #select geneset data for input changepoint 
      gispa.output <- gispa.output[gispa.output[,ncol(gispa.output)]!=1000,]
      gispa.output <- gispa.output[gispa.output[,ncol(gispa.output)]==cpt,]
      genelist <- rownames(gispa.output)
        
      my.settings <- list(
        strip.background=list(col=strip.col[1]),
        strip.border=list(col="black")
      )
      #data type 1
      per_DF <- data.matrix(gispa.output[,c(3:5,ncol(gispa.output))])
      row.names(per_DF) <- NULL
      #compute percentage contribution from each feature profile statistics
      bgfps <- prop.table(per_DF[,1:3], margin=1)*100
      colnames(bgfps)[1]<-"Ref"
      colnames(bgfps)[2]<-"S1"
      colnames(bgfps)[3]<-"S2"
      bgfps <-data.frame(bgfps, per_DF[,ncol(per_DF)])
      colnames(bgfps)[4]<-"changepoints"
      ##genelist <- gispa.output[,2]
      
      DT1 <- likert(~ Ref + S1 + S2 | changepoints, data=bgfps,
                    scales=list(y=list(at=length(genelist):1, 
                                       relation="free", 
                                       labels=genelist, 
                                       cex=input.cex, font=4), 
                                x=list(relation="same", 
                                       at=seq(0, 100, 50), 
                                       labels=seq(0, 100, 50), 
                                       cex=input.cex)),
                    strip=FALSE, between=list(y=c(0,0,0,0)), 
                    par.strip.text=list(cex=input.cex.lab,font=2,lines=1.5), 
                    xlab=list("Proportion",cex=input.cex.lab), 
                    cex.lab=input.cex.lab, ylab="", main="",
                    par.settings=my.settings, Reference=0, 
                    col=c(samp.col[1], samp.col[2], samp.col[3]))
      
      ALL <- update(cbind(DataType1=as.vector(DT1)))
      ALL <- DT1 + layer(panel.abline(v=50, col="black", lty=2))
      
      #data type 2
      if(feature == 2 | feature == 3){
        per_DF <- data.matrix(gispa.output[,c(7:9,ncol(gispa.output))])
        row.names(per_DF) <- NULL
        bgfps <- prop.table(per_DF[,1:3], margin=1)*100
        colnames(bgfps)[1]<-"Ref"
        colnames(bgfps)[2]<-"S1"
        colnames(bgfps)[3]<-"S2"
        bgfps <-data.frame(bgfps, per_DF[,ncol(per_DF)])
        colnames(bgfps)[4]<-"changepoints"
        DT2 <- likert(~ Ref + S1 + S2 | changepoints, data=bgfps,
                      scales=list(y=list(at=length(genelist):1, 
                                         relation="free", 
                                         labels="", 
                                         cex=1.0, font=4)), 
                      layout=c(1,length(cpt)), strip=FALSE, strip.left=FALSE, 
                      between=list(y=c(0,0,0)), 
                      par.strip.text=list(cex=2.0,lines=1.2), 
                      xlab="Proportion", ylab="", main="", Reference=0, 
                      col=c(samp.col[1], samp.col[2], samp.col[3]))
        
        my.settings <- list(
          strip.background=list(col=c(strip.col[1],strip.col[2])),
          strip.border=list(col="black")
        )
        ALL <- update(cbind(DataType1=as.vector(DT1), 
                            DataType2=as.vector(DT2)),
                      scales=list(y=list(relation="free", 
                                        at = list(length(genelist):1,NULL,
                                                  length(genelist):1,NULL,
                                                  length(genelist):1,NULL,
                                                  length(genelist):1,NULL), 
                                         labels=list(genelist,NULL), 
                                         cex=input.cex, font=4),  
                                  x=list(relation="same", 
                                         at=seq(0, 100, 50), 
                                         labels=seq(0, 100, 50), 
                                         cex=input.cex)),
                      between=list(x=input.gap, y=c(0.1,0.3,0.1)), 
                      par.strip.text=list(cex=input.cex.lab,font=2,lines=1.5),
                      xlab=list("Proportion",cex=input.cex.lab), 
                      par.settings=my.settings, Reference=0, 
                      col=c(samp.col[1], samp.col[2], samp.col[3]))
        ALL <- ALL + layer(panel.abline(v=50, col="black", lty=2), under=TRUE)
        dimnames(ALL)[[2]] <-cpt
        dimnames(ALL)[[1]] <-c(type[1],type[2])
        
      }
      #data type 3
      if (feature == 3){
        per_DF <- data.matrix(gispa.output[,c(11:13,ncol(gispa.output))])
        row.names(per_DF) <- NULL
        bgfps <- prop.table(per_DF[,1:3], margin=1)*100
        colnames(bgfps)[1]<-"Ref"
        colnames(bgfps)[2]<-"S1"
        colnames(bgfps)[3]<-"S2"
        bgfps <-data.frame(bgfps, per_DF[,ncol(per_DF)])
        colnames(bgfps)[4]<-"changepoints"
        DT3 <- likert(~ Ref + S1 + S2 | changepoints, data=bgfps,
                      scales=list(y=list(at=length(genelist):1, 
                                         relation="free", 
                                         labels="", 
                                         cex=1.0, font=4)),
                      layout=c(1,length(cpt)), strip=FALSE, strip.left=FALSE, 
                      between=list(y=c(0,0,0)), 
                      par.strip.text=list(cex=2.0,lines=1.2), 
                      xlab="Proportion", ylab="", main="", Reference=0, 
                      col=c(samp.col[1], samp.col[2], samp.col[3]))

        my.settings <- list(
          strip.background=list(col=c(strip.col[1],strip.col[2])),
          strip.border=list(col="black")
        )
        ALL <- update(cbind(DataType1=as.vector(DT1),
                            DataType2=as.vector(DT2),
                            DataType3=as.vector(DT3)),
                      scales=list(y=list(relation="free", 
                                    at = list(length(genelist):1,NULL,NULL,
                                              length(genelist):1,NULL,NULL,
                                              length(genelist):1,NULL,NULL,
                                              length(genelist):1,NULL,NULL),
                                         labels=list(genelist,NULL,NULL),
                                         cex=input.cex, font=4),  
                                  x=list(relation="same", 
                                         at=seq(0, 100, 50), 
                                         labels=seq(0, 100, 50), 
                                         cex=input.cex)), 
                      between=list(x=input.gap, y=c(0.1,0.3,0.1)), 
                      xlab=list("Proportion",cex=input.cex.lab), 
                      par.strip.text=list(cex=input.cex.lab,font=2,lines=1.5),
                      par.settings=my.settings, Reference=0,
                      col=c(samp.col[1], samp.col[2], samp.col[3]))
        ALL <- ALL + layer(panel.abline(v=50, col="black", lty=2), under=TRUE)
        dimnames(ALL)[[2]] <-cpt
        dimnames(ALL)[[1]] <-c(type[1],type[2],type[3])
    }
  
  return(ALL)
}
