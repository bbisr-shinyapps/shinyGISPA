#'@name propBarplot
#'@aliases propBarplot
#'@title A plotting function for proportion by sample
#'@description Given a gene, this function will plot the proportion of each sample over the three samples within each data type
#'@usage propBarplot(gispa.output,feature,cpt,input.cex,input.cex.lab,ft.col,strip.col)
#'@param gispa.output : A data matrix of between gene feature profile statistics for each feature with corrosponding identified changepoints. The row names should corrospond to genes or names to be displayed on y-axis
#'@param feature : Analysis type i.e., one ('1'), two ('2') or three ('3') dimensional feature analysis.
#'@param cpt : Changepoints to be plotted. 
#'@param input.cex : character (or symbol expansion) for the x- and y-axis labels
#'@param input.cex.lab : character (or symbol expansion) for the horizontal and vertial strip labels 
#'@param ft.col : a vector of colors of the bar for the features or data types
#'@param strip.col : color to be used for the vertial strip 
#'@details This function expects the output from the main function of GISPA package, and highlights the gene set in the selected changepoints and their proportion in each of the three sample groups by data type.
#'@return Barplot illustrating each sample proportion for each gene in the selected change point
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@importFrom data.table data.table
#'@importFrom HH likert
#'@importFrom lattice panel.abline
#'@importFrom lattice strip.custom
#'@importFrom latticeExtra layer
#'@examples
#'id <- 20 ## number of rows
#'s <- 3 ## number of columns
#'dm <- matrix(runif(id*s,0,10), nrow=id, ncol=s, 
#'                  dimnames=list(paste("gene", 1:id, sep=""), 
#'                  paste("fs", 1:s, sep="")))
#'changepoints <- sort(sample(1:2, id, replace=TRUE))
#'dm <- cbind(dm,changepoints)
#'propBarplot(gispa.output=dm,feature=2,cpt=1,
#'            input.cex=0.5,input.cex.lab=0.5,
#'            ft.col=c("grey0", "grey60"),strip.col="yellow")
#'@export

propBarplot = function(gispa.output, feature=1, cpt=1, 
                       input.cex=1.5, input.cex.lab=1.5, 
                       ft.col=c("grey0","grey60","grey90"),strip.col="yellow"){
  
  gispa.output <- gispa.output[gispa.output[,ncol(gispa.output)]!=1000,]
  gispa.output <- gispa.output[gispa.output[,ncol(gispa.output)]==cpt,]
              
  if(feature==2){
        prop_one = 1/gispa.output[,1]
        prop_two = 1/gispa.output[,2]
        total = prop_one + prop_two
        Pone = (prop_one/total)*100
        Ptwo = (prop_two/total)*100
        data_to_plot <- data.table(rownames(gispa.output), Pone, Ptwo, 
                                   gispa.output[,ncol(gispa.output)])
        colnames(data_to_plot) <- c("gene", "DT1", "DT2", "changepoints")
        
        bp <- likert(~ DT1 + DT2 | changepoints,  data=data_to_plot, 
                     scales=list(y=list(at=length(data_to_plot$gene):1, 
                                        relation="free", 
                                        labels=data_to_plot$gene, 
                                        cex=input.cex, font=4), 
                                 x=list(relation="same", 
                                        at=seq(0, 100, 50), 
                                        labels=seq(0, 100, 50), 
                                        cex=input.cex)),
                     strip=FALSE, strip.left=strip.custom(bg=strip.col), 
                     par.strip.text=list(cex=input.cex.lab,lines=1.5),
                     between=list(y=c(0,0,0,0)), 
                     xlab=list("Proportion",cex=input.cex.lab), 
                     cex.lab=input.cex.lab, ylab="",  main="",
                     Reference=0, auto.key = list(space ="right",columns=1),
                     col=c(ft.col[1], ft.col[2]))
        bp <- bp + layer(panel.abline(v=50, col="black", lty=2))
  }
  
  if(feature==3){
    prop_one = 1/gispa.output[,1]
    prop_two = 1/gispa.output[,2]
    prop_three = 1/gispa.output[,3]
    total = prop_one + prop_two + prop_three
    Pone = (prop_one/total)*100
    Ptwo = (prop_two/total)*100
    Pthree = (prop_three/total)*100
    data_to_plot <- data.frame(rownames(gispa.output), Pone, Ptwo, Pthree, 
                               gispa.output[,ncol(gispa.output)])
    colnames(data_to_plot) <- c("gene", "DT1", "DT2", "DT3", "changepoints")
    
    bp <- likert(~ DT1 + DT2 + DT3 | changepoints, data=data_to_plot, 
                 scales=list(y=list(at=length(data_to_plot$gene):1, 
                                    relation="free", 
                                    labels=data_to_plot$gene, 
                                    cex=input.cex, font=4), 
                             x=list(at=seq(0, 100, 50),
                                    relation="same",
                                    labels=seq(0, 100, 50), 
                                    cex=input.cex)),
                 strip=FALSE, strip.left=strip.custom(bg=strip.col),
                 par.strip.text=list(cex=input.cex.lab,lines=1.5),
                 between=list(y=c(0,0,0,0)), 
                 xlab=list("Proportion",cex=input.cex.lab), 
                 cex.lab=input.cex.lab, ylab="", main="",
                 Reference=0, auto.key = list(space ="right",columns=1), 
                 col=c(ft.col[1], ft.col[2], ft.col[3]))
    bp <- bp + layer(panel.abline(v=50, col="black", lty=2))
  }

  return(bp)
}
