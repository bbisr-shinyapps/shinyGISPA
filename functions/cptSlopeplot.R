#'@name cptSlopeplot
#'@aliases cptSlopeplot
#'@title Scatterplot representation of identified change points gene set slopes
#'@description This function will plot the average slopes estimated over all gene sets within each change point by data types
#'@usage cptSlopeplot(gispa.output,feature,type,cpt)
#'@param gispa.output : A data matrix of between gene feature profile statistics for each feature with corrosponding identified changepoints. The row names should corrospond to genes or names to be displayed on y-axis
#'@param feature : Analysis type i.e., one ('1'), two ('2') or three ('3') dimensional feature analysis.
#'@param type : Type of data, e.g., EXP  (default) for expression, VAR of variants, CNV for copy number change. 
#'@param cpt : Change point cutoff to be highlighted. The default is 1
#'@details This function expects the output from GISPA function of GISPA package, and highlights the gene set slope profile in the selected changepoints
#'@return Scatterplot illustrating the average slopes by change point to access the best gene set profile
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@import scatterplot3d
#'@importFrom data.table data.table
#'@importFrom graphics plot
#'@importFrom graphics par
#'@importFrom graphics text
#'@importFrom stats lm
#'@examples
#'id <- 200 ## number of rows
#'s <- 3 ## number of columns
#'dm <- matrix(runif(id*s,0,10), nrow=id, ncol=s, 
#'                  dimnames=list(paste("gene", 1:id, sep=""), 
#'                  paste("sample", 1:s, sep="")))
#'changepoints <- sort(sample(1:2, id, replace=TRUE))
#'dm <- cbind(dm,changepoints)
#'cptSlopeplot(gispa.output=dm,feature=1,type="EXP",cpt=1)
#'@export

cptSlopeplot <- function(gispa.output, feature=1, type="EXP", cpt=1){
  
  changepoints <- NULL # Setting the variables to NULL first
  ##select for the changepoint of interest
  gispa.output <- gispa.output[gispa.output[,ncol(gispa.output)]!=1000,]
  
  if(feature==1){
    #data type 1
    subset_data <- gispa.output[,c(3:5,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    colnames(slope) <- "slope"
    intercept <- t(t(lm.r$coeff[2,]))
    colnames(intercept) <- "intercept"
    subset_data <- cbind(subset_data,intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope <- dt_data[,mean(slope),by=changepoints]
    
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) <= cpt] <- "orange"
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) > cpt] <- "grey"
    x <- as.numeric(avg_slope$changepoints)
    y <- as.numeric(avg_slope$V1)
    
    par(bg = "white")
    slopePlot <- plot(x, y,
                      xlim=c(1,max(x)+1),
                      ylim=c(min(y),max(y)),
                      xlab="",
                      ylab=paste("Mean Slope", " (", type[1], ")", sep=""),
                      cex.lab =1.5,
                      pch=16,
                      cex=5,
                      col=avg_slope$cptcolor)
    text(avg_slope$V1, labels=avg_slope$changepoints,cex=1.5,pos=4,offset=0.2)
    
  }
  
  if(feature==2){
    #data type 1
    subset_data <- gispa.output[,c(3:5,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    colnames(slope) <- "slope"
    intercept <- t(t(lm.r$coeff[2,]))
    colnames(intercept) <- "intercept"
    subset_data <- cbind(subset_data,intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope_type_1 <- dt_data[,mean(slope),by=changepoints]
    
    #data type 2
    subset_data <- gispa.output[,c(7:9,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    colnames(slope) <- "slope"
    intercept <- t(t(lm.r$coeff[2,]))
    colnames(intercept) <- "intercept"
    subset_data <- cbind(subset_data,intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope_type_2 <- dt_data[,mean(slope),by=changepoints]
    
    #Merge the data ####
    avg_slope <- merge(avg_slope_type_1,
                       avg_slope_type_2,
                       by=c("changepoints"))
    #plot the data
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) <= cpt] <- "orange"
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) > cpt] <- "grey"
    x <- as.numeric(avg_slope$V1.x)
    y <- as.numeric(avg_slope$V1.y)
    par(bg = "white")
    slopePlot <- plot(x, y,
                      xlim=c(min(x),max(x)+1),
                      ylim=c(min(y),max(y)+1),
                      xlab=paste("Mean Slope", " (", type[1], ")", sep=""),
                      ylab=paste("Mean Slope", " (", type[2], ")", sep=""),
                      cex.lab =1.5,
                      pch=16,
                      cex=5,
                      col=avg_slope$cptcolor)
    text(avg_slope$V1.x, avg_slope$V1.y, 
         labels=avg_slope$changepoints,
         cex=1.5,pos=4,offset=1.2)
  }
  if(feature==3){
    #data type 1
    subset_data <- gispa.output[,c(3:5,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    colnames(slope) <- "slope"
    intercept <- t(t(lm.r$coeff[2,]))
    colnames(intercept) <- "intercept"
    subset_data <- cbind(subset_data,intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope_type_1 <- dt_data[,mean(slope),by=changepoints]
    
    #data type 2
    subset_data <- gispa.output[,c(7:9,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    intercept <- t(t(lm.r$coeff[2,]))
    subset_data[,c("intercept","slope")]<-rbind(intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope_type_2 <- dt_data[,mean(slope),by=changepoints]
    
    #data type 3
    subset_data <- gispa.output[,c(11:13,ncol(gispa.output))]
    #considering we only have three groups
    lm.r <- lm (t(subset_data[, 1:3]) ~ I(1:3) )
    #to get the slope estimates
    slope <- t(t(lm.r$coeff[2,])) 
    colnames(slope) <- "slope"
    intercept <- t(t(lm.r$coeff[2,]))
    colnames(intercept) <- "intercept"
    subset_data <- cbind(subset_data,intercept,slope)
    #take average slope for each changepoint
    dt_data <- data.table(subset_data)
    avg_slope_type_3 <- dt_data[,mean(slope),by=changepoints]
    
    #Merge the data ####
    avg_slope_1_2 <- merge(avg_slope_type_1,
                           avg_slope_type_2,
                           by=c("changepoints"))
    avg_slope <- merge(avg_slope_1_2,
                       avg_slope_type_3,
                       by=c("changepoints"))
    #plot the data
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) <= cpt] <- "orange"
    avg_slope$cptcolor[as.numeric(avg_slope$changepoints) > cpt] <- "grey"
    x <- as.numeric(avg_slope$V1.x)
    y <- as.numeric(avg_slope$V1.y)
    z <- as.numeric(avg_slope$V1)
    par(bg = "white")
    slopePlot <- scatterplot3d(x, y, z,
                               xlim=c(min(x),max(x)+1),
                               ylim=c(min(y),max(y)+1),
                               zlim=c(min(z),max(z)+1),
                               xlab=paste("Mean Slope"," (",type[1],")",sep=""),
                               ylab=paste("Mean Slope"," (",type[2],")",sep=""),
                               zlab=paste("Mean Slope"," (",type[3],")",sep=""),
                               cex.lab =1.5,
                               pch=19,
                               cex.symbols = 5,
                               type="h",
                               color=avg_slope$cptcolor)
    # convert 3D coords to 2D projection
    slopePlot.coords <- slopePlot$xyz.convert(x, y, z)
    text(slopePlot.coords$x, slopePlot.coords$y,
         labels=avg_slope$changepoints,
         cex=1.5,
         pos=4,
         offset=1.2) 
  }
  return (slopePlot)
}