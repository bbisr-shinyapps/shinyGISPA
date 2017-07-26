#'@name cptPlot
#'@aliases cptPlot
#'@title Scatterplot representation of gene sets by change points
#'@description Scatterplot representation of identified change points on the estimated profile statistics within the data
#'@param psv : A data vector of estimated profile statistics on which changepoints are identified
#'@param cut.pts : The estimated profile statistics cutoffs corrosponding to the locations in psv  
#'@param cpt : Change point cutoff to be highlighted. The default is 1
#'@details This function expects 'gispa.output' profile statistics output from GISPA.R main function
#'@return Plot representing all the identified gene sets by change points in the data
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@importFrom data.table setorder
#'@importFrom graphics plot
#'@importFrom graphics par
#'@importFrom graphics text
#'@importFrom graphics abline
#'@examples
#'x <- runif(100, 0.0, 1.0)
#'y <- c(0.2, 0.6, 0.8)
#'cpt.plot <- cptPlot(psv=x,cut.pts=y,cpt=1)
#'@export

cptPlot <- function(psv, cut.pts, cpt){
  if(missing(psv)){-
    stop("data vector is missing!")
  }
  if(missing(cut.pts)){
    stop("changepoint cutoffs to plot are missing!")
  }
  if(missing(cpt)){
    stop("changepoint to highlight are missing!")
  }
  psv <- psv[,ncol(psv)-2]
  psv <- sort(psv,FALSE)
  cut.pts <- sort(cut.pts,TRUE)
  
  #plot the -log10 BGFPS statistics
  par(bg = "white")
  cpt.plot <- plot(psv,
                   type='p', lty=3, lwd=1,
                   main="", xlab="Gene Index", ylab="-Log10(BGFPS)", 
                   cex=2.0, cex.lab =1.5, pch=1, 
                   xlim=c(0,length(psv)))
  
  #plot the change point cut-offs
  for(i in seq_along(cut.pts)) {
    if(i<=cpt) { #user defined number of changepoints
      abline(h=cut.pts[i], 
             lty=2, col='orange')
      text(cut.pts[i], cut.pts[i], cut.pts[i], 
           cex=1.0, pos=4, col="orange")
    }else{
      abline(h=cut.pts[i], 
             lty=2, col='grey')
      text(cut.pts[i], cut.pts[i], cut.pts[i], 
           cex=1.0, pos=4, col="grey")
    }
  }
  
  return(cpt.plot)
}
