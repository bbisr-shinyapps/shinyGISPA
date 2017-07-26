#'@name getBoxPlot
#'@aliases getBoxPlot
#'@title A plotting function for each sample proportion over the changepoints
#'@description Shows the overall distribution of gene set in the selected changepoint for each data type  for each sample group.
#'@usage getBoxPlot(DF,cpt)
#'@param DF : A data frame containing genes as rows followed by sample data points to be plotted with the estimated changepoints.
#'@param cpt : changepoint gene set to be plotted
#'@details This function expects the output from GISPA function of GISPA.NGS package, and highlights the gene set in the selected changepoints and their overall proportion in each of the three sample groups.
#'@return Boxplot pdf illustrating every sample proportion over gene set selected
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@import HH
#'@export

getBoxPlot <- function(DF,cpt=1,ref,cmp1,cmp2,Y_LABEL,samp_col){

  #select for the changepoint of interest
  DF_select <- subset(DF,DF$changepoints<=cpt)
  #generate box-plots of genes for each data tyoe for the selected set of changepoints
  sample.names = c(colnames(DF_select)[ref], colnames(DF_select)[cmp1], colnames(DF_select)[cmp2])
  bb <- boxplot(DF_select[,c(ref,cmp1,cmp2)], col = c(samp_col[1],samp_col[2],samp_col[3]), names=c(sample.names[1],sample.names[2],sample.names[3]), ylab=Y_LABEL, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
  return(bb)
}
