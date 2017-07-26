#'@name GISPA
#'@aliases GISPA
#'@title Gene Integrates Set Profile Analysis
#'@description Identifies gene sets with a similar a prior defined profile using any combination of three feature or data types
#'@usage GISPA(feature,f.sets,g.set,ref.samp.idx,comp.samp.idx,f.profiles,cpt.data,cpt.method,cpt.max)
#'@param feature : Analysis type i.e., one ('1'), two ('2') or three ('3') dimensional feature analysis.
#'@param f.sets : A list of ExpressionSet data objects corrosponding to a data type
#'@param g.set : A GeneSet from an ExpressionSet to subset the f.sets for analysis purposes. 'geneIds' should corrospond to the gene names. Default is null, i.e., genome-wide analysis
#'@param ref.samp.idx : Reference sample column index on which to determine the gene set profile. The default is 3
#'@param comp.samp.idx : The other two relative sample column index against which the profile is being determined. The default is 4 and 5
#'@param f.profiles : A vector of the desired direction of genomic change (or profile) corrosponding to each data type. The values are "up" or "down" for increased and decreased gene set profile, respectively
#'@param cpt.data : Identify changepoints for data using variance (cpt.var), mean (cpt.mean) or both (meanvar). Default is cpt.var.
#'@param cpt.method : Choice of single or multiple changepoint model. Default is "BinSeg".
#'@param cpt.max : The maximum number of changepoints  to search for using "BinSeg" method. Default is 60. This number is dependent on the number of input data points
#'@return The returned value is a data matrix including the original data along with between gene profile statistics and identified changepoints.
#'@import Biobase
#'@import GSEABase
#'@importFrom plyr adply
#'@importFrom data.table data.table
#'@importFrom genefilter rowVars
#'@include computePS.R
#'@include cptModel.R
#'@include cptPlot.R
#'results <- GISPA(feature=1,f.sets=c(exprset),g.set=NULL,
#'                ref.samp.idx=1,comp.samp.idx=c(2,3),
#'                f.profiles=("down"),
#'                cpt.data=var, cpt.method="BinSeg", cpt.max=60)
#'@export

#GISPA Main Function
GISPA = function(feature=1, f.sets, g.set=NULL,
                 ref.samp.idx=1, comp.samp.idx=c(2,3), f.profiles="up", 
                 cpt.data="var", cpt.method="BinSeg", cpt.max=60){
  
  #check the input arguments
  if(feature<1 | feature>3) stop('Number of features should be 1, 2 or 3')
  if(missing(f.sets)) stop("One or more sample values is missing!")
  if(is.na(ref.samp.idx)) stop('Specify reference sample column')
  if(is.na(comp.samp.idx[1])) stop('Specify 1st comparison sample column')
  if(is.na(comp.samp.idx[2])) stop('Specify 2nd comparison sample column')
  #check and assign if a geneSet is provided 
  if(!is.null(g.set)) {
    warning('User geneSet is provided to perform the analysis')
  }
  
  
  for(f in 1:feature){ #number of features
    
    #assign the feature ExpressionSet data sets
    if(is.null(f.sets[f][[1]])) {
      stop(paste('Provide ', f ,' expressionSet', sep=""))
    }
    eset <- f.sets[f]
    
    #subset the data by user-provided geneSet
    if(!is.null(g.set)) {
      eset <- eset[!is.na(geneIds(g.set)),]
    } else{
      
    }
    
    #filter zero samples data columns based on the profile direction
    if(is.na(f.profiles[f])) {
      stop(paste('Provide ', f ,' feature profile', sep=""))
    }
    if(f.profiles[f] == "up"){
      exprs.ftd <- exprs(eset[[1]])[exprs(eset[[1]])[,ref.samp.idx]!=0,]
    } else{
      exprs.ftd <- exprs(eset[[1]])[exprs(eset[[1]])[,comp.samp.idx[1]]!=0 & 
                                    exprs(eset[[1]])[,comp.samp.idx[2]]!=0,]
    }
    
    #compute and filter rows by row variance
    exprs.rvars <- subset(exprs.ftd, rowVars(exprs.ftd)>0)
    
    #transform '0' and '1' values to '0.001' and 0.999
    exprs.rvars[exprs.rvars == 0] <- 0.001
    exprs.rvars[exprs.rvars == 1] <- 0.999
    
    #compute profile statistics(PS)
    ps_val <- adply(exprs.rvars, 1, 
                    function(exprs.rvars) computePS (
                      exprs.rvars[[ref.samp.idx]],
                      exprs.rvars[[comp.samp.idx[1]]],
                      exprs.rvars[[comp.samp.idx[2]]], 
                      profile="up"))
    
    #transform back '0.001' and '0.999' values to 0 and 1
    exprs.rvars[exprs.rvars == 0.001] <- 0
    exprs.rvars[exprs.rvars == 0.999] <- 1
    
    #append the estimate PS to data
    exprs.psval <- cbind(exprs.rvars, ps_val[,2])
    colnames(exprs.psval)[ncol(exprs.psval)] <- paste("ps", f, sep="")
    
    #filter row by non-numeric PS values
    exprs.psftd <- exprs.psval[exprs.psval[,ncol(exprs.psval)]!="Inf" & 
                               exprs.psval[,ncol(exprs.psval)]!="-Inf",]
    exprs.psftd <- exprs.psftd[!is.na(exprs.psftd[,ncol(exprs.psftd)]),]
    
    if (nrow(exprs.psftd)<5){
      stop(paste ('Not enough genes in ', f, ' to proceed'))
    }
    
     #there should be a 'gene' column
      egs <- data.frame(fData(eset[[1]])$'gene')
      colnames(egs) <- c("gene")
      rownames(egs) <- row.names(fData(eset[[1]]))
      
      exprs.tmp <- merge(egs, exprs.psftd,
                           all.y=TRUE, by="row.names", sort=FALSE)
      
      if (nrow(exprs.tmp)<5){
        stop(paste ('row names do not match in ',f,' & feature data',sep=""))
      }

    if(f == 1){
      exprs.eset1 <- exprs.tmp
      rownames(exprs.eset1) <- exprs.eset1$Row.names
      exprs.pm <- cbind(exprs.eset1 [,6])
      colnames(exprs.pm) <- "ps1"
      rownames(exprs.pm) <- exprs.eset1$Row.names
      exprs.psm.gene <- exprs.eset1$gene
      
      if(feature == 1){
        exprs.eset <- exprs.eset1
      }
      
    } else if (f == 2){
      exprs.eset2 <- exprs.tmp
      exprs.eset <- merge(exprs.eset1, exprs.eset2, by="gene", sort=FALSE)
      if (nrow(exprs.eset)<5){
        stop('not enough data points to apply change point model')
      }
      rownames(exprs.eset) <- paste(exprs.eset$Row.names.x,
                                    exprs.eset$Row.names.y,
                                    sep=";")
      exprs.pm <- exprs.eset[,c(6,11)]
      exprs.psm.gene <- exprs.eset$gene
      
    } else if (f == 3){
      exprs.eset3 <- exprs.tmp
      exprs.eset <- merge(exprs.eset, exprs.eset3, by="gene", sort=FALSE)
      if (nrow(exprs.eset)<5){
        stop('not enough data points to apply change point model')
      }
      rownames(exprs.eset) <- paste(exprs.eset$Row.names.x, 
                                    exprs.eset$Row.names.y, 
                                    exprs.eset$Row.names, 
                                    sep=";")
      exprs.pm <- exprs.eset[,c(6,11,16)]
      exprs.psm.gene <- exprs.eset$gene
    } else{
      stop('Number of features allowed are 1, 2 or 3')
    }
    
  }
    
  #estimate within (WGFPS) and between gene feature profile statistics (BGFPS)
  #estimate successive differences of -log10 transformed BGFPS
  #estimate change points using change point model
  #output gene sets by identified change points
  if(is.na(cpt.data)) warning('No cpt.data specified, using default cpt.var')
  if(is.na(cpt.method)) warning('No cpt.method specified, using default BinSeg')
  if(is.na(cpt.max)) warning('No cpt.max specified, using default 60')
  
  if (nrow(exprs.pm) != length(exprs.psm.gene)){
    stop('number of rows do not corrospond to number of genes')
  }
  gispa.output <- cptModel(exprs.pm, exprs.psm.gene, 
                              cpt.data, cpt.method, cpt.max)
  
  ##cptPlot(gispa.output$cpt.dm, gispa.output$locate.cpts)
  
  if(is.null(gispa.output$cpt.dm)){
    stop("No gene sets were identified in the input data sets!")
  }else{
    gispa.output.psm <- merge(exprs.eset, gispa.output$cpt.dm, 
                              by="row.names", all.y=TRUE, sort=FALSE, row.names=TRUE)
    gispa.output.psm <- gispa.output.psm[order(gispa.output.psm[,ncol(gispa.output.psm)-2],decreasing=TRUE),]
    rownames(gispa.output.psm) <- gispa.output.psm[,1]
    return(gispa_results <- list(cpt_out = gispa.output.psm, cpt_plot = gispa.output$cpt.dm, all_cutoffs = gispa.output$locate.cpts))
  }
}