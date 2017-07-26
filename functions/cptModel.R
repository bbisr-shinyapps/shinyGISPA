#'@name cptModel
#'@aliases cptModel
#'@title Computes within and between gene feature profile statistics by feature and amongst features
#'@description Computes the percentiles on the estimated profile statistics within a gene and across genes for one or more combination of feature or data types (expression, methylation, copy-number variation, or variant change)
#'@usage cptModel(psm, genelist, cpt.data, cpt.method, cpt.max)
#'@param psm : A data matrix of estimated gene profile statistics for each feature
#'@param genelist : A vector of gene names or gene symbols corrosponding to the profile statistics
#'@param cpt.data : Identify changepoints in the data using variance (cpt.var), mean (cpt.mean) or both (meanvar). Default is cpt.var.
#'@param cpt.method : Choice of single or multiple changepoint model. Default is "BinSeg".
#'@param cpt.max : The maximum number of changepoints to search for using "BinSeg" method. Default is 60. This number is dependent on the number of input data points
#'@details This function estimates within and between feature profile statistics by gens in addition to the summed percentiles and successive differences
#'@return Estimated change points in the input data set
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@keywords Profile Statistics
#'@importFrom stats ecdf
#'@importFrom data.table data.table
#'@importFrom data.table setorder
#'@importFrom changepoint cpt.var
#'@importFrom changepoint cpt.mean
#'@importFrom changepoint cpt.meanvar
#'@importFrom changepoint cpts
#'@examples
#'id <- 1000 ## number of probes
#'s <- 3 ## number of sample groups
#'dm <- matrix(runif(id*s,0,200), nrow=id, ncol=s, dimnames=list(paste("gene", 1:id, sep="") , paste("fs", 1:s, sep="")))
#'genelist <- rownames(dm)
#'cptModel(dm, genelist, cpt.data="var", cpt.method="BinSeg", cpt.max=60)
#'@export

cptModel = function(psm, genelist,
                    cpt.data="var", cpt.method="BinSeg", cpt.max=60) {

  .I <- NULL
  fs <- NULL
  rn <- NULL
  
  #apply ecdf to each profile statistics by each gene
  psm_dt <- data.table(genelist, psm, keep.rownames = TRUE)
  psm_dt_copy <- psm_dt
  WGFPS <- 0
    for(f in 3:ncol(psm_dt)){
      colnames(psm_dt)[f] <- 'fs'
      f_peg <- psm_dt[, ecdf(fs)(fs), by=genelist]
      colnames(psm_dt)[f] <- colnames(psm_dt_copy)[f]
      WGFPS <- WGFPS + f_peg$V1
    }
  dt_per <- cbind(psm_dt,WGFPS)
  #carry gene combination with minimum within feature profile statistics
  cof_stat <- dt_per[dt_per[,.I[which.min(WGFPS)], by=genelist] [['V1']]]
  row_header <- cof_stat[,genelist,rn]
  #reprep the data
  cof_stat = cof_stat[, -c("WGFPS", "genelist", "rn")]
  
  #apply ecdf to each profile statistics column across genes
  bgfps <- apply(cof_stat, 2, function(c) ecdf(c)(c))
  BGFPS <- rowSums(bgfps)
  #log transformation
  Neg_log10_BGFPS <- -log(BGFPS,10)
  
  
  bgfps_stat_log <- cbind(row_header,bgfps,BGFPS,Neg_log10_BGFPS)
  #sort log transformed values
  bgfps_stat_log10_sort <- bgfps_stat_log[ order(-Neg_log10_BGFPS), ]
  
  #successive differences of sorted log transformed profile statistics
  diff <- c( rep(0,length(bgfps_stat_log10_sort$Neg_log10_BGFPS)) )
  for (i in 1:length(bgfps_stat_log10_sort$Neg_log10_BGFPS)){
    diff[i] <- bgfps_stat_log10_sort$Neg_log10_BGFPS[i] - 
               bgfps_stat_log10_sort$Neg_log10_BGFPS[i+1]
    if(i == length(bgfps_stat_log10_sort$Neg_log10_BGFPS) ){
      diff[i] <- 0    #assign the last element a zero value
    }
  }
  
  #identify change points for the data using {changepoints}
  max <- length(diff)/2
  if (cpt.method == "AMOC") {
      if(cpt.data == "mean"){
          diff.cpts <- cpt.mean(diff,method=cpt.method,Q=1)
      } else if(cpt.data == "var"){
          diff.cpts <- cpt.var(diff,method=cpt.method,Q=1)
      } else{
          diff.cpts <- cpt.meanvar(diff,method=cpt.method,Q=1)
      }
  }else{
      if(max > cpt.max){
        if(cpt.data == "mean"){
            diff.cpts <- cpt.mean(diff,method=cpt.method,Q=cpt.max)
        } else if (cpt.data == "var") {
            diff.cpts <- cpt.var(diff,method=cpt.method,Q=cpt.max)
        } else {
            diff.cpts <- cpt.meanvar(diff,method=cpt.method,Q=cpt.max)
        }
      } else if(max <= cpt.max){
        if(cpt.data == "mean"){
            diff.cpts <- cpt.mean(diff,method=cpt.method,Q=max)
        } else if (cpt.data == "var") {
            diff.cpts <- cpt.var(diff,method=cpt.method,Q=max)
        } else{
            diff.cpts <- cpt.meanvar(diff,method=cpt.method,Q=max)
        }
      }
  }
  
  #check for the number of diff.cpts identified
  if(length(cpts(diff.cpts))<1){
    stop("No gene sets by changepoints were identified in the data set!")
  }else{
    #add the estimated changepoint locations to the data points  
    level <- character()
    l <- 1
    locate.cpts <- cpts(diff.cpts)
    for(i in 1:length(cpts(diff.cpts)) ) {
      idx <- cpts(diff.cpts)[i]
      locate.cpts[i] <- bgfps_stat_log10_sort$Neg_log10_BGFPS[idx]
      if(i==1){
        for(j in 1:cpts(diff.cpts)[i]){
          level[j]<-l
        }
      }else{
        a <- cpts(diff.cpts)[i-1]+1
        for(j in a:cpts(diff.cpts)[i]){
          level[j]<-l
        }
      }
      if(i==length(cpts(diff.cpts))){
        l<-1000
        a <- cpts(diff.cpts)[i]+1
        for(j in a:length(diff)){
          level[j]<-1000
        }
      }
      l<-l+1
    }
    locate.cpts <- as.numeric(locate.cpts)
    locate.cpts <- round(locate.cpts,digits=2)
  }
  
  #append changepoints to the original data frame
  #generate results table
  cpt.out <- data.frame(bgfps_stat_log10_sort,diff,as.numeric(level))
  setorder(cpt.out, -Neg_log10_BGFPS)

  colnames(cpt.out)[ncol(cpt.out)-2] <- "-log10(BGFPS)"
  colnames(cpt.out)[ncol(cpt.out)-1] <- "succ.diff."
  colnames(cpt.out)[ncol(cpt.out)] <- "changepoints"
  cpt.dm <- as.matrix(cpt.out[,-c(1,2)])
  rownames(cpt.dm) <- cpt.out[,1]
  
  return (data = list(cpt.dm = cpt.dm, locate.cpts = locate.cpts))
}
