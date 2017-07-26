#'@name computePS
#'@aliases computePS
#'@title Computes the profile statistics
#'@description Computes the increased or decreased profile statistics for each row (or gene) across the three samples within a feature or data type (expression, methylation, or copy-number variation)
#'@usage computePS(rd1, cd1, cd2, profile)
#'@param rd1 : A numeric value of the reference sample (R) on which to estimate the profile statistics for a given row
#'@param cd1 : A numeric value of the comparison sample 1 (S1) for a given gene, gene probe or gene copy segment
#'@param cd2 : A numeric value of the comparison sample 2 (S2) for a given gene, gene probe or gene copy segment
#'@param profile : The desired direction of genomic change. The values are "up" (default) or "down" to select for increased or decreased gene set profile, respectively
#'@details This function requires three data values corrosponding to three samples for a given gene (or row), respectively
#'@return The returned value is profile statistics computed considering the specified change in the reference sample when compared to the remaining two relative samples. 
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@keywords Profile Statistics
#'@examples
#'rd1 = 40
#'cd1 = 20
#'cd2 = 20
#'computePS(rd1, cd1, cd2, profile="up")
#'@export

computePS = function(rd1, cd1, cd2, profile="up"){

  if(missing(rd1) | missing(cd1) | missing(cd2)){
    stop('One or more sample values is missing!')
  }
  if(is.na(rd1) | is.na(cd1) | is.na(cd2)){
    stop('The sample values are non-numeric!')
  }
  if(length(rd1)>1 | length(cd1)>1 | length(cd2)>1){
    stop('the parameters can not be a list or vector!')
  }

  if(profile=="up"){
    pt1 <- abs( ((rd1 - cd1) / (rd1 - cd2)) - 1 )
    pt2 <- (cd1 * cd2) /  (rd1^2)
    ps <- pt1 + pt2
  }else if(profile=="down"){
    pt1 <- abs ( ((rd1 - cd1) / (rd1 - cd2)) - 1 )
    pt2 <- (rd1^2) /  (cd1 * cd2)
    ps <- pt1 + pt2
  }else{
    stop("the gene set profile is missing!")
  }
  
  return(ps)
}
