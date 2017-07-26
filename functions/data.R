#' Gene Variant (EXP) data
#'
#' A dataset containing the genome-wide gene expression values from 3 multiple myeloma cell line samples.
#' The variables are as follows
#'
#'\itemize{
#'  \item gene
#'  \item gene names
#'  \item sample 1 log2 transformed normalized expression count values
#'  \item sample 2 log2 transformed normalized expression count values
#'  \item sample 3 log2 transformed normalized expression count values
#'  }
#'  @source \url{https://research.themmrf.org/
#'  }
#' @name exprset
#' @docType data
#' @keywords datasets
#' @usage exprset
#' @format A data matrix with 1500 genes and 3 samples
NULL

#' Variant (VAR) data
#'
#' A dataset containing the genome-wide gene variant proportion from 3 multiple myeloma cell line samples.
#' The variables are as follows
#'
#'\itemize{
#'  \item gene
#'  \item gene names
#'  \item sample 1 transformed variant proportion data
#'  \item sample 2 transformed variant proportion data
#'  \item sample 3 transformed variant proportion data
#'  }
#'  @source \url{https://research.themmrf.org/
#'  }
#' @name varset
#' @docType data
#' @keywords datasets
#' @usage varset
#' @format A data matrix with 1101 genes variants and 3 samples
NULL

#' Copy Number Variation (CNV) data
#'
#' A dataset containing the genome-wide gene copy change identified in the 3 multiple myeloma cell line samples.
#' The variables are as follows
#'
#'\itemize{
#'  \item gene
#'  \item copy number variation segment ID
#'  \item sample 1 copy change segment mean value
#'  \item sample 2 copy change segment mean value
#'  \item sample 3 copy change segment mean value
#'  }
#'  @source \url{https://research.themmrf.org/
#'  }
#' @name cnvset
#' @docType data
#' @keywords datasets
#' @usage cnvset
#' @format A data matrix with 534 genes copy change and 3 samples
NULL
