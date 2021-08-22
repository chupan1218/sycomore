#' A comprehensive function for sycomore
#'
#'
#' The sycomore function is to quantify the synergistic and competitive relationship between miRNAs.
#'
#'
#' @param miRtarget miRNA target binding matrix
#' @param miRNAexpression miRNA expression matrix
#' @param mRNAexpression mRNA expression matrix
#'
#' @return sycomore returns a result in the form of a data frame
#' @export
#' @import parallel
#' @importFrom parallel detectCores mclapply makeCluster clusterExport parLapply stopCluster
#'
#' @examples
#' # load datasets
#' load(system.file("extdata/tcga.brca.testdata.Rdata", package="sycomore"))
#' miRNAexpression <- log2(miRNAexpression + 1)
#' mRNAexpression <- log2(mRNAexpression + 1)
#' # call function
#' sycomore(miRtarget, miRNAexpression, mRNAexpression)
#'
sycomore <- function(miRtarget, miRNAexpression, mRNAexpression){

  ## merge the miRtarget
  miRtarget.merge <- merge(miRtarget, miRtarget, by.x = "mRNA", by.y = "mRNA")
  miRtarget.merge <- miRtarget.merge[-which(miRtarget.merge$miRNA.x == miRtarget.merge$miRNA.y), ]

  ## cat("colnames: ", colnames(miRtarget.merge), "\n")
  ## column names: mRNA, miRNA.x, start.x, end.x, miRNA.y start.y, end.y

  ##--- parallel, multicores on Linux, Windows or MacoS---##
  cores <- detectCores(logical = FALSE)
  ## cat("the number of cores is ", cores, "\n")

  cores <- makeCluster(cores)

  clusterExport(cores, c('implement', 'cor', 'discretize3D', 'PID.measure',
                           'miRtarget.merge', 'miRNAexpression', 'mRNAexpression'), envir = environment())
  system.time({
      inforesults <- parLapply(cores, 1:dim(miRtarget.merge)[1], function(x) {implement(miRtarget.merge[x,1], miRtarget.merge[x,2], miRtarget.merge[x,3],
                                                                                        miRtarget.merge[x,4], miRtarget.merge[x,5], miRtarget.merge[x,6],
                                                                                        miRtarget.merge[x,7], mRNAexpression, miRNAexpression)})
    })
  stopCluster(cores)

  ## results
  inforesults <- do.call("rbind", inforesults)

  colnames(inforesults) <- c("mRNA", "miRNA.x", "start.x", "end.x", "miRNA.y", "start.y", "end.y", "interaction", "Pearson.miRNA.pair",
                             "synergy", "unique.miRNA.x", "unique.miRNA.y", "redundancy")

  return(inforesults)
}


