#' A plugin function
#'
#' A plugin function for evaluating binding interaction and information values
#'
#' @param mRNA the mRNA
#' @param miRNA.x the miRNA x
#' @param start.x the start site of miRNA x binding on mRNA
#' @param end.x the end site of miRNA x binding on mRNA
#' @param miRNA.y the miRNA y
#' @param start.y the start site of miRNA y binding on mRNA
#' @param end.y the end site of miRNA binding on mRNA
#' @param mRNAexpression mRNA expression matrix
#' @param miRNAexpression miRNA expression matrix
#'
#' @return implement returns a result vector
#' @export
#' @import Informeasure
#' @importFrom Informeasure discretize3D PID.measure
implement <- function(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y,
                      mRNAexpression, miRNAexpression){

  mRNA.index <- which(rownames(mRNAexpression) == mRNA, arr.ind = F)

  miRNA.x.index <- which(rownames(miRNAexpression) == miRNA.x, arr.ind = F)

  miRNA.y.index <- which(rownames(miRNAexpression) == miRNA.y, arr.ind = F)

  results.vector <- c()

  if(length(mRNA.index) & length(miRNA.x.index) & length(miRNA.y.index)){

    if(as.numeric(end.x) < as.numeric(start.y) | as.numeric(end.y) < as.numeric(start.x))
      interaction <- "neighboring"
    else
      interaction <- "overlapping"

    mRNAexpression <- mRNAexpression[mRNA.index, ]
    miRNAexpression <- miRNAexpression[c(miRNA.x.index, miRNA.y.index), ]

    Pearson.miRNA.pair <- cor(as.numeric(miRNAexpression[1, ]), as.numeric(miRNAexpression[2, ]))

    XYZ <- discretize3D(as.numeric(miRNAexpression[1, ]),
                        as.numeric(miRNAexpression[2, ]),
                        as.numeric(mRNAexpression))

    PID <- PID.measure(XYZ)

    results.vector <- c(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y, interaction, Pearson.miRNA.pair,
                        PID$Synergy, PID$Unique_X, PID$Unique_Y, PID$Redundancy)


    remove(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y, mRNAexpression, miRNAexpression,
           interaction, Pearson.miRNA.pair, PID)
  }
  else{

    results.vector <- NULL

    remove(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y, mRNAexpression, miRNAexpression)
  }

  return(results.vector)
}


