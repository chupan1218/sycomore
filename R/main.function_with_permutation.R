#' A comprehensive function for sycomore
#'
#'
#' The sycomore function is to quantify the synergistic and competitive relationship between miRNAs that share a common target mRNA.
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
sycomore <- function(miRtarget,
                     miRNAexpression,
                     mRNAexpression,
                     permutation = TRUE,
                     B = 1000L,
                     permute_target = c("mRNA", "miRNA.x", "miRNA.y"),
                     one_sided = TRUE,
                     p_adjust_method = "BH",
                     seed = 1L){

  permute_target <- match.arg(permute_target)

  ## merge the miRtarget (pair miRNAs on the same mRNA)
  miRtarget.merge <- merge(miRtarget, miRtarget, by.x = "mRNA", by.y = "mRNA")
  miRtarget.merge <- miRtarget.merge[-which(miRtarget.merge$miRNA.x == miRtarget.merge$miRNA.y), ]

  ## parallel: use PSOCK cluster for Linux/Windows/macOS ---##
  n_cores <- parallel::detectCores(logical = FALSE)
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ## make RNG reproducible across workers
  parallel::clusterSetRNGStream(cl, iseed = as.integer(seed))

  ## local helper: compute PID components for a single triplet
  implement_perm <- function(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y,
                             mRNAexpression, miRNAexpression,
                             permutation, B, permute_target, two_sided){

    mRNA.index <- which(rownames(mRNAexpression) == mRNA, arr.ind = FALSE)
    miRNA.x.index <- which(rownames(miRNAexpression) == miRNA.x, arr.ind = FALSE)
    miRNA.y.index <- which(rownames(miRNAexpression) == miRNA.y, arr.ind = FALSE)

    if(!(length(mRNA.index) && length(miRNA.x.index) && length(miRNA.y.index))){
      return(NULL)
    }

    ## binding interaction type
    if(as.numeric(end.x) < as.numeric(start.y) || as.numeric(end.y) < as.numeric(start.x)){
      interaction <- "neighboring"
    }else{
      interaction <- "overlapping"
    }

    x <- as.numeric(mRNAexpression[mRNA.index, ])
    z1 <- as.numeric(miRNAexpression[miRNA.x.index, ])
    z2 <- as.numeric(miRNAexpression[miRNA.y.index, ])

    Pearson.miRNA.pair <- stats::cor(z1, z2, method = "pearson")

    XYZ <- Informeasure::discretize3D(z1, z2, x)
    PID <- Informeasure::PID.measure(XYZ)

    synergy <- PID$Synergy
    unique.x <- PID$Unique_X
    unique.y <- PID$Unique_Y
    redundancy <- PID$Redundancy
    SR <- synergy - redundancy

    ## default: do permutation test 
    p_SR <- NA_real_

    if(isTRUE(permutation)){
      B <- as.integer(B)
      if(B < 1L) stop("'B' must be >= 1 when permutation = TRUE.")
      
      null_SR <- numeric(B)

      for(b in seq_len(B)){
        if(permute_target == "mRNA"){
          x_perm <- sample(x, replace = FALSE)
          z1_perm <- z1; z2_perm <- z2
        }else if(permute_target == "miRNA.x"){
          z1_perm <- sample(z1, replace = FALSE)
          z2_perm <- z2; x_perm <- x
        }else{ ## "miRNA.y"
          z2_perm <- sample(z2, replace = FALSE)
          z1_perm <- z1; x_perm <- x
        }

        XYZp <- Informeasure::discretize3D(z1_perm, z2_perm, x_perm)
        PIDp <- Informeasure::PID.measure(XYZp)
        null_SR[b] <- PIDp$Synergy - PIDp$Redundancy
      }

      ## empirical p-values
      if(isTRUE(one_sided)){
        p_SR <- (1 + sum(null_SR <= abs(SR)))/(B + 1)
      }else{
        p_SR <- (1 + sum(abs(null_SR) >= abs(SR)))/(B + 1)
      }
    }

    if(isTRUE(permutation)){
      c(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y,
        interaction, Pearson.miRNA.pair,
        synergy, unique.x, unique.y, redundancy, SR, p_SR)
    } else {
      c(mRNA, miRNA.x, start.x, end.x, miRNA.y, start.y, end.y,
        interaction, Pearson.miRNA.pair,
        synergy, unique.x, unique.y, redundancy)
    }
  }

  parallel::clusterExport(
    cl,
    varlist = c("miRtarget.merge", "miRNAexpression", "mRNAexpression",
                "implement_perm", "B", "permutation", "permute_target", "one_sided"),
    envir = environment()
  )

  inforesults <- parallel::parLapply(
    cl,
    X = seq_len(nrow(miRtarget.merge)),
    fun = function(ii){
      implement_perm(miRtarget.merge[ii, 1], miRtarget.merge[ii, 2], miRtarget.merge[ii, 3],
                     miRtarget.merge[ii, 4], miRtarget.merge[ii, 5], miRtarget.merge[ii, 6],
                     miRtarget.merge[ii, 7], mRNAexpression, miRNAexpression,
                     permutation, B, permute_target, one_sided)
    }
  )

  inforesults <- do.call("rbind", inforesults)

  if(isTRUE(permutation)){
    colnames(inforesults) <- c("mRNA", "miRNA.x", "start.x", "end.x", "miRNA.y", "start.y", "end.y",
                               "interaction", "Pearson.miRNA.pair",
                               "synergy", "unique.miRNA.x", "unique.miRNA.y", "redundancy",
                               "SRscore", "p.SRscore")
  } else {
    colnames(inforesults) <- c("mRNA", "miRNA.x", "start.x", "end.x", "miRNA.y", "start.y", "end.y",
                               "interaction", "Pearson.miRNA.pair",
                               "synergy", "unique.miRNA.x", "unique.miRNA.y", "redundancy")
  }

  ## add FDR-adjusted q-values if permutation testing is enabled
  if(isTRUE(permutation)){
    inforesults <- as.data.frame(inforesults, stringsAsFactors = FALSE)
    inforesults$p.SRscore    <- as.numeric(inforesults$p.SRscore)
    inforesults$q.SRscore    <- stats::p.adjust(inforesults$p.SRscore, method = p_adjust_method)
  }

  return(inforesults)
}
