library(testthat)

library(sycomore)

load(system.file("extdata/tcga.brca.testdata.Rdata", package="sycomore"))

miRNAexpression  <- log2(miRNAexpression  + 1)
mRNAexpression   <- log2(mRNAexpression   + 1)


results.unittesting <- sycomore(miRtarget, miRNAexpression, mRNAexpression)


test_that("unit test sycomore", {
  expect_equal(sycomore(miRtarget, miRNAexpression, mRNAexpression),  results.unittesting)
})

