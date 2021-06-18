library(slingshot)
data("crv", package = "traviz")
data("sce", package = "traviz")
counts <- counts(sce)
sce@int_metadata$slingshot <- crv

test_that("plotExpression works",{
  p <- plotExpression(counts, crv, rownames(sce)[1])
  expect_true(is(p, "gg"))
  p2 <- plotExpression(counts, slingshot::as.PseudotimeOrdering(crv), rownames(counts)[1])
  expect_true(is(p2, "gg"))
})

test_that("plotSmoothers works", {
  expect_true(is(plotSmoothers(sce, gene = rownames(counts)[1], counts = counts), "gg"))
  # With all edge case options
  expect_true(is(plotSmoothers(sce, gene = 1, counts = counts, border = FALSE), "gg"))
  expect_true(is(plotSmoothers(sce, gene = 1, counts = counts,
                          pointCol = rep("black", ncol(sce))), "gg"))
  sce$color <- rep("black", ncol(sce))
  expect_true(is(plotSmoothers(sce, gene = 1, counts = counts,
                          pointCol = "color"), "gg"))
  expect_message(plotSmoothers(sce, gene = 1, counts = counts,
                               pointCol = rep("black", 3)))
  expect_true(is(plotSmoothers(sce, gene = 1, counts = counts,
                          curvesCols = rep("black", 2)), "gg"))
  expect_message(plotSmoothers(sce, gene = 1, counts = counts,
                               curvesCol = rep("black", 4)))

})
