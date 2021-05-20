data(crv, package = "traffic")
data(counts, package = "traffic")
data(sce, package = "traffic")

test_that("plotExpression works",{
  expect_type(plotExpression(counts, crv, rownames(counts)[1]), "gg")
  expect_is(plotExpression(counts, as.PseudotimeOrdering(crv), rownames(counts)[1]), "gg")
})

test_that("plotSmoothers works", {
  expect_type(plotSmoothers(sce, gene = rownames(counts)[1], counts = counts), "gg")
  # With all edge case options
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit), border = FALSE), "gg")
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit),
                          pointCol = rep("black", ncol(sdsFit))), "gg")
  sdsFit$color <- rep("black", ncol(sdsFit))
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit),
                          pointCol = "color"), "gg")
  expect_message(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit),
                               pointCol = rep("black", 3)))
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit),
                          curvesCol = rep("black", 2)), "gg")
  expect_message(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit),
                               curvesCol = rep("black", 4)))

})
