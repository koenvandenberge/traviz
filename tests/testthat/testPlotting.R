data(crv, package = "traffic")
data(counts, package = "traffic")
data(sce, package = "traffic")

test_that("plotExpression works",{
  expect_is(plotExpression(counts, crv, rownames(counts)[1]), "gg")
  expect_is(plotExpression(counts, as.PseudotimeOrdering(crv), rownames(counts)[1]), "gg")

})

