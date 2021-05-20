#' @export
#' @title Plot the log-transformed counts and the fitted values for a particular
#'  gene along all lineages
#' @name plotSmoothers
#' @param ... parameters including:
setGeneric(
  name = "plotSmoothers",
  #signature = 'counts',
  def = function(models, ...) {
    standardGeneric("plotSmoothers")
  }
)

#' @export
#' @title Plot gene expression in reduced dimension.
#' @description Plot the gene in reduced dimensional space.
#' @name plotGeneCount
#' @param ... parameters including:
setGeneric(
  name = "plotGeneCount",
  #signature = 'curve',
  def = function(curve, ...) {
    standardGeneric("plotGeneCount")
  }
)


#' @export
#' @title Plot gene expression along pseudotime.
#' @description Plot gene expression along pseudotime.
#' @name plotExpression
#' @param ... parameters including:
setGeneric(
  name = "plotExpression",
  def = function(counts,
                 sds,
                 gene,
                 ...) {
    standardGeneric("plotExpression")
  }
)
