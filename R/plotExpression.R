#' @include utilsTradeSeq.R

#' @title Plot gene expression as function of pseudotime.
#' @description Plots a fast loess smoother of gene expression for each lineage.
#' @param counts The matrix of gene expression counts.
#' @param sds  A \code{SlingshotDataSet} or \code{PseudotimeOrdering} object,
#' typically obtained after trajectory inference using \code{Slingshot}.
#' @param gene Gene name of gene to plot.
#' @param type The type of smoother. Defaults to \code{"loess"}.
#' @param span If \code{type} is \code{"loess"}, the \code{span} of the smoother.
#' See \code{loess} documentation.
#' @param lwd Line width of the smoother. Passed to \code{\link{geom_line}}.
#' @param size Character expansion of the data points. Passed to \code{\link{geom_point}}.
#' @param alpha Numeric between 0 and 1, determines the transparency of data points,
#' see \code{scale_color_viridis_d}.
#' @return A \code{ggplot} object.
#' @examples
#' library(ggplot2)
#' data(crv, package="traviz")
#' data(counts, package="traviz")
#' plotExpression(counts = counts, sds=crv, gene=rownames(counts)[1])
#' @import ggplot2
#' @import viridis
#' @rdname plotExpression
#' @export
setMethod(f = "plotExpression",
          signature = c(counts = "matrix",
                        sds = "SlingshotDataSet",
                        gene = "character"),
          definition = function(counts,
                                 sds,
                                 gene,
                                 type = "loess",
                                 span = 0.75,
                                 alpha = 1,
                                 lwd=1,
                                 size = 2/3){

  if(!gene %in% rownames(counts)){
    stop("The gene name is not present in the count matrix rownames.")
  }

  ### assign cells to lineages
  set.seed(33)
  cw <- slingshot::slingCurveWeights(sds)
  cellAssign <- .assignCells(cw)
  pt <- slingPseudotime(sds, na=FALSE)

  ### estimate loess for each lineage
  y <- counts[gene,]
  smooth <- list()
  for(ll in seq_len(ncol(cw))){
    curY <- y[cw[,ll] == 1]
    curX <- pt[cw[,ll] == 1,ll]
    smooth[[ll]] <- loess(curY ~ curX)
  }

  ### plot raw and smooth
  ptAll <- pt[cbind(seq_len(nrow(pt)), apply(cellAssign,1,function(x) which(x==1)))]
  lineageID <- apply(cellAssign, 1, function(x) which(x==1))
  df <- data.frame(y=y,
                   pt=ptAll,
                   lineage=factor(lineageID))
  p1 <- ggplot(df, aes(x=ptAll, y=y, col=lineage)) +
    theme_classic() +
    geom_point(size = size) +
    scale_color_viridis_d(alpha = alpha) +
    xlab("Pseudotime") +
    ylab("Expression")

  for(ll in seq_len(ncol(cw))){
    ptGrid <- seq(0, max(pt[lineageID == ll,ll]),
                  length.out = 100)
    yhat <- predict(smooth[[ll]],
                    data.frame(curX = ptGrid),
                    type = "response")

    assign(paste0("curDf",ll), data.frame(pt = ptGrid,
                        y = yhat,
                        lineage = factor(ll)))

    if(ll == 1){
      allDf <- get(paste0("curDf",ll))
    } else {
      allDf <- rbind(allDf, get(paste0("curDf",ll)))
    }

    p1 <- p1 + geom_line(data = data.frame("ptAll" = ptGrid,
                                           "y" = yhat,
                                           "lineage" = factor(ll)),
                         lwd = lwd+1)
  }


  return(p1)
})


#' @rdname plotExpression
#' @export
setMethod(f = "plotExpression",
          signature = c(counts = "matrix",
                        sds = "PseudotimeOrdering",
                        gene = "character"),
          definition = function(counts,
                                sds,
                                gene,
                                type = "loess",
                                span = 0.75,
                                alpha = 1,
                                lwd=1,
                                size = 2/3){

            p1 <- plotExpression(counts = counts,
                           sds = as.SlingshotDataSet(sds),
                           gene = gene,
                           type = type,
                           span = span,
                           alpha = alpha,
                           lwd = lwd,
                           size =size)
            return(p1)
})
