#' @include utilsTradeSeq.R

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
  for(ll in 1:ncol(cw)){
    curY <- y[cw[,ll] == 1]
    curX <- pt[cw[,ll] == 1,ll]
    smooth[[ll]] <- loess(curY ~ curX)
  }

  ### plot raw and smooth
  ptAll <- pt[which(as.logical(cellAssign))]
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

  for(ll in 1:ncol(cw)){
    ptGrid <- seq(0, max(pt[lineageID == ll,ll]),
                  length.out = 100)
    yhat <- predict(smooth[[ll]],
                    data.frame(curX = ptGrid),
                    type = "response")

    curDf <- data.frame(pt = ptGrid,
                        y = yhat,
                        lineage = factor(ll))
    p1 <- p1 + geom_line(data = curDf,
                   mapping = aes(x=pt, y=yhat, col=lineage),
                   lwd = lwd+1)
  }
  return(p1)
})


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
