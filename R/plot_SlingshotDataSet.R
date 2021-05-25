#' @title Plot Slingshot output
#' @name plot-SlingshotDataSet
#' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
#'
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#'
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of
#'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
#'   see Details for more.
#' @param linInd integer, an index indicating which lineages should be plotted
#'   (default is to plot all lineages). If \code{col} is a vector, it will be
#'   subsetted by \code{linInd}.
#' @param show.constraints logical, whether or not the user-specified initial
#'   and terminal clusters should be specially denoted by green and red dots,
#'   respectively.
#' @param add logical, indicates whether the output should be added to an
#'   existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:2}).
#' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
#' @param cex numeric, amount by which points should be magnified, see
#'   \code{\link{par}}.
#' @param lwd numeric, the line width, see \code{\link{par}}.
#' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
#' @param ... additional parameters to be passed to \code{\link{lines}}.
#'
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   principal curves will be plotted.
#'
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' library(slingshot)
#' data("slingshotExample", package="slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' sds <- slingshot(rd, cl, start.clus = "1")
#' plot(sds, type = 'b')
#'
#' # add to existing plot
#' plot(rd, col = 'grey50')
#' lines(sds, lwd = 3)
#'
#' @import graphics
#' @import grDevices
#' @import slingshot
#' @export
setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          linInd = NULL,
                          show.constraints = FALSE,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
        col <- rep(col, length(slingLineages(x)))
        curves <- FALSE
        lineages <- FALSE
        if(is.null(type)){
            if(length(slingCurves(x)) > 0){
                type <- 'curves'
            }else if(length(slingLineages(x)) > 0){
                type <- 'lineages'
            }else{
                stop('No lineages or curves detected.')
            }
        }else{
            type <- c('curves','lineages','both')[pmatch(type,
                                                         c('curves','lineages','both'))]
            if(is.na(type)){
                stop('Unrecognized type argument.')
            }
        }

        if(type %in% c('lineages','both')){
            lineages <- TRUE
        }
        if(type %in% c('curves','both')){
            curves <- TRUE
        }

        if(lineages & (length(slingLineages(x))==0)){
            stop('No lineages detected.')
        }
        if(curves & (length(slingCurves(x))==0)){
            stop('No curves detected.')
        }

        if(is.null(linInd)){
            linInd <- seq_along(slingLineages(x))
        }else{
            linInd <- as.integer(linInd)
            if(!all(linInd %in% seq_along(slingLineages(x)))){
                if(any(linInd %in% seq_along(slingLineages(x)))){
                    linInd.removed <-
                        linInd[! linInd %in% seq_along(slingLineages(x))]
                    linInd <-
                        linInd[linInd %in% seq_along(slingLineages(x))]
                    message('Unrecognized lineage indices (linInd): ',
                            paste(linInd.removed, collapse = ", "))
                }else{
                    stop('None of the provided lineage indices',
                         ' (linInd) were found.')
                }
            }
        }

        if(lineages){
            X <- reducedDim(x)
            clusterLabels <- slingClusterLabels(x)
            connectivity <- slingAdjacency(x)
            clusters <- rownames(connectivity)
            nclus <- nrow(connectivity)
            centers <- t(vapply(clusters,function(clID){
                w <- clusterLabels[,clID]
                return(apply(X, 2, weighted.mean, w = w))
            }, rep(0,ncol(X))))
            rownames(centers) <- clusters
            X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
            clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                           drop = FALSE]
            linC <- slingParams(x)
            clus2include <- unique(unlist(slingLineages(x)[linInd]))
        }

        if(!add){
            xs <- NULL
            ys <- NULL
            if(lineages){
                xs <- c(xs, centers[,dims[1]])
                ys <- c(ys, centers[,dims[2]])
            }
            if(curves){
                npoints <- nrow(slingCurves(x)[[1]]$s)
                xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                              function(c){ c$s[,dims[1]] }, rep(0,npoints))))
                ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                              function(c){ c$s[,dims[2]] }, rep(0,npoints))))
            }
            plot(x = NULL, y = NULL, asp = asp,
                 xlim = range(xs), ylim = range(ys),
                 xlab = colnames(reducedDim(x))[dims[1]],
                 ylab = colnames(reducedDim(x))[dims[2]])
        }

        if(lineages){
            for(i in seq_len(nclus-1)){
                for(j in seq(i+1,nclus)){
                    if(connectivity[i,j]==1 &
                       all(clusters[c(i,j)] %in% clus2include)){
                        lines(centers[c(i,j), dims],
                              lwd = lwd, col = col[1], ...)
                    }
                }
            }
            points(centers[clusters %in% clus2include, dims],
                   cex = cex, pch = 16, col = col[1])
            if(show.constraints){
                if(any(linC$start.given)){
                    points(centers[clusters %in%
                                       linC$start.clus[linC$start.given], dims,
                                   drop=FALSE], cex = cex / 2,
                           col = 'green3', pch = 16)
                }
                if(any(linC$end.given)){
                    points(centers[clusters %in%
                                       linC$end.clus[linC$end.given] &
                                       clusters %in% clus2include,
                                   dims,drop=FALSE], cex = cex / 2,
                           col = 'red2', pch = 16)
                }
            }

        }
        if(curves){
            for(ii in seq_along(slingCurves(x))[linInd]){
                c <- slingCurves(x)[[ii]]
                lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
            }
        }
        invisible(NULL)
    }
)

#' @rdname plot-SlingshotDataSet
#' @import slingshot
#' @export
setMethod(
    f = "lines",
    signature = "SlingshotDataSet",
    definition = function(x,
                          type = NULL,
                          dims = seq_len(2),
                          ...) {
        plot(x, type = type, add = TRUE, dims = dims, ...)
        invisible(NULL)
    }
)

#' @name plot3d-SlingshotDataSet
#' @title Plot Slingshot output in 3D
#'
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#'
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of
#'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
#'   Details for more.
#' @param linInd integer, an index indicating which lineages should be plotted
#'   (default is to plot all lineages). If \code{col} is a vector, it will be
#'   subsetted by \code{linInd}.
#' @param add logical, indicates whether the output should be added to an
#'   existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:3}).
#' @param aspect either a logical indicating whether to adjust the aspect ratio
#'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
#' @param size numeric, size of points for MST (default is \code{10}), see
#'   \code{\link[rgl:plot3d]{plot3d}}.
#' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
#' @param ... additional parameters to be passed to \code{lines3d}.
#'
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   principal curves will be plotted.
#'
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' library(rgl)
#' library(slingshot)
#' data("slingshotExample", package="slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' rd <- cbind(rd, rnorm(nrow(rd)))
#' sds <- slingshot(rd, cl, start.clus = "1")
#' plot3d(sds, type = 'b')
#'
#' # add to existing plot
#' plot3d(rd, col = 'grey50', aspect = 'iso')
#' plot3d(sds, lwd = 3, add = TRUE)
#' }
# #' @importFrom rgl plot3d
#' @import slingshot
#' @export
plot3d.SlingshotDataSet <- function(x,
                                    type = NULL,
                                    linInd = NULL,
                                    add = FALSE,
                                    dims = seq_len(3),
                                    aspect = 'iso',
                                    size = 10,
                                    col = 1,
                                    ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package 'rgl' is required for 3D plotting.",
             call. = FALSE)
    }
    col <- rep(col, length(slingLineages(x)))
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
        if(length(slingCurves(x)) > 0){
            type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
            type <- 'lineages'
        }else{
            stop('No lineages or curves detected.')
        }
    }else{
        type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                            'both'))]
        if(is.na(type)){
            stop('Unrecognized type argument.')
        }
    }

    if(type %in% c('lineages','both')){
        lineages <- TRUE
    }
    if(type %in% c('curves','both')){
        curves <- TRUE
    }

    if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
    }

    if(is.null(linInd)){
        linInd <- seq_along(slingLineages(x))
    }else{
        linInd <- as.integer(linInd)
        if(!all(linInd %in% seq_along(slingLineages(x)))){
            if(any(linInd %in% seq_along(slingLineages(x)))){
                linInd.removed <-
                    linInd[! linInd %in% seq_along(slingLineages(x))]
                linInd <-
                    linInd[linInd %in% seq_along(slingLineages(x))]
                message('Unrecognized lineage indices (linInd): ',
                        paste(linInd.removed, collapse = ", "))
            }else{
                stop('None of the provided lineage indices',
                     ' (linInd) were found.')
            }
        }
    }

    if(lineages){
        X <- reducedDim(x)
        clusterLabels <- slingClusterLabels(x)
        connectivity <- slingAdjacency(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
            w <- clusterLabels[,clID]
            return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                       drop = FALSE]
        clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }

    if(!add){
        xs <- NULL
        ys <- NULL
        zs <- NULL
        if(lineages){
            xs <- c(xs, centers[,dims[1]])
            ys <- c(ys, centers[,dims[2]])
            zs <- c(zs, centers[,dims[3]])
        }
        if(curves){
            npoints <- nrow(slingCurves(x)[[1]]$s)
            xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){
                c$s[,dims[1]] }, rep(0,npoints))))
            ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){
                c$s[,dims[2]] }, rep(0,npoints))))
            zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){
                c$s[,dims[3]] }, rep(0,npoints))))
        }
        rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                    xlim = range(xs), ylim = range(ys), zlim = range(zs),
                    xlab = colnames(reducedDim(x))[dims[1]],
                    ylab = colnames(reducedDim(x))[dims[2]],
                    zlab = colnames(reducedDim(x))[dims[3]])
    }

    if(lineages){
        for(i in seq_len(nclus-1)){
            for(j in seq(i+1,nclus)){
                if(connectivity[i,j]==1 &
                   all(clusters[c(i,j)] %in% clus2include)){
                    rgl::lines3d(x = centers[c(i,j),dims[1]],
                                 y = centers[c(i,j),dims[2]],
                                 z = centers[c(i,j),dims[3]],
                                 col = col[1], ...)
                }
            }
        }
        rgl::points3d(centers[clusters %in% clus2include, dims],
                      size = size, col = col[1])
    }
    if(curves){
        for(ii in seq_along(slingCurves(x))[linInd]){
            c <- slingCurves(x)[[ii]]
            rgl::lines3d(c$s[c$ord,dims], col = col[ii], ...)
        }
    }
    invisible(NULL)
}
