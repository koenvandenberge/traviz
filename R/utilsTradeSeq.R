# Manipulate the objects to extract meaningul values ----
### lpmatrix given X and design
predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes
  # conditions is the vector of conditions, if present.

  # if pseudotime is vector, make it a matrix.
  if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)

  condPresent <- !is.null(conditions)
  if(condPresent) nConditions <- nlevels(conditions)

  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "t[1-9]):l[1-9]")

  if(!condPresent){
    lineages <- as.numeric(substr(x = colnames(lpmatrix[,allBs]),
                                  start = 4, stop = 4))
    nCurves <- length(unique(lineages))
    for (ii in seq_len(nCurves)) {
      assign(paste0("id",ii), allBs[which(lineages == ii)])
    }
  } else if(condPresent){
    lineages <- as.numeric(substr(x = colnames(lpmatrix[,allBs]),
                                  start = 8, stop = 8))
    nLineages <- length(unique(lineages))
    curves <- as.numeric(substr(x = colnames(lpmatrix[,allBs]),
                                  start = 8, stop = 9))
    nCurves <- length(unique(curves))
    for (ii in seq_len(nLineages)) {
      for(kk in seq_len(nConditions))
      assign(paste0("id",ii, kk), allBs[which(curves == paste0(ii, kk))])
    }
  }


  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  if(!condPresent){
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nCurves)) {
        if (!all(x[get(paste0("id", ii))] == 0)) {
          return(ii)
        }
      }
    })
  } else if(condPresent){
    # first number is lineage, second number is condition.
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nLineages)) {
        # loop over lineages
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (!all(x[get(paste0("id", ii, kk))] == 0)) {
            return(as.numeric(paste0(ii, kk)))
          }
        }
      }
    })
  }


  # fit splinefun for each basis function based on assigned cells
  if(!condPresent) {
    for (ii in seq_len(nCurves)) { # loop over curves
      for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
        assign(paste0("l",ii,".",jj),
               stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                y = lpmatrix[lineageID == ii, #only cells for lineage
                                             get(paste0("id", ii))[jj]],
                                ties = mean)) #basis function
      }
    }
  } else if(condPresent) {
    for (ii in  seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
          #within curve, loop over basis functions
          assign(paste0("l",ii,kk,".",jj),
                 stats::splinefun(
                   x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                   y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                get(paste0("id", ii, kk))[jj]],
                   ties = mean)) #basis function
        }
      }
    }
  }


  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
  if(!condPresent){
    for (ii in seq_len(nCurves)) { # loop over curves
      if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
        for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
          f <- get(paste0("l", ii, ".", jj))
          Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
        }
      }
    }
  } else if(condPresent){
    # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
    for (ii in seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        # loop over conditions
        if (all(df[, paste0("l", ii, kk)] != 0)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
            # within curve, loop over basis functions
            f <- get(paste0("l", ii, kk, ".", jj))
            Xout[, get(paste0("id", ii, kk))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    }
  }


  # add fixed covariates as in df
  dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
  dfOffsetID <- grep(x = colnames(df), pattern = "offset")
  Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]

  # return
  colnames(Xout) <- colnames(lpmatrix)
  return(Xout)
}

# get the first non-errored fit in models
.getModelReference <- function(models){
  for (i in seq_len(length(models))) {
    m <- models[[i]]
    if (is(m)[1] != "try-error") return(m)
  }
  stop("All models errored")
}

# get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId))
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, conditionId))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}


