library(rbounds)
library(rgenoud)
library(Matching)
library(sensitivitymv)
library(foreign)
library(segmented)
library(MatchingFrontier)
library(svMisc)
library(plotly)
library(scales)
library(data.table)
library(stats)

#################################################
######### model dependence utilities ############
########## credit: Vinícius Miranda #############

trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}
getCutpoint_ <- function(dataset, base.form, cov, median = TRUE){
  if(!median){
    mod.form <- as.formula(paste(as.character(base.form[2]),
                                 as.character(base.form[1]),
                                 cov))
    
    if(length(unique(dataset[[as.character(base.form[2])]])) == 2){
      base.mod <- glm(mod.form, data = dataset, family = 'binomial')
    } else{
      base.mod <- lm(mod.form, data = dataset)
    }
    
    ### CHANGE : When the median is zero, this leads segmented() to break
    psi0 <-median(dataset[[cov]])
    if (psi0 == 0) psi0 <- mean(dataset[[cov]])
    ###
    
    seg.reg <- segmented(base.mod, seg.Z=mod.form[c(1,3)], 
                         psi = psi0, control = seg.control(it.max = 10000))
    cutpoint <- seg.reg$psi[2]
    
  }
  else{
    cutpoint <- median(dataset[[cov]])
    
    ### CHANGE : When the median is zero, this leads the lm.fit() in 
    # modelDependence() to break
    ###
    if (cutpoint == 0) cutpoint <- 0.001
  }
  return(cutpoint)
}


# Get model dependence
modelDependence_ <- 
  function(dataset, treatment, base.form, verbose = TRUE, 
           cutpoints = NA, median = TRUE){
    
    base.form <- as.formula(base.form)
    
    covs <- strsplit(as.character(base.form[3]), '\\+')
    covs <- unlist(lapply(covs, trim))
    
    base.theta <- lm(base.form, data = dataset)$coefficients[[treatment]]
    
    if(verbose){
      cat(paste('Estimate from base model:', round(base.theta, 2), '\n'))
    }
    
    N <- nrow(dataset)
    # estimate theta_p
    
    theta.Ps <- c()
    
    for(cov in covs){
      if(cov == treatment){next}
      
      # Formula for this iteration
      this.form <- paste(as.character(base.form[2]),
                         as.character(base.form[1]),
                         paste(covs[!(covs %in% cov)], collapse = ' + '))
      
      
      base.mod <- lm(base.form, data = dataset)
      
      # Split data
      if(length(unique(dataset[[cov]])) == 2){
        split.inds <- dataset[[cov]] == unique(dataset[[cov]])[1]            
        dat1 <- dataset[split.inds,]
        dat2 <- dataset[!split.inds,]
      }else{
        if(cov %in% names(cutpoints)){
          cutpoint <- cutpoints[names(cutpoints) == cov]
        }else{
          cutpoint <- getCutpoint_(dataset, base.form, cov, median)
        }
        split.inds <- dataset[[cov]] < cutpoint
        dat1 <- dataset[split.inds,]
        dat2 <- dataset[!split.inds,]
      }
      
      # Get theta_ps
      dat1.est <- lm(this.form, data = dat1)$coefficients[[treatment]]
      dat2.est <- lm(this.form, data = dat2)$coefficients[[treatment]]
      
      this.theta.p <- dat1.est * (nrow(dat1) / N) + dat2.est * (nrow(dat2) / N)        
      
      if(verbose){
        cat(paste('Estimate from', cov, 'partition:', round(this.theta.p, 2), '\n'))
      }
      theta.Ps <- c(theta.Ps, this.theta.p)      
    }
    
    covs <- covs[!(covs %in% treatment)]
    failed.covs <-covs[is.na(theta.Ps)]
    
    theta.Ps <- theta.Ps[!is.na(theta.Ps)]
    
    sigma.hat.theta <- sqrt(sum((theta.Ps - base.theta) ^ 2) / length(theta.Ps))
    
    return(sigma.hat.theta)
  }


#######################################################
############ various fit functions ####################

# fit function #0 is Genmatch default: pvals
### Requirements
# Balance matrix should NOT contain the outcome or treatment vector

# fit function #1 is Chris's biggest linear combination discrepancy metric

myfit.imb <- function(matches, BM) {
  # ASSUME
  # no outcome and no treatment column in BM
  # BM is scaled to mean=0, sd=1

  # prepare inputs for Chris's metric
  trt.X <- BM[matches[, 1], ]
  ctrl.X <- BM[matches[, 2], ]
  
  # find the difference between treated and control units
  # this returns a matrix with element-wise differences
  t <- trt.X - ctrl.X
  # sum each column, and square the value & then sum each element
  return(sum(colSums(t)^2))
  
}

# fit function #2 is adjusted Vini's robust.fitfunc.plus
####################### change to OUTCOME being FIRST column

myfit.md <- function(matches, BM) {
  ### Requirements
  # treatment column should be named 'treat'
  # Balance matrix should have the FIRST column as the outcome
  # Balance matrix should also contain the treatment vector
  # Balance matrix should have appropriate column names
  
  # Get auxiliary objects
  vars <- colnames(BM) 
  if (is.null(vars)) stop("Balance Matrix should have appropriate column names.") 
  base.form <- paste(
    vars[1],
    "~",
    paste(vars[-1], collapse="+")
  )

  # Compute robustness
  df <- as.data.frame(rbind(BM[matches[,1],], BM[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = 'treat', 
                         verbose=FALSE, base.form = base.form, median=TRUE)
  
  return(md)
}

# fit function #3 is for hidden bias

myfit.gamma <- function(matches, BM) {
  # ASSUME
  # last column of BM is outcome

  outcomes <- BM[, ncol(BM)]

  # prepare inputs for psens
  trt.out <- outcomes[matches[, 1]]
  ctrl.out <- outcomes[matches[, 2]]
  
  # run psens
  ps <- psens(trt.out, ctrl.out, Gamma=2, GammaInc=.05)
  gamm <- ps$bounds[(ps$bounds['Upper bound'] > 0.05),][1,1]
  
  if (is.na(gamm)) {gamm <- 2}
  
  return(-gamm)
}

# fit function #4 is for Standard Error

myfit.se <- function(matches, BM) {
  # ASSUME
  # last column of BM is outcome
  
  outcomes <- BM[, ncol(BM)]
  
  # prepare inputs for calculating se
  trt.out <- outcomes[matches[, 1]]
  ctrl.out <- outcomes[matches[, 2]]
  weights <- matches[, 3]
  
  # copied from Match by Sekhon
  mest  <- sum((trt.out-ctrl.out)*weights)/sum(weights)
  v1  <- trt.out - ctrl.out
  varest  <- sum( ((v1-mest)^2)*weights)/(sum(weights)*sum(weights))
  se.standard  <- sqrt(varest)
  
  return(se.standard)
}

###############################################################
################### ultimatch utilities #######################

generate.weights <- function(n.vars, n.weights) {
  weights <- list(1:n.weights)
  for (i in 1:n.weights) {
    nms <- runif(n.vars)
    weights[[i]] <- diag(nms)
  }
  return(weights)
}

eval.weights <- function(dta, y, tr, M, wghts, mb.form, md.form, BM.md, BM.imb) {
  
  # get data the way we need it
  Tr = dta[[tr]]
  Y = dta[[y]]
  X = dta[, !(names(dta) %in% c(tr, y))]
  
  mout <- Match(Y=Y, Tr=Tr, X=X, M=M, estimand="ATT", Weight.matrix=wghts)
  mb <- MatchBalance(mb.form, data = dta, match.out=mout, nboots=500, ks=TRUE, print.level = 0)
  ps <- psens(mout, Gamma=2, GammaInc=.1)
  gamm <- ps$bounds[(ps$bounds['Upper bound'] > 0.05),][1,1]
  if (is.na(gamm)) {gamm <- 2}
  
  # model dependence
  matches <- as.matrix(cbind(mout$index.treated, mout$index.control, mout$weights))
  df <- as.data.frame(rbind(BM.md[matches[,1],], BM.md[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = tr, 
                         verbose=FALSE, base.form = md.form, median=TRUE)
  
  # Chris's imbalance
  m <- BM.imb[mout$index.treated,] - BM.imb[mout$index.control,]
  imb <- sum(colSums(m)^2)
  
  stor <- c(mb$AMsmallest.p.value, gamm, md, mout$se.standard, imb, mout$est, mout)
  
  return(stor)
}


#############################################################################
#############################################################################
################################ ULTIMATCH ##################################
#############################################################################
#############################################################################


ulti.match <- function(dta, y = 'y', tr = 'treat', M = 1,
                       form.pval = NULL, form.md = NULL, form.imb = NULL,
                       match.sim = TRUE, match.num = 100, genmatch.sim = TRUE, genmatch.num = 3, 
                       pop.size = 50, wait.generations = 5, max.generations = 10,
                       opt.pval = TRUE, opt.gamma = TRUE, opt.md = TRUE, opt.se = TRUE, opt.imb = TRUE,
                       verbose = FALSE) {
  
  ### DATA PREPROCESSING
  # extract important columns and variable names
  Tr <- dta[[tr]]
  Y <- dta[[y]]
  vars <- names(dta)[!(names(dta) %in% c(tr, y))]
  vars.form <- paste(vars, collapse = '+')
  X <- dta[, vars]
  
  # ? is this very memory demanding ?
  ### PREPARING ALL BALANCE MATRICES
  # p-value (MatchBalance) matrix
  if (is.null(form.pval)) {form.pval <- paste('treat', '~', vars.form)}
  BM.pval <- model.frame(form.pval, dta)
  # gamma matrix - actually just the outcome vector
  # model dependence matrix
  if (is.null(form.md)) {form.md <- paste(y, '~', tr, '+', vars.form)}
  BM.md <- model.frame(form.md, dta)
  # standard error matrix - also just the outcome vector
  # linear discrepancy imbalance matrix
  # formula should be of the form '~ colnames'
  if (is.null(form.imb)) {BM.imb <- scale(X)} else {BM.imb <- scale(model.frame(form.imb, dta))}

  ### PREPARING GENMATCH NUMBERS
  if (opt.pval == TRUE) {opt.pval <- genmatch.num}
  if (opt.gamma == TRUE) {opt.gamma <- genmatch.num}
  if (opt.md == TRUE) {opt.md <- genmatch.num}
  if (opt.se == TRUE) {opt.se <- genmatch.num}
  if (opt.imb == TRUE) {opt.imb <- genmatch.num}
  genmatch.fits <- c(rep('pvals', opt.pval), 
                     rep('myfit.gamma', opt.gamma), 
                     rep('myfit.md', opt.md), 
                     rep('myfit.se', opt.se),
                     rep('myfit.imb', opt.imb))
  
  ### PREPARING MAIN OUTPUT
  # total number of outputs to be stored
  tot.weights <- match.num + length(genmatch.fits)
  stor <- data.frame(matrix(NA, nrow = tot.weights, ncol = 6))
  mouts <- 1:tot.weights
  names(stor) <- c('pvalue', 'gamma', 'md', 'se', 'imb', 'ATT')
  
  #####################
  ### RUNNING MATCH ###
  
  # setup - generate large set of random weights
  num.vars <- ncol(X)
  weights <- generate.weights(num.vars, match.num)
  

  ### MAIN MATCHING LOOP
  if (verbose == TRUE) {show('#### PHASE 1: MATCHING ####')}
  pb <- txtProgressBar(min = 0, max = match.num, style = 3) 
  for (i in 1:match.num) {
    
    eval.out <- eval.weights(dta, y, tr, M, weights[[i]], form.pval, form.md, BM.md, BM.imb)
    stor[i, ] <- eval.out[1:6]
    mouts[i] <- eval.out[7]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # display matching results in a neat table
  if (verbose == TRUE) {
    print('#### PHASE 1: SUMMARY ####')
    tbl <- data.frame(min=rep(NA, 5), mean=rep(NA, 5), max=rep(NA, 5))
    tbl[1,] <- list(min(stor$pvalue, na.rm=TRUE), mean(stor$pvalue, na.rm=TRUE), max(stor$pvalue, na.rm=TRUE))
    tbl[2,] <- list(min(stor$gamma, na.rm=TRUE), mean(stor$gamma, na.rm=TRUE), max(stor$gamma, na.rm=TRUE))
    tbl[3,] <- list(min(stor$md, na.rm=TRUE), mean(stor$md, na.rm=TRUE), max(stor$md, na.rm=TRUE))
    tbl[4,] <- list(min(stor$se, na.rm=TRUE), mean(stor$se, na.rm=TRUE), max(stor$se, na.rm=TRUE))
    tbl[5,] <- list(min(stor$imb, na.rm=TRUE), mean(stor$imb, na.rm=TRUE), max(stor$imb, na.rm=TRUE))
    row.names(tbl) <- c('p value', 'gamma', 'model dep', 'st error', 'imbalance')
    print(tbl)
  }
  
  ########################
  ### RUNNING GENMATCH ###
  
  # prep genmatch
  if (verbose == TRUE) {show('#### PHASE 2: GENMATCH ####')}
  
  # select match results that are already best in each category
  gm.pvalue <- weights[order(stor$pvalue, decreasing=TRUE)[1:opt.pval]]
  gm.gamma <- weights[order(stor$gamma, decreasing=TRUE)[1:opt.gamma]]
  gm.md <- weights[order(stor$md, decreasing=FALSE)[1:opt.md]]
  gm.se <- weights[order(stor$se, decreasing=FALSE)[1:opt.se]]
  gm.imb <- weights[order(stor$imb, decreasing=FALSE)[1:opt.imb]]
  genmatch.weights <- c(gm.pvalue, gm.gamma, gm.md, gm.se, gm.imb)

  ### MAIN GENMATCH LOOP
  pb <- txtProgressBar(min = 0, max = length(genmatch.fits), style = 3) 
  for (i in 1:length(genmatch.fits)) {
    
    j <- match.num+i
    
    # use the genmatch.fits list as a source of the correct balance matrix
    fn <- genmatch.fits[i]
    if (fn == 'pvals') {
      BM <- BM.pval
    } else if (fn == 'myfit.md') {
      BM <- BM.md
    } else if (fn == 'myfit.imb') {
      BM <- BM.imb
    } else {
      BM <- Y
    }
    if (fn != 'pvals') {fn <- get(fn)}
    
    # find the weights from the prepped list
    w <- genmatch.weights[[i]]

    genout <- GenMatch(Tr = Tr, X = X, BalanceMatrix = BM, pop.size = pop.size,
                       max.weight = 1, starting.values = diag(w), M = M, max.generations = max.generations,
                       wait.generations = wait.generations, fit.func = fn, print.level = 0)
    
    eval.out <- eval.weights(dta, y, tr, M, diag(genout$par), form.pval, form.md, BM.md, BM.imb)
    
    # store the outcomes
    stor[j, ] <- eval.out[1:6]
    mouts[j] <- eval.out[7]
    weights[[j]] <- diag(genout$par)
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # display best genmatch results
  if (verbose == TRUE) {
    print('#### PHASE 2: SUMMARY ####')
    if (opt.pval != FALSE) {print(paste('best p value increased to ', max(stor$pvalue)))}
    if (opt.gamma != FALSE) {print(paste('best gamma increased to ', max(stor$gamma)))}
    if (opt.md != FALSE) {print(paste('best model dependence decreased to ', min(stor$md)))}
    if (opt.se != FALSE) {print(paste('best standard error decreased to ', min(stor$se)))}
    if (opt.imb != FALSE) {print(paste('best imbalance decreased to', min(stor$imb)))}
  }

  ulti.out <- list(results = stor, weights = weights, 
                   match.num = match.num, genmatch.fits = genmatch.fits)
  return(ulti.out)
}

###############################################################
############### ultimatch viz pipeline fns ####################


ulti.hull <- function(dta, n = 1) {
  # the hulls are being built for
  # the following pairs of variables
  # 1) imb x Gamma
  # 2) pval x Gamma
  # 3) imb x mod dep
  # 4) pval x mod dep
  # 5) imb x pval
  # 6) Gamma x mod dep
  
  npairs <- 6
  pairs <- list(
    cbind(dta$gamma, dta$imb),
    cbind(dta$gamma, dta$pvalue),
    cbind(dta$md, dta$imb),
    cbind(dta$md, dta$pvalue),
    cbind(dta$pvalue, dta$imb),
    cbind(dta$md, dta$gamma)
  )
  
  # n can be a vector of values separate for each hull number
  if (length(n) == 1) {n <- rep(n, npairs)}
  
  hulls <- list()
  hulls.sep <- list()
  
  for (i in 1:npairs) {
    # select the right columns
    pair <- pairs[[i]]
    idx <- c()
    hulls.sep[[i]] <- list()
    # run as many hulls as the user asked for
    for (j in 1:n[i]) {
      
      # calculate the chull
      h <- chull(pair)
      
      # calculate the means of the columns
      # but without the hull units
      pair[h,] <- NA
      m1 <- mean(pair[,1], na.rm = TRUE)
      m2 <- mean(pair[,2], na.rm = TRUE)
      
      # store the indices of this hull
      idx <- c(idx, h)
      hulls.sep[[i]][[j]] <- h
      
      # replace the hull indices by the means
      # this will make them not appear in further hulls
      # but won't remove them from the data (which would shift indices)
      pair[idx,1] <- m1
      pair[idx,2] <- m2
      
    }
    hulls[[i]] <- idx
  }
  
  # put all indices together
  bighull <- sort(unique(unlist(hulls)))
  
  return(list(hulls = hulls, bighull = bighull, hulls.sep = hulls.sep))
}


ulti.prune <- function(dta, perc = 0.1, ignore.se = TRUE,
                       perc.pval = NA, perc.gamma = NA, perc.md = NA,
                       perc.se = NA, perc.imb = NA,
                       ret.pval = NA, ret.gamma = NA, ret.md = NA, 
                       ret.se = NA, ret.imb = NA) {
  
  # retain largest pvalues
  if (is.na(ret.pval)) {
    if (is.na(perc.pval)) {perc.pval <- perc}
    ct.pval <- quantile(dta$pvalue, 1-perc.pval, names = FALSE)
    idx.pval <- which(dta$pvalue >= ct.pval)
  }
  else {
    idx.pval <- which(dta$pvalue >= ret.pval)
    ct.pval <- ret.pval
  }
  
  # retain largest gammas
  if (is.na(ret.gamma)) {
    if (is.na(perc.gamma)) {perc.gamma <- perc}
    ct.gamma <- quantile(dta$gamma, 1-perc.gamma, names = FALSE)
    idx.gamma <- which(dta$gamma >= ct.gamma)
  }
  else {
    idx.gamma <- which(dta$gamma >= ret.gamma)
    ct.gamma <- ret.gamma
  }
  
  # retain lowest model dependence values
  if (is.na(ret.md)) {
    if (is.na(perc.md)) {perc.md <- perc}
    ct.md <- quantile(dta$md, perc.md, names = FALSE)
    idx.md <- which(dta$md <= ct.md)
  }
  else {
    idx.md <- which(dta$md <= ret.md)
    ct.md <- ret.md
  }
  
  # if desired, retain lowest SE values
  ct.se <- NULL
  idx.se <- c()
  if (ignore.se == FALSE) {
    if (is.na(ret.se)) {
      if (is.na(perc.se)) {perc.se <- perc}
      ct.se <- quantile(dta$se, perc.se, names = FALSE)
      idx.se <- which(dta$se <= ct.se)
    }
    else {
      idx.se <- which(dta$se <= ret.se)
      ct.se <- ret.se
    }
  }

  # retain lowest imbalance values
  if (is.na(ret.imb)) {
    if (is.na(perc.imb)) {perc.imb <- perc}
    ct.imb <- quantile(dta$imb, perc.imb, names = FALSE)
    idx.imb <- which(dta$imb <= ct.imb)
  }
  else {
    idx.imb <- which(dta$imb <= ret.imb)
    ct.imb <- ret.imb
  }
  
  cts <- list(pvalue = ct.pval, gamma = ct.gamma, md = ct.md, se = ct.se, imb = ct.imb)
  idxs <- sort(unique(c(idx.pval, idx.gamma, idx.md, idx.se, idx.imb)))
  
  return(list(cutoffs = cts, indices = idxs))

}

ulti.2D <- function(dta, prune = NULL, hull = NULL) {
  
  # plot the following pairs of variables
  # 1) imb x Gamma
  # 2) pval x Gamma
  # 3) imb x mod dep
  # 4) pval x mod dep
  # 5) imb x pval
  # 6) Gamma x mod dep
  
  # if prune is passed as an object,
  # plotting highlights retained values
  
  # if hull is passed as an object,
  # plotting also draws the hulls
  
  npairs <- 6
  pairs <- list(
    cbind(dta$gamma, -dta$imb),
    cbind(dta$gamma, dta$pvalue),
    cbind(-dta$md, -dta$imb),
    cbind(-dta$md, dta$pvalue),
    cbind(dta$pvalue, -dta$imb),
    cbind(-dta$md, dta$gamma)
  )
  
  axs <- list(
    c('Gamma', 'imbalance'),
    c('Gamma', 'p-value'),
    c('model dep.', 'imbalance'),
    c('model dep.', 'p-value'),
    c('p-value', 'imbalance'),
    c('model dep.', 'Gamma')
  )
  
  if (!is.null(prune)) {
    ct <- prune$cutoffs
    # the list are in the order
    # x1 x2 x3 y1 y2 y3
    rcts <- list(
      c(1, ct$gamma, max(dta$gamma), -max(dta$imb), -ct$imb, 0),
      c(1, ct$gamma, max(dta$gamma), 0, ct$pvalue, max(dta$pvalue)),
      c(-max(dta$md), -ct$md, 0, -max(dta$imb), -ct$imb, 0),
      c(-max(dta$md), -ct$md, 0, 0, ct$pvalue, max(dta$pvalue)),
      c(0, ct$pvalue, max(dta$pvalue), -max(dta$imb), -ct$imb, 0),
      c(-max(dta$md), -ct$md, 0, 1, ct$gamma, max(dta$gamma))
    )
  }
  
  
  
  plot.list <- list()
  for (i in 1:npairs) {
    
    pair <- pairs[[i]]
    
    fig <- plot_ly(type="scatter", mode="markers")
    fig <- fig %>%
      add_markers(pair[,1], pair[,2], opacity = 0.5, showlegend=F)

    
    # plot the hulls if they have been passed
    if (!is.null(hull)) {
      hulls <- hull$hulls.sep[[i]]
      for (j in 1:length(hulls)) {
        hl <- hulls[[j]]
        
        fig <- fig %>%
          add_markers(pair[hl,1], pair[hl,2], showlegend=F)
        fig <- fig %>%
          add_polygons(pair[hl,1], pair[hl,2], opacity = 0.1, showlegend=F)
      }
    }
    
    
    # plot the retain zones if they have been passed
    if (!is.null(prune)) {
      rct <- rcts[[i]]
      fig <- fig %>% layout(
        shapes = list(
          list(type = "rect",
               fillcolor = "green",
               opacity = 0.1,
               line = list(width = 0),
               x0 = rct[1], x1 = rct[3],
               y0 = rct[5], y1 = rct[6]),
          list(type = "rect",
               fillcolor = "green",
               opacity = 0.1,
               line = list(width = 0),
               x0 = rct[2], x1 = rct[3],
               y0 = rct[4], y1 = rct[5])
        )
      )
    }
    
    fig <- fig %>% layout(
      xaxis = list(title = list(text = axs[[i]][1], 
                                standoff = 0), 
                   side = "bottom",
                   automargin = TRUE),
      yaxis = list(title = axs[[i]][2]))
    
    plot.list[[i]] <- fig
  }
  
  fig <- subplot(plot.list, 
                 nrows = 3, 
                 titleX = TRUE, 
                 titleY = TRUE,
                 margin = c(0.04, 0.04, 0, 0.08),
                 heights = c(0.36,0.36,0.28))
  
  return(fig)
}


ulti.3D <- function(dta, xax, prune = NULL, ptsize = 20,
                    max.md = NA, max.imb = NA, min.gamma = NA, min.pval = NA) {
  
  ######### USER INPUT NOT IMPLEMENTED YET LOL GOTTA DIRECTLY PLAY WITH CODE
  # default behaviour:
  # X - HAS TO BE SPECIFIED either "p-value" or "imbalance"
  # Y - model dependence
  # Z - hidden bias (gamma)
  # point size - standard error # REMOVED
  # colour - treatment effect
  
  # the prune object, in prune$indices, contains points to plot
  
  if (!is.null(prune)) {dta <- dta[prune$indices,]}
  
  if (is.na(max.md)) {max.md <- max(dta$md)}
  if (is.na(max.imb)) {max.imb <- max(dta$imb)}
  if (is.na(min.gamma)) {min.gamma <- min(dta$gamma)}
  if (is.na(min.pval)) {min.pval <- min(dta$pvalue)}

  fig <- plot_ly(type="scatter3d", mode="markers", showlegend=FALSE)
  
  # 50% probability to flip color scale
  # FOR NOW WE DON'T DO THIS
  # if (sample(c(0,1), size=1) == 0) {
  #   dta$ATT <- rescale(dta$ATT, to=c(max(dta$ATT), min(dta$ATT)))}
  
  if (xax == "imbalance") {xdta <- dta$imb
                           xrng <- c(0, max.imb)}
  if (xax == "p-value") {xdta <- dta$pvalue
                         xrng <- c(min.pval, max(dta$pvalue))}
  
  fig <- fig %>%
    add_markers(x = xdta,
                y = dta$md,
                z = dta$gamma,
                marker = list(size = ptsize,
                              color = dta$ATT,
                              colorscale = "Viridis",
                              showscale = FALSE,
                              line = list(width = 0)),
                hovertext = paste("imbalance", dta$imb,
                                  "<br>p-value:", dta$pvalue,
                                  "<br>MD:", dta$md,
                                  "<br>gamma:", dta$gamma),
                hoverinfo = "text")
  
  fig <- fig %>% layout(
    scene = list(
      xaxis = list(
        title = xax,
        range = xrng),
      yaxis = list(
        title = "Model dep.",
        range = c(0, max.md)),
      zaxis = list(
        title = "Gammas",
        range = c(min.gamma, max(dta$gamma)))
    ))
  
  return(fig)
  
}




#######################################################
############## prep data & formulas ###################


# Load Data
lalonde <- read.dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta")[,-1]
dw.treat <- lalonde[lalonde$treat == 1,] 
cps <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls.dta")[,-1]
lalonde.cps <- rbind(dw.treat, cps)

# add dummy variables re74_0 and re75_0
lalonde.cps$re74_0 <- as.numeric(lalonde.cps$re74 > 0)
lalonde.cps$re75_0 <- as.numeric(lalonde.cps$re75 > 0)

# formulas for Match Balance, Model Dependence and Chris's Imabalance
form.pval <- as.formula('treat ~ age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  re74_0 + re75_0 + I(age*education) + I(age*black) + I(age*hispanic) + I(age*married) + 
  I(age*nodegree) + I(age*re74) + I(age*re75) + I(age*re74_0) + I(age*re75_0) +
  I(education*black) + I(education*hispanic) + I(education*married) + I(education*nodegree) + 
  I(education*re74) + I(education*re75) + I(education*re74_0) + I(education*re75_0) +
  I(black*married) + I(black*nodegree) + I(black*re74) + I(black*re75) + 
  I(black*re74_0) + I(black*re75_0) + I(hispanic*married) + I(hispanic*nodegree) + 
  I(hispanic*re74) + I(hispanic*75) + I(hispanic*re74_0) + I(hispanic*re75_0) +
  I(married*nodegree) + I(married*re74) + I(married*re75) + I(married*re74_0) + I(married*re75_0) +
  I(nodegree*re74) + I(nodegree*re75) + I(nodegree*re74_0) + I(nodegree*re75_0) + I(re74*re75) +
  I(re74*re74_0) + I(re74*re75_0) + I(re75*re75_0) + I(re75*re74_0) + I(re74_0*re75_0)')

form.md <- as.formula('re78 ~ treat + age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  re74_0 + re75_0 + I(age*education) + I(age*black) + I(age*hispanic) + I(age*married) + 
  I(age*nodegree) + I(age*re74) + I(age*re75) + I(age*re74_0) + I(age*re75_0) +
  I(education*black) + I(education*hispanic) + I(education*married) + I(education*nodegree) + 
  I(education*re74) + I(education*re75) + I(education*re74_0) + I(education*re75_0) +
  I(black*married) + I(black*nodegree) + I(black*re74) + I(black*re75) + 
  I(black*re74_0) + I(black*re75_0) + I(hispanic*married) + I(hispanic*nodegree) + 
  I(hispanic*re74) + I(hispanic*75) + I(hispanic*re74_0) + I(hispanic*re75_0) +
  I(married*nodegree) + I(married*re74) + I(married*re75) + I(married*re74_0) + I(married*re75_0) +
  I(nodegree*re74) + I(nodegree*re75) + I(nodegree*re74_0) + I(nodegree*re75_0) + I(re74*re75) +
  I(re74*re74_0) + I(re74*re75_0) + I(re75*re75_0) + I(re75*re74_0) + I(re74_0*re75_0)')

form.imb <- as.formula(' ~ age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  re74_0 + re75_0 + I(age*education) + I(age*black) + I(age*hispanic) + I(age*married) +
  I(age*nodegree) + I(age*re74) + I(age*re75) + I(age*re74_0) + I(age*re75_0) +
  I(education*black) + I(education*hispanic) + I(education*married) + I(education*nodegree) +
  I(education*re74) + I(education*re75) + I(education*re74_0) + I(education*re75_0) +
  I(black*married) + I(black*nodegree) + I(black*re74) + I(black*re75) +
  I(black*re74_0) + I(black*re75_0) + I(hispanic*married) + I(hispanic*nodegree) +
  I(hispanic*re74) + I(hispanic*75) + I(hispanic*re74_0) + I(hispanic*re75_0) +
  I(married*nodegree) + I(married*re74) + I(married*re75) + I(married*re74_0) + I(married*re75_0) +
  I(nodegree*re74) + I(nodegree*re75) + I(nodegree*re74_0) + I(nodegree*re75_0) + I(re74*re75) +
  I(re74*re74_0) + I(re74*re75_0) + I(re75*re75_0) + I(re75*re74_0) + I(re74_0*re75_0)')


###########
## setup ##
# 
# u <- ulti.match(lalonde.cps, y = 're78', tr = 'treat',
#                 form.pval = form.pval, form.md = form.md, form.imb = form.imb, 
#                 match.num = 3, genmatch.num = 1,
#                 pop.size = 10, wait.generations = 2, max.generations = 3, 
#                 verbose = TRUE, M = 1)
# 
# save(u, file='u.RData')
# 
# ulti.plot(u)


# utility function for loading data into a variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# e.g.
# M1C <- loadRData("M1C.RData")

r <- loadRData("u_M1_1500_5_100_10_50_1.RData")
r1 <- loadRData("FINAL-m1.RData")
r2 <- loadRData("FINAL-m2.RData")
r3 <- loadRData("FINAL-m3.RData")

dta <- rbind(r$results, r1$results, r2$results, r3$results)

# using ret.gamma means all gammas over 1.4
# using perc.pval means top 4% of pvalues
prune <- ulti.prune(dta)

# create the 6 2D plots
ulti.2D(dta, prune)

# create the 3D plot with imbalance
ulti.3D(dta, prune, x = "imbalance")
# scuffed
ulti.3D(dta, prune, x = "imbalance", max.imb = 40000, max.md = 1000)


# create the 3D plot with p-value
ulti.3D(dta, prune, x = "p-value")
# scuffed
ulti.3D(dta, prune, x = "p-value", max.md = 1000)

##### CHANGELOG
# DONE: change shape to diamonds
# DONE: use match results for genmatch starting weights
# DONE: setup M=1 / M=2
# DONE: instead of genmatch pvalue <- use Chris's metric

# 50/50 color swap
# enable cropping but ensure color scale remains the same
# enable user choice for the number of diamonds
# enable hiding matching results
# enable user choice of circle transparency

##### checks on user input to be added
# any NaN in BM.imb
# match.num > genmatch.num
# treatment and outcome columns correct


################
################
# this function currecntly unused
ulti.plot <- function(stor, pval.limit = 1e-6) {
  # default behaviour:
  # X - p values
  # Y - model dependence
  # Z - hidden bias
  # point size - standard error
  # colour - treatment effect
  
  # changelog:
  # if genmatch included, plot it as diamonds
  # for now this is always true
  
  # check if pvalues are 0 and remove from log plots
  res <- stor$results
  res$pvalhigh <- res$pvalue > pval.limit
  res$resizedse <- rescale(res$se, to = c(40, 10))
  
  circles <- res[1:stor$match.num, ]
  circles <- circles[which(circles$pvalhigh == TRUE), ]
  diamonds <- res[(1+stor$match.num):nrow(res), ]
  diamonds <- diamonds[which(diamonds$pvalhigh == TRUE), ]
  
  fig <- plot_ly(type="scatter3d", mode="markers", showlegend=FALSE)
  
  fig <- fig %>%
    add_markers(x = circles$pvalue,
                y = circles$md,
                z = circles$gamma,
                marker = list(size = circles$resizedse,
                              color = circles$ATT,
                              colorscale = "Viridis",
                              showscale = TRUE,
                              line = list(color="white")),
                hovertext = paste("p-value:", circles$pvalue,
                                  "<br>MD:", circles$md,
                                  "<br>gamma:", circles$gamma,
                                  "<br>SE:", circles$se),
                hoverinfo = "text")
  
  fig <- fig %>%
    add_markers(x = diamonds$pvalue,
                y = diamonds$md,
                z = diamonds$gamma,
                marker = list(size = diamonds$resizedse, 
                              symbol = "diamond", 
                              color = diamonds$ATT,
                              colorscale = "Viridis",
                              showscale = TRUE,
                              line = list(color="white")),
                hovertext = paste("p-value:", diamonds$pvalue,
                                  "<br>MD:", diamonds$md,
                                  "<br>gamma:", diamonds$gamma,
                                  "<br>SE:", diamonds$se),
                hoverinfo = "text")
  
  fig <- fig %>% layout(
    title = "Layout options in a 3d scatter plot",
    scene = list(
      xaxis = list(title = "p values", type = "log"),
      yaxis = list(title = "Model dependence"),
      zaxis = list(title = "Gammas")
    ))
  
  return(fig)
}

