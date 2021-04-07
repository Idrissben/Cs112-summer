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
########## credit: Vinicius Miranda #############

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

myfit.md <- function(matches, BM) {
  ### Requirements
  # treatment column should be named 'treat'
  # Balance matrix should have the last column as the outcome
  # Balance matrix should also contain the treatment vector
  # Balance matrix should have appropriate column names
  
  # Get auxiliary objects
  vars <- colnames(BM) 
  nvars <- length(vars)
  if (is.null(vars)) stop("Balance Matrix should have appropriate column names.") 
  base.form <- paste(
    vars[nvars],
    "~",
    paste(vars[-nvars], collapse="+")
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
################### ultimatch #################################

generate.weights <- function(n.vars, n.weights) {
  weights <- list(1:n.weights)
  for (i in 1:n.weights) {
    nms <- runif(n.vars)
    weights[[i]] <- diag(nms)
  }
  return(weights)
}

eval.weights <- function(dta, y, tr, M, wghts, mb.form, md.form, BM.imb) {
  
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
  df <- as.data.frame(rbind(dta[matches[,1],], dta[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = tr, 
                         verbose=FALSE, base.form = md.form, median=TRUE)
  
  # Chris's imbalance
  BM.imb <- BM.imb[matches[,1],] - BM.imb[matches[,2],]
  imb <- sum(colSums(BM.imb)^2)
  
  stor <- c(mb$AMsmallest.p.value, gamm, md, mout$se.standard, imb, mout$est, mout)
  
  return(stor)
}

ulti.match <- function(dta, y = 'y', tr = 'treat', M = 1,
                       form = NULL, form.pval = NULL, form.md = NULL, form.imb = NULL,
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
  if (is.null(form.md)) {form.md <- paste('treat', '~', vars.form, '+', y)}
  BM.md <- model.frame(form.md, dta)
  # standard error matrix - also just the outcome vector
  # linear discrepancy imbalance matrix
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
    
    eval.out <- eval.weights(dta, y, tr, M, weights[[i]], form.pval, form.md, BM.imb)
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
    row.names(tbl) <- c('p value', 'gamma', 'model dep', 'st error', 'lin. disc.')
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
    
    eval.out <- eval.weights(dta, y, tr, M, diag(genout$par), form.pval, form.md, BM.imb)
    
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
  res$resizedse <- rescale(res$se, to = c(40, 5))
  
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


#######################################################
############## prep data & formulas ###################


# Load Data
lalonde <- read.dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta")[,-1]
dw.treat <- lalonde[lalonde$treat == 1,] 
cps <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls.dta")[,-1]
lalonde.cps <- rbind(dw.treat, cps)

# formulas for Match Balance, Model Dependence and Chris's Imabalance
form.pval <- as.formula('treat ~ age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75)')

form.md <- as.formula('re78 ~ treat + age + education + black + 
                      hispanic + married + nodegree + re74 + re75')

form.imb <- as.formula('treat ~ age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75)')


###########
## setup ##

ut <- ulti.match(lalonde.cps, y = 're78', tr = 'treat',
                form.pval = form.pval, form.md = form.md, form.imb = form.imb, 
                match.num = 1500, genmatch.num = 5, 
                pop.size = 500, wait.generations = 5, max.generations = 30, 
                verbose = TRUE, M = 1)

# save(u, file='u.RData')

ulti.plot(ut)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# M1 <- loadRData("M1.RData")

# changelog notes
# DONE: change shape to diamonds
# DONE: use match results for genmatch starting weights
# DONE: setup M=1 / M=2
# DONE: instead of genmatch pvalue <- use Chris's metric

# 50/50 color swap
# enable cropping but ensure color scale remains the same
# enable user choice for the number of diamonds
# enable hiding matching results
# enable user choice of circle transparency

