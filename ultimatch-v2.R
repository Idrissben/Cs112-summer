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

#################################################
######### model dependence utilities ############

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

# fit function #1 is Genmatch default: pvals
### Requirements
# Malance matrix should NOT contain the outcome or treeatment vector

# fit function #2 is Vini's robust.fitfunc.plus

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

generate.weights <- function(vars, nm) {
  weights <- list(1:nm)
  for (i in 1:nm) {
    nms <- runif(vars)
    weights[[i]] <- diag(nms)
  }
  return(weights)
}

eval.weights <- function(dta, y, tr, wghts, mb.form, md.form) {
  
  # get data the way we need it
  Tr = dta[[tr]]
  Y = dta[[y]]
  X = dta[, !(names(dta) %in% c(tr, y))]
  
  mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=wghts)
  mb <- MatchBalance(mb.form, data = dta, match.out=mout, nboots=500, ks=TRUE, print.level = 0)
  ps <- psens(mout, Gamma=2, GammaInc=.1)
  gamm <- ps$bounds[(ps$bounds['Upper bound'] > 0.05),][1,1]
  if (is.na(gamm)) {gamm <- 2}
  
  # model dependence
  matches <- as.matrix(cbind(mout$index.treated, mout$index.control, mout$weights))
  df <- as.data.frame(rbind(dta[matches[,1],], dta[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = tr, 
                         verbose=FALSE, base.form = md.form, median=TRUE)
  
  stor <- c(mb$AMsmallest.p.value, gamm, md, mout$se.standard, mout$est)
  
  return(stor)
}

ulti.match <- function(dta, mb.form, md.form, y = 'y', tr = 'treat', verbose = FALSE,
                       match.sim = TRUE, match.num = 100, M = 1,
                       genmatch.sim = TRUE, genmatch.num = 3, genmatch.pop = 50, genmatch.gen = 5,
                       opt.pval = TRUE, opt.gamma = TRUE, opt.md = TRUE, opt.se = TRUE) {
  
  # get data the way we need it
  Tr = dta[[tr]]
  Y = dta[[y]]
  X = dta[, !(names(dta) %in% c(tr, y))]
  
  # set up fit function numbers for genmatch
  if (opt.pval == TRUE) {opt.pval <- genmatch.num}
  if (opt.gamma == TRUE) {opt.gamma <- genmatch.num}
  if (opt.md == TRUE) {opt.md <- genmatch.num}
  if (opt.se == TRUE) {opt.se <- genmatch.num}
  genmatch.fits <- c(rep('pvals', opt.pval), 
                     rep('myfit.gamma', opt.gamma), 
                     rep('myfit.md', opt.md), 
                     rep('myfit.se', opt.se))
  
  # total number of outputs to be stored
  tot.weights <- match.num + length(genmatch.fits)

  # first running all match with randomized weights
  # setup
  num.vars <- ncol(X)
  weights <- generate.weights(num.vars, match.num)
  
  stor <- data.frame(matrix(NA, nrow = tot.weights, ncol = 5))
  names(stor) <- c('pvalue', 'gamma', 'md', 'se', 'ATT')
  
  # MAIN MATCHING LOOP
  if (verbose == TRUE) {show('#### PHASE 1: MATCHING ####')}
  pb = txtProgressBar(min = 0, max = match.num, style = 3) 
  for (i in 1:match.num) {
    
    eval.out <- eval.weights(dta, y, tr, weights[[i]], mb.form, md.form)
    stor[i, ] <- eval.out
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # display matching results in a neat table
  if (verbose == TRUE) {
    print('#### PHASE 1: SUMMARY ####')
    tbl <- data.frame(min=rep(NA, 4), mean=rep(NA, 4), max=rep(NA, 4))
    tbl[1,] <- list(min(stor$pvalue, na.rm=TRUE), mean(stor$pvalue, na.rm=TRUE), max(stor$pvalue, na.rm=TRUE))
    tbl[2,] <- list(min(stor$gamma, na.rm=TRUE), mean(stor$gamma, na.rm=TRUE), max(stor$gamma, na.rm=TRUE))
    tbl[3,] <- list(min(stor$md, na.rm=TRUE), mean(stor$md, na.rm=TRUE), max(stor$md, na.rm=TRUE))
    tbl[4,] <- list(min(stor$se, na.rm=TRUE), mean(stor$se, na.rm=TRUE), max(stor$se, na.rm=TRUE))
    row.names(tbl) <- c('p value', 'gamma', 'model dep', 'st error')
    print(tbl)
  }
  
  # prep genmatch
  if (verbose == TRUE) {show('#### PHASE 2: GENMATCH ####')}
  
  # select match results that are already best in each category
  # but always allow default weights to be one of the choices (the last one)
  gm.pvalue <- weights[order(stor$pvalue, decreasing=TRUE)[1:(opt.pval-1)]]
  gm.gamma <- weights[order(stor$gamma, decreasing=TRUE)[1:(opt.gamma-1)]]
  gm.md <- weights[order(stor$md, decreasing=FALSE)[1:(opt.md-1)]]
  gm.se <- weights[order(stor$se, decreasing=FALSE)[1:(opt.se-1)]]
  genmatch.weights <- c(gm.pvalue, 0, gm.gamma, 0, gm.md, 0, gm.se, 0)

  BM <- cbind(X, dta[tr], dta[y])
  
  # MAIN GENMATCH LOOP
  pb = txtProgressBar(min = 0, max = length(genmatch.fits), style = 3) 
  for (i in 1:length(genmatch.fits)) {
    
    j <- match.num+i
    
    # use the genmatch.fits list as a source of optimization function name
    fn <- genmatch.fits[i]
    if (fn == 'pvals') {BaM <- X} else {
      BaM <- BM
      fn <- get(fn)
    }
    
    # find the weights from the prepped list
    if (genmatch.weights[[i]] == 0) {w <- diag(num.vars)} else {w <- genmatch.weights[[i]]}
    
    genout <- GenMatch(Tr = Tr, X = X, BalanceMatrix = BaM, pop.size = genmatch.pop,
                       max.weight = 1, diag(w), M = M,
                       wait.generations = genmatch.gen, fit.func = fn, print.level = 0)
    
    eval.out <- eval.weights(dta, y, tr, diag(genout$par), mb.form, md.form)
    
    # store the outcomes
    stor[j, ] <- eval.out
    weights[[j]] <- diag(genout$par)
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # display best genmatch results
  if (verbose == TRUE) {
    print('#### PHASE 2: SUMMARY ####')
  }
  if (verbose == TRUE && opt.pval != FALSE) {
    print(paste('best p value increased to ', max(stor$pvalue)))
  }
  if (verbose == TRUE && opt.gamma != FALSE) {
    print(paste('best gamma increased to ', max(stor$gamma)))
  }
  if (verbose == TRUE && opt.md != FALSE) {
    print(paste('best model dependence decreased to ', min(stor$md)))
  }
  if (verbose == TRUE && opt.se != FALSE) {
    print(paste('best standard error decreased to ', min(stor$se)))
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
  res$resizedse <- rescale(res$se, to = c(50, 10))
  
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
                              showscale = FALSE,
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
                              showscale = FALSE,
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

# formulas for MD and MB
mb.form <- as.formula('treat ~ age + education + black + 
                      hispanic + married + nodegree + re74 + re75')

md.form <- as.formula('re78 ~ treat + age + education + black + 
                      hispanic + married + nodegree + re74 + re75')


###########
## setup ##

u <- ulti.match(lalonde.cps, y = 're78', tr = 'treat', mb.form = mb.form, md.form = md.form, 
           match.num = 1000, genmatch.num = 25, genmatch.pop = 50, genmatch.gen = 20, verbose = TRUE,
           M = 1)

save(u, file='u.RData')

ulti.plot(u)

# changelog notes
# DONE: change shape to diamonds
# DONE: use match results for genmatch starting weights
# DONE: setup M=1 / M=2
# - instead of genmatch pvalue <- use Chris's metric





