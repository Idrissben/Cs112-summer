# Load imports
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
           seed = 12345, cutpoints = NA, median = TRUE){
    set.seed(seed)
    
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

generate_weights <- function(vars, nm = 100) {
  weights <- list(1:nm)
  for (i in 1:nm) {
    nms <- runif(vars)
    weights[[i]] <- diag(nms)
  }
  return(weights)
}

ulti.match <- function(X, Y, Tr, dta, mb.form, md.form, treatment = 'treat', num_weights = 10) {
  
  # setup
  num_vars <- ncol(X)
  starts <- num_weights
  wghts <- generate_weights(num_vars, starts)
  stor <- list(pvalues = 1:starts, gammas = 1:starts, 
               mds = 1:starts, ses = 1: starts, ATTs = 1:starts)
  
  # start
  st.time <- Sys.time()

  for (i in 1:starts) {
    mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=wghts[[i]])
    mb <- MatchBalance(mb.form, match.out=mout, nboots=500, ks=TRUE, print.level = 0)
    ps <- psens(mout, Gamma=2, GammaInc=.1)
    gamm <- ps$bounds[(ps$bounds['Upper bound'] > 0.05),][1,1]
    if (is.na(gamm)) {gamm <- 2}
    
    # model dependence
    matches <- as.matrix(cbind(mout$index.treated, mout$index.control, mout$weights))
    df <- as.data.frame(rbind(dta[matches[,1],], dta[matches[,2],]), stringsAsFactors = FALSE)
    md <- modelDependence_(dataset = df, treatment = treatment, 
                           verbose=FALSE, base.form = md.form, median=TRUE)
    
    stor$pvalues[i] <- mb$AMsmallest.p.value
    stor$gammas[i] <- gamm
    stor$mds[i] <- md
    stor$ses[i] <- mout$se
    stor$ATTs[i] <- mout$est
    
    progress(i, starts, progress.bar = TRUE)
  }
  
  print("")
  cat("Total time taken", Sys.time() - st.time)
  
  return(stor)
}


# Load Data
lalonde <- read.dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta")[,-1]
dw.treat <- lalonde[lalonde$treat == 1,] 
cps <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls.dta")[,-1]
lalonde.cps <- rbind(dw.treat, cps)

# Create data objects
attach(lalonde.cps)
Y <- re78
Tr <- treat
X <- cbind(age, education, black, hispanic, married, nodegree, re74, re75)
BalanceMat <- cbind(age, I(age^2), education, I(education^2), black,
                    hispanic, married, nodegree, re74 , I(re74^2), re75, I(re75^2),
                    I(re74*re75), I(age*nodegree), I(education*re74), I(education*re75))

mb.form <- as.formula('treat ~ age + I(age^2) + education + I(education^2) + black +
  hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
  I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75)')

md.form <- as.formula('re78 ~ treat + age + education + black + 
                      hispanic + married + nodegree + re74 + re75')

# CALL THE FUNCTION
stor100 <- ulti.match(X, Y, Tr, lalonde.cps, mb.form, md.form, num_weights = 100)

plot(stor$pvalues, stor$mds)


# plot where SE is size and md on one axis
sesizes <- rescale(stor.1e6$ses, to = c(40, 5))

fig <- plot_ly(x=stor.1e6$pvalues, 
               y=stor.1e6$mds, 
               z=stor.1e6$gammas, 
               type="scatter3d", 
               mode="markers", 
               color = stor.1e6$ATTs,
               marker = list(size = sesizes),
               hovertext = paste("p value:", stor.1e6$pvalues,
                                 "<br>MD:", stor.1e6$mds,
                                 "<br>Gamma:", stor.1e6$gammas,
                                 "<br>SE:", stor.1e6$ses),
               hoverinfo = "text")

fig <- fig %>% layout(
  title = "Layout options in a 3d scatter plot",
  scene = list(
    xaxis = list(title = "p values", type = "log"),
    yaxis = list(title = "Model dependence"),
    zaxis = list(title = "Gammas")
  ))

fig



# plot where md is size
sizes <- rescale(stor.1e6$mds, to = c(10, 30))

fig <- plot_ly(x=stor.1e6$pvalues, 
               y=stor.1e6$ses, 
               z=stor.1e6$gammas, 
               type="scatter3d", 
               mode="markers", 
               color = stor.1e6$ATTs,
               marker = list(size = sizes),
               hovertext = paste("p value:", stor.1e6$pvalues,
                                 "<br>SE:", stor.1e6$ses,
                                 "<br>Gamma:", stor.1e6$gammas,
                                 "<br>MD:", stor.1e6$mds),
               hoverinfo = "text")

fig <- fig %>% layout(
  title = "Layout options in a 3d scatter plot",
  scene = list(
    xaxis = list(title = "p values", type = "log"),
    yaxis = list(title = "Standard errors"),
    zaxis = list(title = "Gammas")
  ))

fig

idxs <- stor$pvalues > 1e-10

stor.1e6 <- list('pvalues' = stor$pvalues[idxs], 
                 'gammas' = stor$gammas[idxs], 
                 'mds' = stor$mds[idxs], 
                 'ses' = stor$ses[idxs],
                 'ATTs' = stor$ATTs[idxs])

stor.1k <- stor

