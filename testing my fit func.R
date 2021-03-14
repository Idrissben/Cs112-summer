library(Matching)

data(lalonde)

#The covariates we want to match on
X = cbind(lalonde$age, lalonde$educ, lalonde$black, lalonde$hisp, lalonde$married, lalonde$nodegr, 
          lalonde$u74, lalonde$u75, lalonde$re75, lalonde$re74)

#The covariates we want to obtain balance on
BalanceMat <- cbind(lalonde$age, lalonde$educ, lalonde$black, lalonde$hisp, lalonde$married, lalonde$nodegr, 
                    lalonde$u74, lalonde$u75, lalonde$re75, lalonde$re74,
                    I(lalonde$re74*lalonde$re75))

library(rbounds)

# [p-value, psens, TE]

generate_weights <- function(vars, nm = 100, mn = 0, mx = 1000) {
  weights <- list(1:nm)
  for (i in 1:nm) {
    nms <- runif(vars, min = mn, max = mx)
      weights[[i]] <- diag(nms)
  }
  return(weights)
}

# setup
starts <- 1000
wghts <- generate_weights(10, starts)
stor <- list(pvalues = 1:starts, psenses = 1:starts, ATTs = 1:starts)

for (i in 1:starts) {
  mout <- Match(Y=Y, Tr=lalonde$treat, X=X, estimand="ATT", Weight.matrix=wghts[[i]])
  mb <- MatchBalance(lalonde$treat~lalonde$age +lalonde$educ+lalonde$black+ lalonde$hisp+ lalonde$married+ lalonde$nodegr+ 
                       lalonde$u74+ lalonde$u75+ lalonde$re75+ lalonde$re74+ I(lalonde$re74*lalonde$re75),
                     match.out=mout, nboots=500, ks=TRUE, print.level = 0)
  ps <- psens(mout, Gamma=1.5, GammaInc=.1)
  stor$pvalues[i] <- mb$AMsmallest.p.value
  stor$psenses[i] <- ps$pval
  stor$ATTs[i] <- mout$est
}

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
rbPall <- rbPal(10)[as.numeric(cut(stor$ATTs,breaks = 10))]

plot(stor$pvalues, stor$psenses, col = rbPall, pch = 20, log = 'xy')

## AT THIS POINT, WE HAVE A 2D PLOT WITH COLOR THAT REPRESENTS P-SENS, WITH PVAL.
