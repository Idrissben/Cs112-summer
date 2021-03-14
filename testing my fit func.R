library(Matching)

data(lalonde)

#The covariates we want to match on
X = cbind(lalonde$age, lalonde$educ, lalonde$black, lalonde$hisp, lalonde$married, lalonde$nodegr, 
          lalonde$u74, lalonde$u75, lalonde$re75, lalonde$re74)

#The covariates we want to obtain balance on
BalanceMat <- cbind(lalonde$age, lalonde$educ, lalonde$black, lalonde$hisp, lalonde$married, lalonde$nodegr, 
                    lalonde$u74, lalonde$u75, lalonde$re75, lalonde$re74,
                    I(lalonde$re74*lalonde$re75))

#
#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#so the population size has been set to be only 16 via the 'pop.size'
#option. This is *WAY* too small for actual problems.
#For details see http://sekhon.berkeley.edu/papers/MatchingJSS.pdf.
#

my.fitfunc <- function(matches, BM) 
{
  # my.fitfunc is a function to calculate balance by the mean
  # standardized difference in the eQQ plot for each variable.
  # Minimize the maximum of these differences across variables.
  # Lexical optimization is conducted.  This is the same as
  # fit.func="qqmean.max"
  
  # note: "matches" has three columns:
  # column 1: index of treated obs
  # column 2: index of control obs
  # column 3: weights for matched-pairs
  
  # note: BM is the BalanceMatrix the user passed into GenMatch
  
  index.treated <- matches[,1]
  index.control <- matches[,2]
  
  nvars <- ncol(BM)
  qqsummary   <- c(rep(NA,nvars))
  
  for (i in 1:nvars)
  {    
    
    qqfoo <- qqstats(BM[,i][index.treated], BM[,i][index.control], standardize=TRUE)
    qqsummary[i] <- qqfoo$meandiff
  } #end of for loop
  
  return(sort(qqsummary, decreasing=TRUE))
  # [5, 4, 4, 2, 1, 1, 1]
}

genout <- GenMatch(Tr=lalonde$treat, X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=1,
                   pop.size=16, max.generations=10, wait.generations=1,
                   fit.func=my.fitfunc)

#The outcome variable
Y=lalonde$re78/1000

#
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=lalonde$treat, X=X, estimand="ATT", Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#                        
mb <- MatchBalance(lalonde$treat~lalonde$age +lalonde$educ+lalonde$black+ lalonde$hisp+ lalonde$married+ lalonde$nodegr+ 
                     lalonde$u74+ lalonde$u75+ lalonde$re75+ lalonde$re74+ I(lalonde$re74*lalonde$re75),
                   match.out=mout, nboots=500, ks=TRUE)

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


