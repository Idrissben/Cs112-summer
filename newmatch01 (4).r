library(Matching)
library(rgenoud)

#
#
#
#
# THE IMPORTANT FUNCTION
NewMatch <- function(Tr, X, pop.size = 10, max.generations = 5, domains = c(1.01, 1.99)) {
  start.time <- Sys.time()
  
  Xs <- X
  
  if (length(Tr) != nrow(Xs)) {
    return('ERROR - INCOMPATIBLE LENGTHS')
  }
  
  n.obs <- length(Tr)
  n.var <- ncol(Xs)
  
  # run an initial GenMatch to get starting parameters
  genout.start <- GenMatch(Tr = Tr, X = Xs, pop.size = pop.size*100, max.generations = max.generations*5, BalanceMatrix = X)
  n.p.values <- length(genout.start$values)
  
  for (i in c(1:n.var)) {
    Xs[, i] <- Xs[, i]*genout.start$par[i]
  }
  
  # prepare domains
  dom <- cbind(rep(domains[1], n.var), rep(domains[2], n.var))
  
  GenMatchWrapper <- function(exponents) {
    
    # print(exponents)
    XN <- Xs
    
    for (i in c(1:n.var)) {
      XN[, i] <- XN[, i]^exponents[i]
    }
    
    # print(head(XN))
    
    genout <- GenMatch(Tr = Tr, X = XN, print.level = 1, project.path = paste(tempdir(), "/genoud.txt", sep = ""), pop.size = pop.size, max.generations = max.generations, BalanceMatrix = X)
    
    return(genout$value[1]) # = highest lowest p-value
  }
  
  genoudout <- genoud(GenMatchWrapper, nvars = n.var, max = TRUE, pop.size = pop.size, max.generations = max.generations, Domains = dom, boundary.enforcement = 2)
  
  # parse the file to find the best result's weights
  file_data <- read.delim(paste(tempdir(), "/genoud.txt", sep = ""), skip = 1, header = FALSE, nrows = 1)
  best.weights <- file_data[1, (2+n.p.values):(1+n.p.values+n.var)]
  best.weights <- as.numeric(best.weights[1,])
  
  XM <- Xs
  for (i in c(1:n.var)) {
    XM[, i] <- XM[, i]^genoudout$par[i]
  }
  
  print(head(XM))
  
  genout.fin <- GenMatch(Tr = Tr, X = XM, pop.size = pop.size*100, max.generations = max.generations*5, BalanceMatrix = X, starting.values = best.weights)
  mout.fin <- Match(Tr = Tr, X = XM, Weight.matrix = genout.fin)
  
  # prepare output
  end.time <- Sys.time()
  outputlist <- list(mout=mout.fin, genout=genout.fin, matches=genout.fin$matches, pvalues=genout.fin$value, weights=genout.fin$par, exponents=genoudout$par, time=end.time - start.time)
  
  print("#############")
  print(genout.fin$matches)
  print(sprintf("p-values: %s", genout.fin$value))
  print(sprintf("weights: %s", genout.fin$par))
  print(sprintf("exponents: %s", genoudout$par))
  print(sprintf("time taken: %s", end.time - start.time))
  return(outputlist)
}
