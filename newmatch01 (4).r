library(Matching)
library(rgenoud)

#
#
#
#
#
# THE IMPORTANT FUNCTION
NewMatch <- function(Tr, Xs, mode = 'exp', pop.size = 10, max.generations = 5) {
  start.time <- Sys.time()
  
  if (length(Tr) != nrow(Xs)) {
    return('ERROR - INCOMPATIBLE LENGTHS')
  }
  
  # patch for zeros in log mode
  if (mode == 'log') {
    Xs[Xs == 0] <- 0.0001
  }
  
  n.obs <- length(Tr)
  n.var <- ncol(Xs)
  
  genout_start <- GenMatch(Tr = Tr, X = Xs, pop.size = pop.size*100)
  
  
  for (i in 1:length(genout_start$par)){
    Xs[,i] = Xs[,i] * genout_start$par[i]}  

  GenMatchWrapper <- function(exponents) {
    
    print(exponents)
    XN <- Xs
    
    # MODE SWITCH BETWEEN EXPONENTIATION AND LOGARITHM
    if (mode == 'exp') {
      for (i in c(1:n.var)) {
        XN[, i] <- XN[, i]^exponents[i]
      }
    } else if (mode == 'log') {
      for (i in c(1:n.var)) {
        XN[, i] <- log(XN[, i], exponents[i])
      }
    } else {
      return('ERROR - UNKNOWN MODE')
    }
    
    print(head(XN))
    
    genout <- GenMatch(Tr = Tr, X = XN, print.level = 0, pop.size = pop.size, max.generations = max.generations)

    return(genout$value[1])
  }
  
  # prepare domains
  dom <- cbind(rep(1.01, n.var), rep(1.99, n.var))
  print(dom)
  
  
  genoudout <- genoud(GenMatchWrapper, nvars = n.var, max = TRUE, pop.size = pop.size, max.generations = max.generations, Domains = dom, boundary.enforcement = 2)
  
  XM <- Xs
  # MODE SWITCH
  if (mode == 'exp') {
    for (i in c(1:n.var)) {
      XM[, i] <- XM[, i]^genoudout$par[i]
    }
  } else if (mode == 'log') {
    for (i in c(1:n.var)) {
      XM[, i] <- log(XM[, i], genoudout$par[i])
    }
  } else {
    return('ERROR - UNKNOWN MODE')
  }
  
  print(head(XM))
  
  genout_fin <- GenMatch(Tr = Tr, X = XM, pop.size = 11)
  mout_fin <- Match(Tr = Tr, X = XM, Weight.matrix = genout_fin)
  
  end.time <- Sys.time()
  return(c(mout_fin, XM, end.time - start.time))
    
}




