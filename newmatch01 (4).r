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
  
  genout_start <- GenMatch(Tr = Tr, X = Xs, pop.size = 11)
  
  
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


# FAKE DATASET TESTING
cities = c("City A","City B","City C","City D","City E")
populations = c(10000,100000,200000,1000000,1200000,1500000)
area = c(11,20,30,50,55,60)
Tre = c(0,1,0,0,1,0)

foo2 = data.frame(cbind(cities,populations,area,Tre))
foo2$Tre <- foo2$Tre == 1
foo2$populations = as.numeric(as.character(populations))
foo2$area = as.numeric(as.character(area))
head(foo2)

# delete columns
foo2X <- foo2[, c(-1, -4)]
head(foo2X)
################### CITIES TEST
OUT2 <- NewMatch(foo2$Tre, foo2X)

#########################
# PREPARING THE PEACEKEEPING DATASET
foo <- read.csv("https://tinyurl.com/y2zt2hyc")
foo <- foo[, c(6:8, 11:16, 99, 50, 114, 49, 63, 136, 109, 126, 48, 160, 142, 10)]
foo <- foo[c(-19, -47), ]
which(is.na(foo) == TRUE)
head(foo)

# match balance before matching; smallest p value
mb_before <- MatchBalance(untype4 ~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + treaty + develop + exp + decade, data=foo, nboots=1000, print.level = 0)
mb_before$BMsmallestVarName
mb_before$BMsmallest.p.value

# THIS MIGHT NEED A REWRITE / DATAFRAME PROBLEMS
treat <- foo$untype4
X <- data.frame(cbind(foo$wartype, foo$logcost, foo$wardur, foo$factnum, foo$factnum2, foo$treaty, foo$develop, foo$exp, foo$decade))
head(X)
class(X$X6)
X[4, 6]^8

genout <- GenMatch(Tr = treat, X = X, print.level = 0)
mout <- Match(Tr = treat, X = X, Weight.matrix = genout)
mb <- MatchBalance(untype4 ~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + treaty + develop + exp + decade, match.out = mout, nboots = 500, print.level = 0, data = foo)
mb$AMsmallestVarName
mb$AMsmallest.p.value

OUTPEACE <- NewMatch(treat, X, pop.size = 10, max.generations = 3)
OUTPEACE2 <- NewMatch(treat, X, pop.size = 10, max.generations = 5) # 1 - 15
OUTPEACE3 <- NewMatch(treat, X, pop.size = 10, max.generations = 5) # -1 - 5
OUTPEACE4 <- NewMatch(treat, X, pop.size = 10, max.generations = 5) # 1 - 10
OUTPEACE5 <- NewMatch(treat, X, pop.size = 10, max.generations = 5) # 1 - 10

OUTPEACELOG <- NewMatch(treat, X, pop.size = 10, max.generations = 5, mode = 'log') # 1 - 10
OUTPEACELOG2 <- NewMatch(treat, X, pop.size = 10, max.generations = 5, mode = 'log') # 1 - 10
OUTPEACELOG3 <- NewMatch(treat, X, pop.size = 10, max.generations = 5, mode = 'log') # 1 - 10


