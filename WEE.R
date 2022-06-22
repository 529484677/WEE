
#  An example to implement the WEE method

remove(list = ls())
timestart<-Sys.time()


library(MASS)  
library(maxLik)
library(foreach)
library(doParallel)
library("survival")
library('smcure')
source('Beran.R')

set.seed(1)
beta.true = c(.5, .5)
gamma.true =  c(0.3,1,1)
tau = 4  
lambdaC = 0.2 # censoring rate
n = 100
repetition = 400
B = 100


H_0 = function(t){
  f = t
  return(f)
}
Get.Data = function(){
  ## get  Y, delta, Z
  
  Get.random <- function() {
    Z <- matrix(1, nrow = n, ncol = length(beta.true))
    Z1 = runif(n, min = 0, max = 1)
    Z2 = runif(n, min = 0, max = 2)
    
    Z[,1] <- Z1
    Z[,2] <- Z2
    
    # get T0, the survival time of uncure subgroup
    T0 = NULL # survival time of uncured person
    u <- runif(n, min = 0, max = 1)
    for (i in 1:n) {
      # model Cox: S(t|z) = [exp( -H0*e ) -  exp( -tau*e )] / (1 -  exp( -tau*e )), e = exp(bata1*Z1 + beta2*Z2)
      S0 =  function(t){
        f = ( exp(-H_0(t) * exp(Z[i,]%*%beta.true)) - exp(-tau * exp(Z[i,]%*%beta.true)) ) / ( 1- exp(-tau * exp(Z[i,]%*%beta.true)) ) *  ifelse(0<=t&t<=tau,1,0) } 
      T0[i] <- uniroot(function(t) {1 - S0(t) - u[i]}, c(0, tau))$root
    }
    
    
    Z.multi = cbind(1, Z)  # add the intercept
    u = gamma.true[1]*Z.multi[,1] + (Z.multi[,2]*gamma.true[2]) + (Z.multi[,3]*gamma.true[3])
    g = function(u){
      ## P(eta=1|u)
      f = 1 / (1 + exp(-u))  # logit, incidence model 1
    }
    peta = g(u)  # P(eta=1|z)
    eta <- rbinom(n, 1, peta)  # 1 represents subject i is uncured, 0 is cured
    Large = 100
    TT = eta*T0 + (1-eta)*Large  # T, transformed failure time of the entire population
    
    
    ## get C
    # C <- runif(n, 0.01, c0*max(T0))
    C = rexp(n, lambdaC)
    Y = pmin(TT,C)
    delta = rep(0, times = n)
    delta[TT <= C] <- 1
    T1.max = max(Y[delta==1])
    
    ## get true w for omniscient method
    # w = delta + 1-delta * pi(z)*S.u(Y|Z)/(1-pi(z)+pi(z*)S.u(Y|Z))
    pi.z = peta
    S.u = ( exp(-H_0(Y) * exp(Z%*%beta.true)) - exp(-tau * exp(Z%*%beta.true)) ) / ( 1- exp(-tau * exp(Z%*%beta.true)) ) *  ifelse(0<=Y&Y<=tau,1,0)
    w.true = delta + (1-delta) * (pi.z*S.u)/(1-pi.z+pi.z*S.u)*ifelse(Y<=T1.max, 1, 0)
    
    
    # summarize all simulation data
    data = data.frame(Y=Y, delta=delta, Z=Z, eta=eta, w.true=w.true)
    data = data[order(data[,1]),]  # order by Y
    
    
    censor.rate = sum(delta==0)/n
    cure.rate = sum(eta==0)/n  # eta=1 means uncured
    
    return(list(data=data, censor.rate=censor.rate, cure.rate=cure.rate))
  }
  
  data = vector('list', repetition)
  censor.rate = NULL
  cure.rate = NULL
  for (i in 1:repetition) {
    temp = Get.random()
    data[[i]] = temp$data
    censor.rate[i] = temp$censor.rate
    cure.rate[i] = temp$cure.rate
  }
  
  
  return(list(data=data, censor.rate=censor.rate, cure.rate=cure.rate))
}
Transform = function(x) {
  ## change a array to a list
  r.num = nrow(x)  
  y = list()
  for(r in 1:r.num){
    y[[r]] = x[r,]
  }
  return(y)
}
GetSummary = function(){
  
  beta.True <- as.numeric(beta.true)
  n.beta = length(beta.True)
  options(digits=3)  # decimal point
  
  beta.est.our = matrix(0, nrow = repetition, ncol = n.beta)
  cens.rate = cure.rate = NULL
  beta.se.our = matrix(0, nrow = repetition, ncol = n.beta)
  
  
  for(i in 1:repetition){
    beta.est.our[i,] = MyRsult[[i]]$beta.our
    cens.rate[i] = MyRsult[[i]]$cens.rate
    cure.rate[i] = MyRsult[[i]]$cure.rate
    beta.se.our[i,] = MyRsult[[i]]$beta.se.our
  }
  beta.est.avrge.our = apply(beta.est.our, 2, mean)
  
  # get bias, SD, SE, MSE
  bias.our = beta.est.avrge.our - beta.True
  beta.SD.our = apply(beta.est.our, 2, sd)
  beta.SE.our = apply(beta.se.our, 2, mean)
  
  
  # CP
  get.CP = function(TrueVaule, est, se){
    # TrueVaule: dim is p*1
    # est: dim is repetition*P
    # se: dim is repetition*p
    p = length(TrueVaule)
    TrueVaule.mat = matrix(TrueVaule, nrow = repetition, ncol = p, byrow = TRUE)
    CP = apply((est - qnorm(0.975,0,1)*se)<=TrueVaule.mat & TrueVaule.mat<=(est + qnorm(0.975,0,1)*se), 2, mean)
    return(CP)
  }
  beta.CP.our = get.CP(beta.true, beta.est.our, beta.se.our)
  
  
  cens.Rate = round(mean(cens.rate), digits = 3)
  cure.Rate = round(mean(cure.rate), digits = 3)
  
  # summary for beta
  betaname = c('beta1.WEE', 'beta2.WEE')
  our = cbind(bias.our, beta.SD.our, beta.SE.our, beta.CP.our)
  rownames(our) = betaname
  colnames(our) = c('Bias', 'SD', 'SE', 'CP')
  
  
  runningtime = difftime(Sys.time(), timestart, units="mins")
  cat('===========================================================', '\n',
      'Run Time               :', runningtime, 'minus', '\n',
      'Sample size            :', n,   '\n',
      'Right-censoring rate   :', cens.Rate, '\n',
      'Cure rate              :', cure.Rate, '\n'
  )
  print(our)
  cat('===========================================================', '\n','\n','\n')
  
  
  return(all)
}
Get.cureRate = function(){
  cure.rate = NULL
  cens.rate = NULL
  for(i in 1:repetition){
    cens.rate[i] = MyRsult[[i]]$cens.rate
    cure.rate[i] = MyRsult[[i]]$cure.rate
  }
  cens.Rate = round(mean(cens.rate), digits = 3)
  cure.Rate = round(mean(cure.rate), digits = 3)
  return(list(cure.Rate=cure.Rate, cens.Rate=cens.Rate))
}
Repeat = function(r){
  data = Data$data[[r]]  # the data of r-th repetation
  cens.rate = Data$censor.rate[r]
  cure.rate = Data$cure.rate[r]
  data.record = data
  
  Y= data$Y
  delta = data$delta
  w.true = data$w.true; eta = data$eta; 
  data = subset(data, select = -c(eta,w.true))
  Z = as.matrix(data[,-(1:2)])
  T1.max = max(Y[delta==1])
  
  
  ## our 
  p = length(beta.true)
  h = rep(nrow(data)^(-2/7), p)
  wght = NWweight2(Z = Z, Z.sample = Z, hn = h)  # NW weight B_nk(Z)
  S.beran = S.Ti(data = data, w = wght)  # Beran estimator of S(Y.i|Z.i)
  P.T_inf = P.T(data, wght) # cure parbobality, P(T=infinity|z) = 1 - P(eta=1|z) = 1 - pi(z) =  Beran estimator of S(T1.max|Z) 
  case.wght = delta + (1-delta) * (S.beran-P.T_inf)/S.beran *  ifelse(Y<=T1.max, 1, 0)
  
  
  ## algrithm 3: max our weighted partial likelihood function and using NR of maxLik()  
  loglk = function(beta){
    # beta: the parameter to be estimated, dim is p*1
    # loglk = log Prod_i=1^n [fenzi/ fenmu]^delta.i
    # febzi = exp(Z.i*beta)
    # fenmu = sum_{k in R.i} w.k*exp(Z.k*beta) = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta)
    
    # get a upper triangular matrix A
    # the 1th row is (1,1,1,...,1), means I(Y1>=Y1), I(Y2>=Y1), I(Y3>=Y1),..., I(Y_n>=Y_1); 2th row is (0,1,1,...,1), means I(Y1>=Y2), I(Y2>=Y2), I(Y3>=Y2),..., I(Y_n>=Y_2) 
    A = matrix(1, nrow = n, ncol = n)
    A[lower.tri(A)]=0  
    fenzi = exp(Z%*%beta)  # dim is n*1
    fenmu = A %*% (case.wght*exp(Z%*%beta))  # dim is n*1, the i-th is equal to sum_k=1^n I(Y_k >= Y_i) * w_k*exp(Z_k*beta). Note that Y.k have already been ordered 
    L.i = (fenzi / fenmu)^(delta*case.wght)
    f = sum(log(L.i)) 
    return(f)
  }
  myfit = maxLik(loglk, start = rep(0, length(beta.true)), method = 'NR')
  beta.our = myfit$estimate
  
  # standard error of weight (purterbation) bootstrap method for our
  SE.bootstrap = function(b){
    ## compute bootstrap estimates gamma.est and beta.est
    G = G.set[r,b,]
    wght = NWweight2(Z = Z, Z.sample = Z, hn = h, G)  # NW weight B_nk(Z) using input perturbed value G
    S.beran = S.Ti(data = data, w = wght)  # Beran estimator of S(Y.i|Z.i)  using perturbed NW weight B_nk(Z)
    P.T_inf = P.T(data, wght) # cure parbobality, P(T=infinity|z) = 1 - P(eta=1|z) = 1 - pi(z) =  Beran estimator of S(T1.max|Z) using perturbed NW weight B_nk(Z)
    # case.wght = delta + (1-delta) * (S.beran-P.T_inf)/S.beran  # using perturbed S.beran and P.T_inf
    case.wght = delta + (1-delta) * (S.beran-P.T_inf)/S.beran *  ifelse(Y<=T1.max, 1, 0)
    
    
    ## algrithm 3: max log Ln using our partial likelihood function and using maxLik()  using perturbed lam * log Ln
    loglk = function(beta){
      # beta: the parameter to be estimated, dim is p*1
      # loglk = log Prod_i=1^n [fenzi/ fenmu]^delta.i
      # febzi = exp(Z.i*beta)
      # fenmu = sum_{k in R.i} w.k*exp(Z.k*beta) = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta)
      
      # get a upper triangular matrix A
      # the 1th row is (1,1,1,...,1), means I(Y1>=Y1), I(Y2>=Y1), I(Y3>=Y1),..., I(Y_n>=Y_1); 2th row is (0,1,1,...,1), means I(Y1>=Y2), I(Y2>=Y2), I(Y3>=Y2),..., I(Y_n>=Y_2) 
      A = matrix(1, nrow = n, ncol = n)
      A[lower.tri(A)]=0  
      fenzi = exp(Z%*%beta)  # dim is n*1
      fenmu = A %*% (case.wght*exp(Z%*%beta))  # dim is n*1, the i-th is equal to sum_k=1^n I(Y_k >= Y_i) * w_k*exp(Z_k*beta). Note that Y.k have already been ordered 
      # L.i = (fenzi / fenmu)^delta
      L.i = (fenzi / fenmu)^(delta*case.wght)
      f = sum(log(L.i)*G)  # input perturbed value G
      # f = log( prod(L.i) )
      return(f)
    }
    myfit = maxLik(loglk, start = rep(0, length(beta.true)),method = 'NR')
    beta.our = myfit$estimate
    return(beta.our)
  }
  bootstrap.est = sapply(1:B, SE.bootstrap)
  beta.se.our= apply(bootstrap.est, 1, sd)
  
  
  
  return(list(beta.our=beta.our
              , beta.se.our = beta.se.our
              , cens.rate=cens.rate, cure.rate=cure.rate))
  
}


# parallel
no_cores <- detectCores(logical = FALSE)  # the number of core
cl <- makeCluster(no_cores)  # 
registerDoParallel(cl)  # registe parallel
packagelist = c( 'maxLik', 'survival','smcure','MASS')
Data = Get.Data()  # the simulated data over 400 repetations 
G.set = array(0, dim = c(repetition, B, n))  # perturbed value
for (i in 1:repetition) { G.set[i,,] = matrix(data = rexp(B*n, 1), nrow = B, ncol = n) }
MyRsult = foreach(r = 1:repetition, .combine = rbind, .multicombine = TRUE, .packages = packagelist)  %dopar%  Repeat(r)
MyRsult = Transform(MyRsult)
MySummary = GetSummary()
stopImplicitCluster()  # end parallel










