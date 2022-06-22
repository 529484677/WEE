
# An example to implete the SWEE method

remove(list = ls())
timestart<-Sys.time()


library(foreach)
library(doParallel)
source('Beran.R')


seed = 1 
beta.true = c(.5, .5)
gamma.true = c(-1,1,1)  # cure rate c(.3,1,1) , high rate c(-1,1,1)) 
tau = 4 
lambdaC = 0.2 # censoring rate
repetition = 400
B = 100  # bootstrap times


n = 10000  # total sample size
k = 20  # sample size of a batch
m = n/k  # number of batches


r1 = 0.6  # learn rate
alp = 0.7 # learn rate


H_0 = function(t){
  f = t
  return(f)
}
Get.Data = function(n){
  ## get  Y, delta, Z
  
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
      f = ( exp(-H_0(t) * exp(Z[i,]%*%beta.true)) - exp(-H_0(tau) * exp(Z[i,]%*%beta.true)) ) / ( 1- exp(-H_0(tau) * exp(Z[i,]%*%beta.true)) ) *  ifelse(0<=t&t<=tau,1,0) } 
    T0[i] <- uniroot(function(t) {1 - S0(t) - u[i]}, c(0, tau))$root
  }
  
  ## model P(eta=1|Z) = ..., u = gamma0 + gamma1*Z1 + gamma2*Z2
  Z.multi = cbind(1, Z)  # add the intercept
  u = gamma.true[1]*Z.multi[,1] + (Z.multi[,2]*gamma.true[2]) + (Z.multi[,3]*gamma.true[3])
  g = function(u){
    ## P(eta=1|u)
    f = 1 / (1 + exp(-u))  # logit, incidence model 1
  }
  peta = g(u)  # P(eta=1|z)
  eta <- rbinom(n, 1, peta)  # 1 represents that subject i is uncured, 0 is cured
  Large = 100
  TT = eta*T0 + (1-eta)*Large  # T, transformed failure time of the entire population
  
  
  ## get C
  C = rexp(n, lambdaC)
  Y = pmin(TT,C)
  delta = rep(0, times = n)
  delta[TT <= C] <- 1
  T1.max = max(Y[delta==1])
  
  ## get true w for the omniscient method
  # w = delta + 1-delta * pi(z)*S.u(Y|Z)/(1-pi(z)+pi(z*)S.u(Y|Z))
  pi.z = peta
  S.u = ( exp(-H_0(Y) * exp(Z%*%beta.true)) - exp(-H_0(tau) * exp(Z%*%beta.true)) ) / ( 1- exp(-H_0(tau) * exp(Z%*%beta.true)) ) *  ifelse(0<=Y&Y<=tau,1,0)
  w.true = delta + (1-delta) * (pi.z*S.u)/(1-pi.z+pi.z*S.u)*ifelse(Y<=T1.max, 1, 0)
  
  
  # summarize all simulation data
  data = data.frame(Y=Y, delta=delta, Z=Z, eta=eta, w.true=w.true)
  data = data[order(data[,1]),]  # order by Y
  
  
  censor.rate = sum(delta==0)/n
  cure.rate = sum(eta==0)/n  # eta=1 means uncured
  
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
  
  beta.est.our = beta.se.our = matrix(0, nrow = repetition, ncol = n.beta)
  cens.rate = cure.rate = NULL
  
  
  for(i in 1:repetition){
    beta.est.our[i,] = MyRsult[[i]]$beta.our
    beta.se.our[i,] = MyRsult[[i]]$beta.se.our
    cens.rate[i] = MyRsult[[i]]$censor.rate
    cure.rate[i] = MyRsult[[i]]$cure.rate
  }
  
  cens.Rate = round(mean(cens.rate), digits = 3)
  cure.Rate = round(mean(cure.rate), digits = 3)
  beta.est.avrge.our = apply(beta.est.our, 2, mean)
  
  # get bias, SD, SE
  bias.our = beta.est.avrge.our - beta.True
  beta.SD.our = apply(beta.est.our, 2, sd)
  beta.SE.our = apply(beta.se.our, 2, mean)
  
  
  # CP
  get.CP = function(TrueVaule, est, se){
    # TrueVaule: dim is p*1
    # est: dim is repetition*P
    # se: dim is repetition*p
    p = length(TrueVaule)
    if(!is.matrix(se)){ se = matrix(se, nrow = repetition, ncol = p, byrow = TRUE)}  # make se to be a matrix
    TrueVaule.mat = matrix(TrueVaule, nrow = repetition, ncol = p, byrow = TRUE)
    CP = apply((est - qnorm(0.975,0,1)*se)<=TrueVaule.mat & TrueVaule.mat<=(est + qnorm(0.975,0,1)*se), 2, mean)
    return(CP)
  }
  beta.CP.our = get.CP(beta.true, beta.est.our, beta.se.our)
  
  
  # summary for beta
  betaname = c('beta1.SWEE', 'beta2.SWEE')
  our = cbind(bias.our, beta.SD.our, beta.SE.our, beta.CP.our); 
  rownames(our) = betaname
  colnames(our) = c('Bias', 'SD', 'SE', 'CP')
  
  runningtime = difftime(Sys.time(), timestart, units="mins")
  cat('===========================================================', '\n',
      'Run Time                :', runningtime, 'minus', '\n',
      'Nample size  n          :', n,   '\n',
      'Sub dataset size k      :', k,   '\n',
      'Number of subdatasets m :', m,   '\n',
      'Right-censoring rate    :', cens.Rate, '\n',
      'Cure rate               :', cure.Rate, '\n'
  )
  print(our)
  cat('===========================================================', '\n','\n','\n')
  
}
Un = function(beta, data, case.wght){
  # sore function Un(beta) = 1/n * sum_i=1^n w.i*delta.i*(Z - fenzi/fenmu)
  # febzi = sum_{k in R.i} w.k*exp(Z.k*beta) * Z = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta) * Z, dim is n*p
  # fenmu = sum_{k in R.i} w.k*exp(Z.k*beta) = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta), dim is n*1
  # beta: dim is p*1, and p >= 2
  # data: dim is n*(2+p), data = (Y,delta,Z)
  
  
  Y= data$Y
  delta = data$delta
  Z = as.matrix(data[,-(1:2)])
  
  n = nrow(data)
  A = matrix(1, nrow = n, ncol = n)
  A[lower.tri(A)]=0  
  
  
  # partial derivative using Z1
  fenzi = A %*% (case.wght*exp(Z%*%beta)*Z[,1])  # dim is n*1, the i-th row is equal to sum_k=1^n I(Y_k >= Y_i) * { w_k*exp(Z_k*beta) * Z_k1}, the the 1-th component of Z. Note that Y.k have already been ordered 
  fenmu = A %*% (case.wght*exp(Z%*%beta))  # dim is n*1, the i-th is equal to sum_k=1^n I(Y_k >= Y_i) * w_k*exp(Z_k*beta). Note that Y.k have already been ordered 
  f.1 =  (delta*case.wght) * (Z[,1] - replace(fenzi/as.vector(fenmu), is.na(fenzi/as.vector(fenmu)),0))  # dim is n*1
  f.1 = sum(na.omit(f.1))  # the cured subject with w=0 may be NaN caused by fenzi/fenmu = 0/0, but in fact Un = w*delta * (Z - fenzi/fenmu). So w=0 shuld have been removed in Un 
  
  # partial derivative using Z2
  fenzi = A %*% (case.wght*exp(Z%*%beta)*Z[,2])  # dim is n*1, the i-th row is equal to sum_k=1^n I(Y_k >= Y_i) * { w_k*exp(Z_k*beta) * Z_k2}, the the 2-th component of Z. Note that Y.k have already been ordered 
  fenmu =  fenmu
  f.2 = (delta*case.wght) * (Z[,2] - replace(fenzi/as.vector(fenmu), is.na(fenzi/as.vector(fenmu)), 0))  # dim is n*1
  f.2 = sum(na.omit(f.2))  # the cured subject with w=0 may be NaN caused by fenzi/fenmu = 0/0, but in fact Un = w*delta * (Z - fenzi/fenmu). So w=0 shuld have been removed in Un 
  
  f = c(f.1, f.2) 
  return(-f)
}
Un.array = function(beta, data, case.wght.boot, G){
  # sore function Un(beta) = 1/n * sum_i=1^n w.i*delta.i*(Z - fenzi/fenmu), dim is p*B
  # febzi = sum_{k in R.i} w.k*exp(Z.k*beta) * Z = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta) * Z, dim is n*p
  # fenmu = sum_{k in R.i} w.k*exp(Z.k*beta) = sum_{k=1}^n I(Y.k>=Y.i) * w.k*exp(Z.k*beta), dim is n*1
  # beta: dim is p*1, and p >= 2
  # data: dim is n*(2+p), data = (Y,delta,Z)
  #
  # beta: dim is p*B
  # case.wght: dim is n*B, B is the bootstrap times
  # G: dim is n*B, the perturbation value, n is the size of a sub dataset
  
  B = ncol(G)
  Y= data$Y
  delta = data$delta
  Z = as.matrix(data[,-(1:2)])
  
  n = nrow(data)
  A = matrix(1, nrow = n, ncol = n)
  A[lower.tri(A)]=0  
  
  
  # partial derivative using Z1
  delta.mat = matrix(delta, nrow = n, ncol = B, byrow = FALSE)
  Z1.mat = matrix(Z[,1], nrow = n, ncol = B, byrow = FALSE)
  exp.Zbeta = exp(apply(beta, 2, function(x) Z%*%x))  # dim is n*B, exp(Z%*%beta)
  fenzi = apply(case.wght.boot*exp(Z%*%beta)*Z1.mat, 2, function(x) A%*%x)  # dim is n*B
  fenmu = apply(case.wght.boot*exp(Z%*%beta), 2, function(x) A%*%x)  # dim is n*B
  # the cured subject with w=0 may be NaN caused by fenzi/fenmu = 0/0, but in fact Un = w*delta * (Z - fenzi/fenmu). So w=0 shuld have been removed in Un 
  f.1 =  G * (delta.mat*case.wght.boot) * (Z1.mat - replace(fenzi/fenmu,is.na(fenzi/fenmu), 0))  # if fenzi/fenmu=NaN, then make fenzi/fenmu=0
  f.1 = apply(f.1, 2, na.omit)
  if(is.list(f.1)){
    f.1 = sapply(f.1, sum)
  }else{
    f.1 = apply(f.1, 2, sum)
  }
  
  # partial derivative using Z2
  Z2.mat = matrix(Z[,2], nrow = n, ncol = B, byrow = FALSE)
  fenzi = apply(case.wght.boot*exp(Z%*%beta)*Z2.mat, 2, function(x) A%*%x)  # dim is n*B
  fenmu = fenmu
  # the cured subject with w=0 may be NaN caused by fenzi/fenmu = 0/0, but in fact Un = w*delta * (Z - fenzi/fenmu). So w=0 shuld have been removed in Un 
  f.2 =  G * (delta.mat*case.wght.boot) * (Z2.mat - replace(fenzi/fenmu,is.na(fenzi/fenmu), 0))  # if fenzi/fenmu=NaN, then make fenzi/fenmu=0
  f.2 = apply(f.2, 2, na.omit)
  if(is.list(f.2)){
    f.2 = sapply(f.2, sum)
  }else{
    f.2 = apply(f.2, 2, sum)
  }
  
  f = rbind(f.1, f.2)  # dim is p*B
  return(-f)
}
Repeat = function(rr){
  
  set.seed(seed+rr)
  Num.sub = m  # the Number of sub datasets
  n.sub = k  # the size of a sub dataset
  p = length(beta.true)
  
  beta.all =  matrix(0, nrow = Num.sub, ncol = p)  # beta.est of all sub datasets
  beta.est = rep(0, p)  # iterative initial value beta.est.i and i=0
  censor.rate = cure.rate = NULL
  
  
  beta.est.boot = matrix(0, nrow = p, ncol = B)  # dim is p*B, iterative initial value beta.est.i and i=0
  beta.boot.all = array(0, dim = c(Num.sub, p, B))  # dim is i*p*B, beta.est of all sub datasets under B times of bootstrap 
  
  
  for (i in 1:Num.sub) {
    
    # learn rate r.i
    r = r1*i^(-alp)
    
    ## get the i-th dataset , i=1,...,n
    data.fit = Get.Data(n.sub)  #  get the i-th dataset of smaple size being n.sub
    data = data.fit$data
    censor.rate[i] = data.fit$censor.rate
    cure.rate[i] = data.fit$cure.rate
    
    Y= data$Y
    delta = data$delta
    w.true = data$w.true; eta = data$eta; 
    data = subset(data, select = -c(eta,w.true))
    Z = as.matrix(data[,-(1:2)])
    T1.max = max(Y[delta==1])
    
    
    ## our method
    h = rep(nrow(data)^(-2/7), p)
    wght = NWweight2(Z = Z, Z.sample = Z, hn = h, lam = rep(1, n.sub))  # NW weight B_nk(Z)
    S.beran = S.Ti(data = data, w = wght)  # Beran estimator of S(Y.i|Z.i)
    P.T_inf = P.T(data, wght) # cure parbobality, P(T=infinity|z) = 1 - P(eta=1|z) = 1 - pi(z) =  Beran estimator of S(T1.max|Z) 
    case.wght = delta + (1-delta) * (S.beran-P.T_inf)/S.beran *  ifelse(Y<=T1.max, 1, 0)
    beta.est = beta.est - r*Un(beta.est, data, case.wght)  # beta.i = beta.{i-1} - r.i * s.i(beta.{i-1})
    beta.all[i,] = beta.est  # bete.est of all sub dataset i, i=1,...,n
    # get beta.est under perturbed random number to compute beta.se, dim is i*p*B
    if(B<=1){
      # do not conduct boostrap if B<=1
      beta.est.boot = matrix(0, nrow = p, ncol = B)  # dim is p*B, iterative initial value beta.est.i and i=0
      beta.boot.all[i,,] = beta.est.boot  # dim is i*p*B, bete.est.boot of all sub dataset i, i=1,...,n
    }else{
      # G = matrix(1, nrow = n.sub, ncol = B, byrow = TRUE)
      G = matrix(rexp(n.sub*B, 1), nrow = n.sub, ncol = B, byrow = TRUE)  # perturbation values
      wght.boot = NWweight2.array(Z = Z, Z.sample = Z, hn = h, G)  # dim is n*n*B, NW weight B_nk(Z) using input perturbed value G
      S.beran.boot = S.Ti.array(data = data, w = wght.boot)  # dim is n*n*B, Beran estimator of S(Y.i|Z.i)  using perturbed NW weight B_nk(Z)
      P.T_inf.boot = P.T.array(data, wght.boot) # dim is n*B, cure parbobality, P(T=infinity|z) = 1 - P(eta=1|z) = 1 - pi(z) =  Beran estimator of S(T1.max|Z) using perturbed NW weight B_nk(Z)
      case.wght.boot = apply((S.beran.boot-P.T_inf.boot)/S.beran.boot, 2, function(x) delta+(1-delta)*x*ifelse(Y<=T1.max, 1, 0) )  # dim is n*B, perturbed weight w
      beta.est.boot = beta.est.boot - r*Un.array(beta.est.boot, data, case.wght.boot, G)  # dim is p*B,  using perturbation values G
      beta.boot.all[i,,] = beta.est.boot  # dim is i*p*B, bete.est.boot of all sub dataset i, i=1,...,n
    }
    
  }
  censor.rate = mean(censor.rate)
  cure.rate = mean(cure.rate)
  
  
  beta.our = apply(beta.all, 2, mean)  # beta.bar = 1/n sum_i=1^n beta.est.i
  beta.our.boot = apply(beta.boot.all, 3, function(x)  rep(1,Num.sub)%*%x/Num.sub)  # dim is p*B, beta.bar under perturbed bootstrap
  beta.se.our = apply(beta.our.boot, 1, sd)
  
  
  remove(beta.est, beta.est.boot)  # make initial value of beta to be zero again
  remove(beta.boot.all, beta.all)  # save memory
  
  return(list(beta.our=beta.our, beta.se.our=beta.se.our
              , censor.rate=censor.rate, cure.rate=cure.rate))
}


# parallel
no_cores <- detectCores(logical = FALSE)   # the number of core
cl <- makeCluster(no_cores)  # 
registerDoParallel(cl)  # registe parallel
MyRsult = foreach(r = 1:repetition, .combine = rbind, .multicombine = TRUE)  %dopar%  Repeat(r)
MyRsult = Transform(MyRsult)
MySummary = GetSummary()
stopImplicitCluster()  # end parallel













