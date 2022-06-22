


#==================================================#
####               Sub functions                 ####
#==================================================#


KK = function(x,hn){
  # kernel funcition
  f = 3/4 * (1-(x/hn)^2) * (abs(x/hn)<=1)  # Epanechnikov kernel
  return(f)
}
KK2 = function(x,hn){
  # the indicater function for a binary covariate 
  f = hn^( ifelse(x!=0, 1,0) )
  return(f)
}
NWweight2 = function(Z, Z.sample, hn, lam=rep(1,n)){
  ## calculate the Nadaraya-Watson weights B_nk(Z) of dim m*k, Z is m*p, k=1,..,n of n observations of Z
  #
  # hn: bandwitdth dim is p*1
  # Z: the Argument Z in  B_nk(Z), Z could be 1*p (1 person) or m*p (m persons)
  # Z.sample: dim is n*p, the n observations of Z 
  # lam: the disturb value from bootstrap (to compute standar err), if do not computer SE, set lam = rep(1,n)
  # output B_nk(Z): dim is 1*k or m*k,k=1,...,n, the cols are (B_n1,B_nk(Z),...,B_nn(Z)), the i-th row is the i person in Z, i=1,...,m
  if(is.matrix(Z)){m = dim(Z)[1]; p1 = dim(Z)[2]  # Z is m*p using a matrix format
  }else{
    m = 1; p1 = length(Z)  # Z is 1*p using a vector format
    }
  n = dim(Z.sample)[1]
  Z = as.data.frame(Z); rownames(Z)= paste("person.",1:m,sep=""); colnames(Z)=paste("Z_",1:p1,sep="")
  Z.sample = as.data.frame(Z.sample); rownames(Z.sample)= paste("Obser",1:n,sep=""); colnames(Z.sample)=paste("Z_",1:p1,sep="")
  
  # Z-Z.k of dim m*k,k=1,...n
  MZ=array(apply(Z,2,function(x) rep(x,n)),dim=c(m,n,p1))  #  each fixed p1(Z is two dims, p1=2), dim is m*n in which each row denotes a person
  tMZ=aperm(array(apply(Z.sample,2,function(x) rep(x,m)),dim=c(n,m,p1)), c(2,1,3))   # each fixed p1(Z is two dims, p1=2), dim is m*n in which each row denotes a person, the k-th col denotes Z_k
  dMZ=MZ-tMZ # each fixed p1(Z is two dims, p1=2), dim is m*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
  
  Kp=array(NA,dim=c(m,n,p1))
  for(i in 1:p1){
    ## Kp[,,i] denotes for fixed Z.i, in the the n*n matrix, the k column denotes a person's Z minus the k person's Z
    # Kp[,,i]=KK(dMZ[,,i],hn[i])  # here, i deonts Z.i, i=1,2,...,p1'dim
    if(sum(dMZ[,,i]==1|dMZ[,,i]==-1|dMZ[,,i]==0)==(n*n) & any(dMZ[,,i]!=0)){
      Kp[,,i]=KK2(dMZ[,,i],hn[i]) # if Z is a binay covariate, use the Indicator function
    }
    else{
      Kp[,,i]=KK(dMZ[,,i],hn[i]) # if Z is a continous covariate, use a kernel function
    }
    
  }
  
  prodKp=apply(Kp,c(1,2),prod)  # K(x) = K(Z.1)*...* K(Z.p), since Z is p dims. K(x) is m*n.
  Lam=matrix(lam,nrow=m,ncol=n,byrow=T) 
  lamKp=prodKp*Lam  # each Z.k, times a disturbe value !!!
  B_nk=lamKp/matrix(apply(lamKp,1,sum),nrow=m,ncol=n,byrow=F)  # dim is m*n, the i row is the i person, the k column denbote Z-Z.k, namley a person's Z minus the k person's Z 
  return(B_nk)  
}
F.C = function(data, w){
  ## compute F.C([t,inifint|x) = product_i {1 - b.ni(x)/sum_k B.nk(x)I(Y.k>=Y.i)}^{I(Y.i<=t and delta.i=0)}, at t=Y.i Z=Z.j, i=1,...,n, j = 1,...,m
  #
  # data:
  # w: NWeight = B.nk(x), dim is m*k, k =1,...,n
  # output:  F.C([Y.i,inifint|Z.j), dim is m*n, the j row is the j person, the i col is F.C([Y.i,inifint|Z.j)
  
  m = nrow(w)  # m rows
  n = ncol(w)  # n cols
  if(m>1){k= cbind(rep(0,m),w[,1:(n-1)]) # add 0 at first col and delete the last col of B.nk(x)
  }else{k=t(c(rep(0,m), w[,1:(n-1)]))  }
  # k = cbind(rep(0,m), w[,1:(n-1)])  # add 0 at first col and delete the last col of B.nk(x)
  cumsum.k = t(apply(k, 1, cumsum))  # sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  cum.w = 1 - cumsum.k  # denominator: H([t,infinit)|x]) = sum_i B.ni(x)*I(Y.(i)>=t), the  1 col is t=Y.(1), the 2 col is t=Y.(2), ...., the n col is t=Y.(n) 
  # numerator: w = B.ni(x)
  w = replace(w, cum.w==0, 0)  # if the denominator is 0, then the numerator is 0
  cum.w  = replace(cum.w, cum.w==0, 1)  # if the denominator is 0, then the denominator is 1
  q = 1 - w/cum.w   # 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t), the j row is the j person, the 1 col is t=Y.(1), the i col is t=Y.(i), ...., the n col is t=Y.(n), dim is n*i
  delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  q = replace(q, delta==1, 1)  # if delta=1, then q = 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t) is q=1
  q = t(apply(q, 1, cumprod))  # F.c([t,infinit)|x) = product_i {1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t)}^{I(Y.i<=t and delta.i=0)}, dim is n*i, the j row is the j person, the i col is F.c([t,infinit)|x) at t=Y.(i)   
  
  
  return(q)
}
S.Ti = function(data, w){
  ## compute S.T(Y.i|x.i) where S.T(t|x) = product_i {1 - H1(ds|x)/(H([Y.i,infinit)|x)}^I(Y.i<=t,delta=1)  
  ## and compute f.T0(t|x) = F.T0(dt|x) = S.T0(Y.(i-1)|x.i) - S.T0(Y.i|x.i), at dt=Y.i and x=x.i
  #
  # w: NWeight = B.nk(x), dim is n*k, k =1,...,n
  #
  # output: 
  # S.T0(Y.i|x.i), i=1,...,n , dim is i*1
  # f.T0(t|x) = F.T0(dt|x) at dt=Y.i, x=x.i, d=1,...,n, dim is i*1
  #
  # Note: we can compute F.T0((t,infinit)|x) = S.T0(Y.i|x.j), the j row is the j person, the i col is F.T0((Y.i,inifint|x.j), dim is j*i
  # ..... but we only output S.T0(Y.i|x.i)
  
  # S.T0(Y.i|x.j)  dim is j*i, j=1,...,n, i=1,...,n
  m = nrow(w)
  n = ncol(w)
  k = cbind(rep(0,m), w[,1:(n-1)])  # add 0 at first col and delete the last col of B.nk(x)
  cumsum.k = t(apply(k, 1, cumsum))  # sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  cum.w = (1 - cumsum.k)   # H([t,infinit)|x), the j roe is the j person, the k col is this value H-P()*F.c at t=Y.(k)   
  cum.w = cum.w*(cum.w>w+0.00001) + (w+0.0001)*(cum.w<=w+0.00001)  # ??? if the dinominator<numerator, then dinominator<numerato, then 1 - dinominator/nummerator = 0, then S.T0(t|x)=0
  q = 1 - w/cum.w  # q = 1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x)), the j row is the j person, the i col is this value at t=Y.i 
  delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  q = replace(q, delta==0, 1)  # if delta=0, then q=1
  q = t(apply(q, 1, cumprod))  # F.T0((t,infinit)|x) = S.T0(Y.i|x.j) = product_i {1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x))}, the j row is the j person, the i col is  F.T0((t,infinit)|x) at t=Y.(i)  
  
  # S.T0(Y.i|x.i)  dim is i*1, i=1,...,n
  S.T0 = diag(q) # F.T0((t,infinit)|x) = S.T0(Y.i|x.i)
  
  return(S.T0)
  
  ## f.T0(Y.i|x.i) dim is i*1, i=1,...,n
  # qq = cbind(rep(1,n),q)  ## F.T0((t,infinit)|x), and add the value S(0|x)=1, each row denotes each person, each column denotes S(Y.(i)|x)
  # qqq = c(t(qq))  # all n persons to be a vector
  # r = apply(cbind(rep(n,n),rank(data[,1])),1,min)  # 1,2,...n
  # s = seq(0,n-1,1)  
  # S.Yi1 = qqq[(n+1)*s+r]  # S.T0(Y.(i-1)|x.(i))
  # S.Yi = qqq[(n+1)*s+r+1]  # S.T0(Y.(i)|x.(i))
  # f.T0 = S.Yi1 - S.Yi  # f.T0(Y.i|x.i) = S.T0(Y.(i-1)|x.i) - S.T0(Y.i|x.i), i=1,...,n, dim is i*1
  # return(list(S.T0=S.T0, f.T0=f.T0))
  
}
P.T = function(data, w){
  ## compute P(T=inifint|x.j) = product_i {1 - b.ni(x)/sum_k B.nk(x)I(Y.k>=t)}^{I(delta.i=1)}, at t=infinit=Y.n, x=x.j, j=1,...,n
  ## namely  P(T=infinity|z) = 1 - P(eta=1|z) = 1-pi(z) = S.T(t|x) at t>=Y.(n)
  #
  # data:
  # w: NWeight = B.nk(x), dim is n*k, k =1,...,n
  # output: P(T=inifint|x.j), dim is j*1 , j=1,...,n
  
  m = nrow(w)  # m rows
  n = ncol(w)  # n cols
  k = cbind(rep(0,m), w[,1:(n-1)])  # add 0 at first col and delete the last col of B.nk(x)
  cumsum.k = t(apply(k, 1, cumsum))  # sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  cum.w = 1 - cumsum.k  # denominator: H([t,infinit)|x]) = sum_i B.ni(x)*I(Y.(i)>=t), the  1 col is t=Y.(1), the 2 col is t=Y.(2), ...., the n col is t=Y.(n) 
  # numerator: w = B.ni(x)
  w = replace(w, cum.w==0, 0)  # if the denominator is 0, then the numerator is 0
  cum.w = replace(cum.w, cum.w==0, 1)  # if the denominator is 0, then the denominator is 1
  q = 1 - w/cum.w  # 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t), the j row is the j person, the 1 col is t=Y.(1), the i col is t=Y.(i), ...., the n col is t=Y.(n), dim is n*i 
  delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  q = replace(q, delta==0, 1)   # if delta=0, then q = 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t) is q=1
  q = t(apply(q, 1, prod))  # P(T=infinit|x.j) = S(Y.(n)|x.j) = product_i {1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t)} at x=x.j, dim is n*1, the j row is the j person at S(Y.(n)|x.j)
  
  
  return(as.numeric(q))
}
## for Onlie data using bootstrap method of using array to speed up and avoid 'for' circle
NWweight2.array = function(Z, Z.sample, hn, lam){
  ## calculate the Nadaraya-Watson weights B_nk(Z) of dim m*k*B, Z is m*p, k=1,..,n of n observations of Z
  #
  # hn: bandwitdth dim is p*1
  # Z: the Argument Z in  B_nk(Z), Z could be 1*p (1 person) or m*p (m persons)
  # Z.sample: dim is n*p, the n observations of Z 
  # lam: the disturb value from bootstrap (to compute standar err), dim is n*B
  # output B_nk(Z): dim is m*k*B, k=1,...,n, b=1,...,B, for fixed b, the cols are (B_n1,B_nk(Z),...,B_nn(Z)), the i-th row is the i person in Z, i=1,...,m
  
  
  B = ncol(lam)  # the bootstrap times B
  if(is.matrix(Z)){m = dim(Z)[1]; p1 = dim(Z)[2]  # Z is m*p using a matrix format
  }else{
    m = 1; p1 = length(Z)  # Z is 1*p using a vector format
  }
  n = dim(Z.sample)[1]
  Z = as.data.frame(Z); rownames(Z)= paste("person.",1:m,sep=""); colnames(Z)=paste("Z_",1:p1,sep="")
  Z.sample = as.data.frame(Z.sample); rownames(Z.sample)= paste("Obser",1:n,sep=""); colnames(Z.sample)=paste("Z_",1:p1,sep="")
  
  # Z-Z.k of dim m*k,k=1,...n
  MZ=array(apply(Z,2,function(x) rep(x,n)),dim=c(m,n,p1))  #  each fixed p1(Z is two dims, p1=2), dim is m*n in which each row denotes a person
  tMZ=aperm(array(apply(Z.sample,2,function(x) rep(x,m)),dim=c(n,m,p1)), c(2,1,3))   # each fixed p1(Z is two dims, p1=2), dim is m*n in which each row denotes a person, the k-th col denotes Z_k
  dMZ=MZ-tMZ # each fixed p1(Z is two dims, p1=2), dim is m*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
  
  Kp=array(NA,dim=c(m,n,p1))
  for(i in 1:p1){
    ## Kp[,,i] denotes for fixed Z.i, in the the n*n matrix, the k column denotes a person's Z minus the k person's Z
    # Kp[,,i]=KK(dMZ[,,i],hn[i])  # here, i deonts Z.i, i=1,2,...,p1'dim
    if(sum(dMZ[,,i]==1|dMZ[,,i]==-1|dMZ[,,i]==0)==(n*n) & any(dMZ[,,i]!=0)){
      Kp[,,i]=KK2(dMZ[,,i],hn[i]) # if Z is a binay covariate, use the Indicator function
    }
    else{
      Kp[,,i]=KK(dMZ[,,i],hn[i]) # if Z is a continous covariate, use a kernel function
    }
    
  }
  
  prodKp=apply(Kp,c(1,2),prod)  # K(x) = K(Z.1)*...* K(Z.p), since Z is p dims. K(x) is m*n.
  # Lam=matrix(lam, nrow=m, ncol=n, byrow=T) 
  # lamKp=prodKp*Lam  # each Z.k, times a disturbe value !!!, dim is m*n
  # B_nk=lamKp/matrix(apply(lamKp,1,sum),nrow=m,ncol=n,byrow=F)  # dim is m*n, the i row is the i person, the k column denbote Z-Z.k, namley a person's Z minus the k person's Z 
  Lam = array(apply(lam, 2, function(x) outer(rep(1,m), x)), dim = c(m,n,B))  # for fixed b, Lam is dim of m*n,  b=1,...B
  lamKp = array(apply(Lam, 3, function(x) x*prodKp), dim = c(m,n,B))  #  for fixed b, lamKp is dim of m*n, lamKp = prodKp*Lam
  temp = apply(lamKp, 3, function(x) x%*%rep(1,n))
  temp2 = apply(array(apply(temp, 2, function(x) outer(rep(1,n), x)), dim = c(n,m,B)), 3, t)
  fenmu = array(temp2, dim = c(m,n,B))  # sum_k=1^n epsilon.k*K(x-x.k)
  B_nk = lamKp / fenmu  # dim is m*n*b, for fix b,  the i row is the i person, the k column denbote Z-Z.k, namley a person's Z minus the k person's Z 
  return(B_nk)  
}
S.Ti.array = function(data, w){
  ## compute S.T(Y.i|x.i) where S.T(t|x) = product_i {1 - H1(ds|x)/(H([Y.i,infinit)|x)}^I(Y.i<=t,delta=1)  
  ## and compute f.T0(t|x) = F.T0(dt|x) = S.T0(Y.(i-1)|x.i) - S.T0(Y.i|x.i), at dt=Y.i and x=x.i
  #
  # w: NWeight = B.nk(x), dim is m*k*B, k =1,...,n
  #
  # output: 
  # S.T0(Y.i|x.i), i=1,...,n , dim is i*B
  #
  # Note: we can compute F.T0((t,infinit)|x) = S.T0(Y.i|x.j), the j row is the j person, the i col is F.T0((Y.i,inifint|x.j), dim is j*i
  # ..... but we only output S.T0(Y.i|x.i)
  
  
  B = dim(w)[3]
  m = dim(w)[1]
  n = dim(w)[2]
  
  # k = cbind(rep(0,m), w[,1:(n-1)])  # add 0 at first col and delete the last col of B.nk(x)
  # cumsum.k = t(apply(k, 1, cumsum))  # sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  # cum.w = (1 - cumsum.k)   # H([t,infinit)|x), the j roe is the j person, the k col is this value H-P()*F.c at t=Y.(k)   
  # cum.w = cum.w*(cum.w>w+0.00001) + (w+0.0001)*(cum.w<=w+0.00001)  # ??? if the dinominator<numerator, then dinominator<numerato, then 1 - dinominator/nummerator = 0, then S.T0(t|x)=0
  # q = 1 - w/cum.w  # q = 1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x)), the j row is the j person, the i col is this value at t=Y.i 
  # delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  # q = replace(q, delta==0, 1)  # if delta=0, then q=1
  # q = t(apply(q, 1, cumprod))  # F.T0((t,infinit)|x) = S.T0(Y.i|x.j) = product_i {1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x))}, the j row is the j person, the i col is  F.T0((t,infinit)|x) at t=Y.(i)  
  # 
  
  k = array(apply(w, 3, function(x) cbind(0, x[,1:(n-1)])), dim = c(m,n,B))  # for fixed b, add 0 at first col and delete the last col of B.nk(x)
  cumsum.k = array(apply(k, 3, function(x) t(apply(x, 1, cumsum))), dim = c(m,n,B))  # for fixed b, sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  cum.w = 1 - cumsum.k
  cum.w = cum.w*(cum.w>w+0.00001) + (w+0.0001)*(cum.w<=w+0.00001)  # ??? if the dinominator<numerator, then dinominator<numerato, then 1 - dinominator/nummerator = 0, then S.T0(t|x)=0
  q = 1 - w/cum.w  # q = 1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x)), the j row is the j person, the i col is this value at t=Y.i 
  delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  q = array(apply(q, 3, function(x) replace(x, delta==0, 1)), dim = c(m,n,B))
  q = array(apply(q, 3, function(x) t(apply(x, 1, cumprod))), dim = c(m,n,B))
  
  # S.T0(Y.i|x.i)  dim is i*1, i=1,...,n
  # S.T0 = diag(q) # F.T0((t,infinit)|x) = S.T0(Y.i|x.i)
  S.T0 = apply(q, 3, function(x) diag(x))  # dim is n*B
  
  return(S.T0)
}
P.T.array = function(data, w){
  ## compute P(T=inifint|x.j) = product_i {1 - b.ni(x)/sum_k B.nk(x)I(Y.k>=t)}^{I(delta.i=1)}, at t=infinit=Y.n, x=x.j, j=1,...,n
  ## namely  P(T=infinity|z) = 1 - P(eta=1|z) = 1-pi(z) = S.T(t|x) at t>=Y.(n)
  #
  # data:
  # w: NWeight = B.nk(x), dim is m*k*B, k =1,...,n
  # output: P(T=inifint|x.j), dim is j*B , j=1,...,n
  
  # m = nrow(w)  # m rows
  # n = ncol(w)  # n cols
  # k = cbind(rep(0,m), w[,1:(n-1)])  # add 0 at first col and delete the last col of B.nk(x)
  # cumsum.k = t(apply(k, 1, cumsum))  # sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  # cum.w = 1 - cumsum.k  # denominator: H([t,infinit)|x]) = sum_i B.ni(x)*I(Y.(i)>=t), the  1 col is t=Y.(1), the 2 col is t=Y.(2), ...., the n col is t=Y.(n) 
  # # numerator: w = B.ni(x)
  # w = replace(w, cum.w==0, 0)  # if the denominator is 0, then the numerator is 0
  # cum.w = replace(cum.w, cum.w==0, 1)  # if the denominator is 0, then the denominator is 1
  # q = 1 - w/cum.w  # 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t), the j row is the j person, the 1 col is t=Y.(1), the i col is t=Y.(i), ...., the n col is t=Y.(n), dim is n*i 
  # delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  # q = replace(q, delta==0, 1)   # if delta=0, then q = 1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t) is q=1
  # q = t(apply(q, 1, prod))  # P(T=infinit|x.j) = S(Y.(n)|x.j) = product_i {1 - B.ni(x)/ sum_k B.nk(x)*I(Y.(k)>=t)} at x=x.j, dim is n*1, the j row is the j person at S(Y.(n)|x.j)
  # 
  
  B = dim(w)[3]
  m = dim(w)[1]
  n = dim(w)[2]
  k = array(apply(w, 3, function(x) cbind(0, x[,1:(n-1)])), dim = c(m,n,B))  # for fixed b, add 0 at first col and delete the last col of B.nk(x)
  cumsum.k = array(apply(k, 3, function(x) t(apply(x, 1, cumsum))), dim = c(m,n,B))  # for fixed b, sum_i B.ni(x)*I(Y.(i)<=t), the 1 col is t=Y.(0)=0, the 2 col is t=Y.(1), ...., the n col is t=Y.(n-1) 
  cum.w = 1 - cumsum.k
  cum.w = cum.w*(cum.w>w+0.00001) + (w+0.0001)*(cum.w<=w+0.00001)  # ??? if the dinominator<numerator, then dinominator<numerato, then 1 - dinominator/nummerator = 0, then S.T0(t|x)=0
  q = 1 - w/cum.w  # q = 1 - H1(dt|x)/(H([t,infinit)|x)-P(T=infinit|x)*F.C([t,infinit)|x)), the j row is the j person, the i col is this value at t=Y.i 
  delta = matrix(data$delta, nrow = m, ncol = n, byrow = TRUE)
  q = array(apply(q, 3, function(x) replace(x, delta==0, 1)), dim = c(m,n,B))
  q = apply(q, 3, function(x) t(apply(x, 1, prod)))  # dim is n*B
  
  return(q)
}




