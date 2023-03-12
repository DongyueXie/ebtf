#'@title Fit zeroth order trend filtering/dynamic linear model via variational inference.
#'@description Assuming q(b,z) = q(b)q(z)
#'@param x,v data vector and variance
#'@param sk grid of prior standard deviations
#'@param w_init init of prior weights
#'@param m_init init of smoothed curve/posterior mean
#'@param maxiter,tol criteria to stop the iteration
#'@param printevery how often to print results
#'@param calc_obj_every how often to calculate elbo
#'@param plot_every how often to plot results
#'@param lambda0 The parameter promoting further sparsity/smoothness
#'@return A list of: posterior, fitted_g
#'@examples
#'library(ebtf)
#'set.seed(123)
#'n=2^9
#'xx = DJ.EX(n,rsnr = 3)
#'mu = xx$blocks
#'v = rep(1,n)
#'v[sample(1:n,round(n*0.1))] = 10
#'v[sample(1:n,round(n*0.1))] = 0.1
#'x = mu + rnorm(n,0,sqrt(v))
#'
#'fit = ebtf0(x,v)

#'plot(x,col='grey80')
#'lines(mu,col='grey50')
#'lines(fit$posterior$mean)

#'rmse = function(x,y){sqrt(mean((x-y)^2))}
#'rmse(mu,smashr::ti.thresh(x,sqrt(v)))
#'rmse(mu,smashr::smash.gaus(x,sqrt(v)))
#'rmse(mu,fit$posterior$mean)
#'
#'@importFrom smashr smash.gaus
#'@importFrom smashr ti.thresh
#'@importFrom wavethresh DJ.EX
#'@importFrom Rfast rowsums
#'@importFrom Rfast colsums
#'@importFrom Rfast rowMaxs
#'@export

ebtf0 = function(x,v,
                 sk = NULL,
                 w_init=NULL,
                 m_init=NULL,
                 maxiter = 100,
                 tol=1e-8,
                 printevery=Inf,
                 calc_obj_every = 1,
                 plot_every = Inf,
                 lambda0 = 100){
  # This function updates m,s2 and alpha using matrix operation and as blocks


  n = length(x)
  if(length(v)==1){
    v = rep(v,n)
  }
  if(is.null(sk)){
    sk = choose_sk(x,v)
  }
  K = length(sk)
  vec1 = rep(1,n-1)
  S = outer(vec1,sk^2,"*")
  logS = log(S)
  if(is.null(m_init)){
    m = ti.thresh(x,sqrt(v))
  }else{
    m = m_init
  }



  if(is.null(w_init)){
    alpha2n = -(m[-n]^2+m[-1]^2-2*m[-n]*m[-1])/S/2-logS/2
    alpha2n = exp(alpha2n - rowMaxs(alpha2n,value = TRUE))
    alpha2n = alpha2n / rowsums(alpha2n)
    w = colsums(alpha2n) + c(lambda0,rep(0,K-1))
    w = w/sum(w)
  }else{
    w = w_init
    alpha2n = matrix(rep(w,each=n-1),nrow=n-1)
  }

  #s2 = c()
  x1 = x[1]
  x2n = x[2:n]
  v2n = v[2:n]
  v1 = v[1]
  mu = m[1]
  tau2 = 1
  m1 = m[1]
  m2n = m[2:n]
  ml = c(m1,m2n[1:(n-2)])
  mr = c(m2n[-1],0)
  obj_old = -Inf
  obj_vec = c()
  logW = log(outer(vec1,w)+.Machine$double.eps)

  start_time = Sys.time()
  for(iter in 1:maxiter){

    if(iter%%plot_every==0){
      plot(x,col='grey80')
      lines(c(m1,m2n))
    }

    temp0 = sum(alpha2n[1,]/sk^2)
    temp1 = c(rowSums(alpha2n/S))
    temp2 = c(temp1[-1], 0)

    s2_1 = 1/(1/v1+1/tau2+temp0)
    s2_2n = 1/(1/v2n+temp1+temp2)
    s2_l = c(s2_1,s2_2n[1:(n-2)])

    m1 = s2_1*(x1/v1+mu/tau2+temp0*m2n[1])
    m2n = s2_2n*(x2n/v2n+ml*temp1+mr*temp2)

    ml = c(m1,m2n[1:(n-2)])
    mr = c(m2n[-1],0)

    #update alpha
    alpha2n = -(m2n^2+s2_2n+ml^2+s2_l-2*m2n*ml)/S/2-logS/2+logW
    alpha2n = exp(alpha2n - rowMaxs(alpha2n,value = TRUE))
    alpha2n = alpha2n / rowsums(alpha2n)
    alpha2n = pmax(alpha2n,1e-10)

    w = colsums(alpha2n) + c(lambda0,rep(0,K-1))
    w = w/sum(w)
    logW = log(outer(vec1,w)+.Machine$double.eps)
    mu = m1
    tau2 = s2_1
    if(iter%%calc_obj_every==0){
      obj_new = ebtf0_obj(x1,x2n,v1,v2n,m1,m2n,s2_1,s2_2n,ml,s2_l,alpha2n,S,logS,logW,w,lambda0)
      obj_vec[iter/calc_obj_every] = obj_new
      if(abs(obj_new-obj_old)/n<tol){
        break
      }
      obj_old = obj_new
    }
    if(iter%%printevery==0){
      print(paste("Done iter",iter,"obj =",obj_new))
    }
  }
  return(list(posterior=list(mean = c(m1,m2n), sd = sqrt(c(s2_1,s2_2n)),alpha=alpha2n),
              fitted_g = list(w=w,mean=rep(0,K),sd=sk),
              run_time = difftime(Sys.time(),start_time,units = 'secs')))
}

ebtf0_obj = function(x1,x2n,v1,v2n,m1,m2n,s2_1,s2_2n,ml,s2_l,alpha,S,logS,logW,w,lambda0){
  #n = length(x)
  return(sum(-1/2/v2n*(m2n^2+s2_2n-2*x2n*m2n))-1/2/v1*(m1^2+s2_1-2*x1*m1)+
           sum(alpha*logW)-sum(alpha*logS)/2-
           1/2*sum((m2n^2+s2_2n-2*m2n*ml+ml^2+s2_l)*alpha/S)+
           1/2*sum(log(s2_2n))-sum(alpha*log(alpha+.Machine$double.eps))+lambda0*log(w[1]))
}

#'@title choose default grid
choose_sk = function(x,v,grid_mult = sqrt(2)){
  n = length(x)
  z = x[-1]-x[-n]
  sz = sqrt(v[-1] + v[-n])
  smin = min(sz)/10
  smax = 2*sqrt(max(z^2-sz^2))
  sk = exp(seq(log(smin),log(smax),by=log(grid_mult)))
  if(sk[1]>1e-4){
    sk = c(1e-4,sk)
  }
  sk
}





