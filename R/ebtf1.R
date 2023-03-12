#'@title Fit first order trend filtering/dynamic linear model via variational inference.
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
#'fit = ebtf1(x,v)

#'plot(x,col='grey80')
#'lines(mu,col='grey50')
#'lines(fit$posterior$mean)

#'@importFrom smashr smash.gaus
#'@importFrom smashr ti.thresh
#'@importFrom wavethresh DJ.EX
#'@importFrom Rfast rowsums
#'@importFrom Rfast colsums
#'@importFrom Rfast rowMaxs
#'@export

ebtf1 = function(x,v,
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
  vec1 = rep(1,n-2)
  S = outer(vec1,sk^2,"*")
  logS = log(S)
  if(is.null(m_init)){
    m = ti.thresh(x,sqrt(v))
  }else{
    m = m_init
  }


  #s2 = c()
  x1 = x[1]
  x2 = x[2]
  x3n = x[3:n]
  v1 = v[1]
  v2 = v[2]
  v3n = v[3:n]

  mu = m[1]
  tau2 = 1
  m1 = m[1]
  m2 = m[2]
  m3n = m[3:n]

  ml1 = c(m2,m3n[1:(n-3)])
  ml2 = c(m1,m2,m3n[1:(n-4)])
  mr1 = c(m3n[-1],0)
  mr2 = c(mr1[-1],0)


  if(is.null(w_init)){
    if(is.null(m_init)){
      w  = rep(1/K,K)
      alpha3n = matrix(1,n-2,K)/K
      alpha2 = w
    }else{
      alpha2 = -(m2^2-2*m2*m1+m1^2)/sk^2/2-log(sk^2)/2
      alpha2 = exp(alpha2 - max(alpha2))
      alpha2 = alpha2 / sum(alpha2)
      alpha3n = -((m3n^2)-2*m3n*(2*ml1-ml2)+4*(ml1^2) + (ml2^2) - 4*ml1*ml2)/S/2-logS/2
      alpha3n = exp(alpha3n - rowMaxs(alpha3n,value = TRUE))
      alpha3n = alpha3n / rowsums(alpha3n)

      w = colSums(alpha3n)
      w =w + alpha2 + c(lambda0,rep(0,K-1))
      w = w/sum(w)
    }
  }else{
    w = w_init
    alpha3n = matrix(rep(w,each=n-2),nrow=n-2)
    alpha2 = w
  }


  obj_old = -Inf
  obj_vec = c()
  logW = log(outer(vec1,w)+.Machine$double.eps)

  s2 = rep(0,n+2)
  m = c(m,0,0)
  alpha = rbind(rep(0,K),alpha2,alpha3n,rep(0,K),rep(0,K))
  start_time = Sys.time()
  for(iter in 1:maxiter){

    if(iter%%plot_every==0){
      plot(x,col='grey80')
      lines(m[1:n])
    }

    for(i in 1:n){
      if(i==1){
        s2[1] = 1/(1/v1+1/tau2+sum(alpha[2,]/sk^2) + sum(alpha[3,]/sk^2))
        m[1] = s2[1]*(x1/v1+mu/tau2 + m[2]*sum(alpha[2,]/sk^2) + (2*m[2]-m[3])*sum(alpha[3,]/sk^2))
      }else if(i==2){
        s2[2] = 1/(1/v1+sum(alpha[2,]/sk^2)+4*sum(alpha[3,]/sk^2)+sum(alpha[4,]/sk^2))
        m[2] = s2[2]*(x2/v2+m[1]*sum(alpha[2,]/sk^2)+2*(m[1]+m[3])*sum(alpha[3,]/sk^2) + (2*m[3]-m[4])*sum(alpha[4,]/sk^2))##
        alpha[2,] = -(m[2]^2+s2[2]-2*m[2]*m[1]+m[1]^2+s2[1])/sk^2/2-log(sk^2)/2+log(w)#
      }else{
        s2[i] = 1/(1/v[i]+sum(alpha[i,]/sk^2)+4*sum(alpha[i+1,]/sk^2)+sum(alpha[i+2,]/sk^2))
        m[i] = s2[i]*(x[i]/v[i]+(2*m[i-1]-m[i-2])*sum(alpha[i,]/sk^2)+2*(m[i-1]+m[i+1])*sum(alpha[i+1,]/sk^2) + sum(alpha[i+2,]/sk^2)*(2*m[i+1]-m[i+2]))
        alpha[i,] = -((m[i]^2+s2[i])-2*m[i]*(2*m[i-1]-m[i-2])+4*(m[i-1]^2+s2[i-1]) + (m[i-2]^2+s2[i-2]) - 4*m[i-1]*m[i-2])/sk^2/2-log(sk^2)/2+log(w)
      }
    }

    # alpha[2,] = -(m[2]^2+s2[2]-2*m[2]*m[1]+m[1]^2+s2[1])/sk^2/2-log(sk^2)/2+log(w)
    # for(i in 3:n){
    #   alpha[i,] = -((m[i]^2+s2[i])-2*m[i]*(2*m[i-1]-m[i-2])+4*(m[i-1]^2+s2[i-1]) + (m[i-2]^2+s2[i-2]) - 4*m[i-1]*m[i-2])/sk^2/2-log(sk^2)/2+log(w)
    # }

    #update alpha
    alpha = exp(alpha - rowMaxs(alpha,value = TRUE))
    alpha = alpha / rowsums(alpha)
    alpha = pmax(alpha,1e-10)
    #browser()
    w = colSums(alpha[2:n,])
    w =w + c(lambda0,rep(0,K-1))
    w = w/sum(w)
    logW = log(outer(vec1,w))
    mu = m[1]
    tau2 = s2[1]

    alpha[1,] = 0
    alpha[n+1,] = 0
    alpha[n+2,] = 0

    if(iter%%calc_obj_every==0){
      obj_new = ebtf1_obj(x1,x2,x3n,v1,v2,v3n,m[1],m[2],m[3:n],
                          s2[1],s2[2],s2[3:n],
                          m[2:(n-1)],m[1:(n-2)],
                          s2[2:(n-1)],s2[1:(n-2)],alpha[2,],alpha[3:n,],
                          S,logS,logW,w,sk,mu,tau2,lambda0)
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
  return(list(posterior=list(mean = m[1:n],sd=sqrt(s2[1:n])),
              fitted_g=list(w=w,mean=rep(0,K),sd=sk),
         run_time = difftime(Sys.time(),start_time,units = 'secs')))
}

ebtf1_obj = function(x1,x2,x3n,v1,v2,v3n,m1,m2,m3n,
                     s2_1,s2_2,s2_3n,
                     ml1,ml2,
                     s2_l1,s2_l2,alpha2,alpha3n,
                     S,logS,logW,w,sk,mu,tau,lambda0){
  #n = length(x)
  Eloglik = (-log(v1)/2-(x1^2-2*x1*m1+m1^2+s2_1)/2/v1 -log(v2)/2-(x2^2-2*x2*m2+m2^2+s2_2)/2/v2 - sum(log(v3n)/2+(x3n^2-2*x3n*m3n+m3n^2+s2_3n)/2/v3n))
  Elogprior = -log(tau)/2-(m1^2+s2_1-2*m1*mu + mu^2)/2/tau + sum(alpha2*(-log(sk^2)/2-(m2^2+s2_2-2*m2*m1+m1^2+s2_1)/2/sk^2)) + sum(alpha2*log(w))
  -sum(alpha3n*logS)/2-1/2*sum((m3n^2+s2_3n-2*m3n*(2*ml1-ml2)+4*(ml1^2+s2_l1) + (ml2^2+s2_l2) - 4*ml1*ml2)*alpha3n/S) + sum(alpha3n*logW)
  Elogpost = (sum(log(s2_3n)) + log(s2_1) + log(s2_2))/2 - sum(alpha2*log(alpha2) - sum(alpha3n*log(alpha3n)))
  return(Eloglik+Elogprior+Elogpost+lambda0*log(w[1]))
}





