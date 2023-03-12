#'@title Fit second order trend filtering/dynamic linear model via variational inference.
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
#'fit = ebtf2(x,v)

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
ebtf2 = function(x,v,w_init=NULL,m_init = NULL,sk=NULL,
                      maxiter=100,tol=1e-5,
                      calc_obj_every = 1,
                      printevery=Inf,
                      plot_every = Inf,
                      lambda0 = 100){


  n = length(x)
  if(length(v)==1){
    v = rep(v,n)
  }
  if(is.null(sk)){
    sk = choose_sk(x,v)
  }
  K = length(sk)
  vec1 = rep(1,n-3)
  S = outer(vec1,sk^2,"*")
  logS = log(S)

  # init
  if(is.null(m_init)){
    m = ti.thresh(x,sqrt(v))
  }else{
    m = m_init
  }


  mu = m[1]
  tau2 = 1
  alpha = matrix(nrow=n+3,ncol=K)
  alpha[2,] = -(m[2]^2-2*m[1]*m[2]+m[1]^2)/sk^2/2 - log(sk^2)/2
  alpha[3,] = -(m[3]^2-2*m[3]*(2*m[2]-m[1])+4*(m[2]^2)+m[1]^2-4*m[2]*m[1])/sk^2/2 - log(sk^2)/2
  for(i in 4:n){
    alpha[i,] = -(m[i]^2 - 2*m[i]*(3*m[i-1]-3*m[i-2]+m[i-3]) + 9*(m[i-1]^2-2*m[i-1]*m[i-2]+m[i-2]^2) + m[i-3]^2 + 6*(m[i-1]-m[i-2])*m[i-3])/sk^2/2-log(sk^2)/2
  }
  alpha = exp(alpha - rowMaxs(alpha,value = TRUE))
  alpha = alpha / rowsums(alpha)
  alpha = pmax(alpha,1e-10)
  #browser()
  if(is.null(w_init)){
    w = colSums(alpha[2:n,])
    w =w + c(lambda0,rep(0,K-1))
    w = w/sum(w)
  }else{
    w = w_init
  }

  logW = log(outer(vec1,w))

  obj_old = -Inf
  obj_vec = c()

  s2 = rep(0,n+3)
  m = c(m,0,0,0)
  alpha[1,] = 0
  alpha[n+1,] = 0
  alpha[n+2,] = 0
  alpha[n+3,] = 0

  start_time = Sys.time()
  for(iter in 1:maxiter){
    if(iter%%plot_every==0){
      plot(x,col='grey80')
      lines(m[1:n])
    }

    for(i in 1:n){
      if(i==1){
        s2[1] = 1/(1/v[1]+1/tau2+sum(alpha[2,]/sk^2)+sum(alpha[3,]/sk^2)+sum(alpha[4,]/sk^2))
        m[1] = s2[1]*(x[1]/v[1] + mu/tau2+sum(alpha[2,]/sk^2)*m[2] + sum(alpha[3,]/sk^2)*(2*m[2]-m[3]) + sum(alpha[4,]/sk^2)*(3*m[2]-3*m[3]+m[4]))
      }else if(i==2){
        s2[2] = 1/(1/v[2]+sum(alpha[2,]/sk^2)+4*sum(alpha[3,]/sk^2)+9*sum(alpha[4,]/sk^2)+sum(alpha[5,]/sk^2))
        m[2] = s2[2]*(x[2]/v[2] + sum(alpha[2,]/sk^2)*m[1]+2*sum(alpha[3,]/sk^2)*(m[1]+m[3]) + 3*sum(alpha[4,]/sk^2)*(m[1]+3*m[3] - m[4]) + sum(alpha[5,]/sk^2)*(3*m[3]-3*m[4]+m[5]))
        alpha[2,] = -(m[2]^2+s2[2]-2*m[1]*m[2]+m[1]^2+s2[1])/sk^2/2 - log(sk^2)/2 + log(w)
      }else if(i == 3){
        s2[3] = 1/(1/v[3]+sum(alpha[3,]/sk^2)+9*sum(alpha[4,]/sk^2)+9*sum(alpha[5,]/sk^2)+sum(alpha[6,]/sk^2))
        m[3] = s2[3]*(x[3]/v[3] + sum(alpha[3,]/sk^2)*(2*m[2]-m[1])+3*sum(alpha[4,]/sk^2)*(3*m[2]-m[1]+m[4])
                      + 3*sum(alpha[5,]/sk^2)*(m[2]+3*m[4] - m[5]) + sum(alpha[6,]/sk^2)*(3*m[4]-3*m[5]+m[6]))
        alpha[3,] = -(m[3]^2+s2[3]-2*m[3]*(2*m[2]-m[1])+4*(m[2]^2+s2[2])+m[1]^2+s2[1]-4*m[2]*m[1])/sk^2/2 - log(sk^2)/2 + log(w)
      }else{
        s2[i] = 1/(1/v[i]+sum(alpha[i,]/sk^2)+9*sum(alpha[i+1,]/sk^2)+9*sum(alpha[i+2,]/sk^2)+sum(alpha[i+3,]/sk^2))
        m[i] = s2[i]*(x[i]/v[i]+sum(alpha[i,]/sk^2)*(3*m[i-1]-3*m[i-2]+m[i-3]) + 3*sum(alpha[i+1,]/sk^2)*(3*m[i-1]-m[i-2]+m[i+1]) + 3*sum(alpha[i+2,]/sk^2)*(m[i-1]+3*m[i+1]-m[i+2]) + sum(alpha[i+3,]/sk^2)*(3*m[i+1]-3*m[i+2]+m[i+3]))
        alpha[i,] = -(m[i]^2+s2[i] - 2*m[i]*(3*m[i-1]-3*m[i-2]+m[i-3]) + 9*(m[i-1]^2+s2[i-1]-2*m[i-1]*m[i-2]+m[i-2]^2+s2[i-2]) + m[i-3]^2 + s2[i-3] + 6*(m[i-1]-m[i-2])*m[i-3])/sk^2/2-log(sk^2)/2+log(w)
      }
    }

    # alpha[2,] = -(m[2]^2+s2[2]-2*m[1]*m[2]+m[1]^2+s2[1])/sk^2/2 - log(sk^2)/2 + log(w)
    # alpha[3,] = -(m[3]^2+s2[3]-2*m[3]*(2*m[2]-m[1])+4*(m[2]^2+s2[2])+m[1]^2+s2[1]-4*m[2]*m[1])/sk^2/2 - log(sk^2)/2 + log(w)
    # for(i in 4:n){
    #   alpha[i,] = -(m[i]^2+s2[i] - 2*m[i]*(3*m[i-1]-3*m[i-2]+m[i-3]) + 9*(m[i-1]^2+s2[i-1]-2*m[i-1]*m[i-2]+m[i-2]^2+s2[i-2]) + m[i-3]^2 + s2[i-3] + 6*(m[i-1]-m[i-2])*m[i-3])/sk^2/2-log(sk^2)/2+log(w)
    # }


    alpha = exp(alpha - rowMaxs(alpha,value = TRUE))
    alpha = alpha / rowsums(alpha)
    alpha = pmax(alpha,1e-8)
    #browser()
    w = colsums(alpha[2:n,])
    w =w + c(lambda0,rep(0,K-1))
    w = w/sum(w)
    logW = log(outer(vec1,w))
    mu = m[1]
    tau2 = s2[1]

    alpha[1,] = 0
    alpha[n+1,] = 0
    alpha[n+2,] = 0
    alpha[n+3,] = 0

    if(iter%%calc_obj_every==0){
      obj_new = ebtf2_obj(x,v,m[1:n],s2[1:n],
                          m[3:(n-1)],m[2:(n-2)],m[1:(n-3)],
                          s2[3:(n-1)],s2[2:(n-2)],s2[1:(n-3)],alpha[1:n,],
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

ebtf2_obj = function(x,v,m,s2,
                     ml1,ml2,ml3,
                     s2_l1,s2_l2,s2_l3,alpha,
                     S,logS,logW,w,sk,mu,tau,lambda0){
  n = length(x)
  Eloglik =  - sum(log(v)/2) - sum((x^2-2*x*m+m^2+s2)/2/v)
  Elogprior = -log(tau)/2-(m[1]^2+s2[1]-2*m[1]*mu + mu^2)/2/tau + sum(alpha[2,]*(log(w)-log(sk^2)/2-(m[2]^2+s2[2]-2*m[2]*m[1]+m[1]^2+s2[1])/2/sk^2))
  + sum(alpha[3,]*(log(w)-log(sk^2)/2-(m[3]^2+s2[3]-2*m[3]*(2*m[2]-m[1])+4*(m[2]^2+s2[2]) + m[1]^2+s2[1] - 4*m[2]*m[1])/sk^2/2))
  +sum(alpha[4:n,]*(logW-logS/2-(m[4:n]^2+s2[4:n]-2*m[4:n]*(3*ml1-3*ml2+ml3)+9*(ml1^2+s2_l1-2*ml1*ml2 + ml2^2+s2_l2) + ml3^2+s2_l3 + 6*(ml1-ml2)*ml3))/S/2)
  Elogpost = sum(log(s2))/2 - sum(alpha[2:n,]*log(alpha[2:n,]))
  return(Eloglik+Elogprior+Elogpost+lambda0*log(w[1]))
}


