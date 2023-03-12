#'@title Fit empirical Bayes trend filtering/dynamic linear model via variational inference.
#'@description Assuming the posterior factorizes as q(b,z) = q(b)q(z)
#'@param x,v data vector and variance
#'@param ord the order of trend filtering, currently can be 0,1, or 2.
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
#'fit = ebtf(x,v)

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
#'
ebtf = function(x,v,
                ord = 0,
                sk = NULL,
                w_init=NULL,
                m_init=NULL,
                maxiter = 500,
                tol=1e-5,
                printevery=Inf,
                calc_obj_every = 1,
                plot_every = Inf,
                lambda0 = 100){
  if(ord==0){
    ebtf0(x=x,v=v,
          sk=sk,
          w_init=w_init,
          m_init=m_init,
          maxiter=maxiter,
          tol=tol,
          printevery=printevery,
          calc_obj_every=calc_obj_every,
          plot_every=plot_every,
          lambda0=lambda0)
  }else if(ord==1){
    ebtf1(x=x,v=v,
          sk=sk,
          w_init=w_init,
          m_init=m_init,
          maxiter=maxiter,
          tol=tol,
          printevery=printevery,
          calc_obj_every=calc_obj_every,
          plot_every=plot_every,
          lambda0=lambda0)
  }else if(ord==2){
    ebtf2(x=x,v=v,
          sk=sk,
          w_init=w_init,
          m_init=m_init,
          maxiter=maxiter,
          tol=tol,
          printevery=printevery,
          calc_obj_every=calc_obj_every,
          plot_every=plot_every,
          lambda0=lambda0)
  }else{
    stop('choose ord from 0,1,2')
  }
}
