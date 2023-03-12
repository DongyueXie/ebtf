library(genlasso)
library(btf)
library(spmrf)
library(rstan)
library(ebtf)
library(genlasso)
rmse = function(x,y){sqrt(mean((x-y)^2))}
mae = function(x,y){mean(abs(x-y))}
simu_study = function(mu,snr=3,nreps=20,seed=12345,ord=0){
  set.seed(seed)
  errors = c()
  runtimes = c()
  sigma2 = var(mu)/snr
  n = length(mu)
  for (i in 1:nreps) {
    print(paste('Running ',i))
    runtime_i = c()
    y = mu + rnorm(n,0,sqrt(sigma2))


    t1_tf=Sys.time()
    fit_tf = trendfilter(y,ord=ord)
    fit_tf_cv = cv.trendfilter(fit_tf)
    fit_tf = fit_tf$beta[,fit_tf$lambda==fit_tf_cv$lambda.min]
    t2_tf=Sys.time()
    runtime_i = c(runtime_i,difftime(t2_tf,t1_tf,units = 'secs'))

    t1_btf=Sys.time()
    fit_btf = btf(y,k=ord)
    t2_btf=Sys.time()
    runtime_i = c(runtime_i,difftime(t2_btf,t1_btf,units = 'secs'))

    fit_ebtf0 = ebtf(y,sigma2,ord = 0,m_init = fit_tf)
    fit_ebtf1 = ebtf(y,sigma2,ord = 1,m_init = fit_tf)
    fit_ebtf2 = ebtf(y,sigma2,ord = 2,m_init = fit_tf)
    runtime_i = c(runtime_i,c(fit_ebtf0$run_time,fit_ebtf1$run_time,fit_ebtf2$run_time))

    # t1_spmrf=Sys.time()
    # spmrf_draw = spmrf(prior="horseshoe", likelihood="normal", order=ord+1,  data=list(y=y))
    # fit_spmrf_ci = extract_theta(spmrf_draw,obstype = 'normal',alpha=0.05)
    # fit_spmrf = fit_spmrf_ci$postmed
    # t2_spmrf=Sys.time()
    # runtime_i = c(runtime_i,difftime(t2_spmrf,t1_spmrf,units = 'secs'))

    # t1_smash = Sys.time()
    # fit_smash = smash.gaus(y,sqrt(sigma2))
    # t2_smash = Sys.time()
    # runtime_i = c(runtime_i,difftime(t2_smash,t1_smash,units = 'secs'))



    runtimes = rbind(runtimes,runtime_i)
    errors = rbind(errors,c(rmse(mu,fit_ebtf0$posterior$mean),
                            rmse(mu,fit_ebtf1$posterior$mean),
                            rmse(mu,fit_ebtf2$posterior$mean),
                            rmse(mu,fit_tf),
                            rmse(mu,apply(fit_btf$f,2,mean))
                            #rmse(mu,fit_smash)
                            ))
    maes = rbind(errors,c(mae(mu,fit_ebtf0$posterior$mean),
                          mae(mu,fit_ebtf1$posterior$mean),
                          mae(mu,fit_ebtf2$posterior$mean),
                          mae(mu,fit_tf),
                          mae(mu,apply(fit_btf$f,2,mean))
                          #mae(mu,fit_smash)
                          ))

  }
  colnames(runtimes) = c('ebtf0','ebtf1','ebtf2','cv.tf','btf','smash')
  colnames(errors) = c('ebtf0','ebtf1','ebtf2','cv.tf','btf','smash')
  colnames(maes) = c('ebtf0','ebtf1','ebtf2','cv.tf','btf','smash')
  return(list(runtimes=runtimes,rmses = errors,maes = maes))
}


# block structure

n=2^9
xx = DJ.EX(n,signal = sqrt(3))
mu = xx$blocks
blocks_res = simu_study(mu,nreps = 20)

# mu = xx$heavi
# heavi_res = simu_study(mu,nreps = 20,ord=2)

mu = smashrgen:::angles.f((1:n)/n)
mu = mu*sqrt(3/var(mu))
angle_res = simu_study(mu,nreps = 20,ord=1)

# mu = smashrgen:::sblocks.f((1:n)/n)
# mu = mu*sqrt(3/var(mu))
# sblocks_res = simu_study(mu,nreps = 20,ord=0)
n=2^10
xx = DJ.EX(n,signal = sqrt(3))
mu = xx$blocks
plot(xx$blocks,type='l',xlab='',ylab='')
sigma2 = 1
y = mu + rnorm(n,0,sqrt(sigma2))

# fit_tf = trendfilter(y,ord=0)
# fit_tf_cv = cv.trendfilter(fit_tf)
# fit_tf = fit_tf$beta[,fit_tf$lambda==fit_tf_cv$lambda.min]
# fit_btf = btf(y,k=0)
# fit_ebtf = ebtf(y,sigma2,ord=0,tol=1e-8,m_init=fit_tf,lambda0 = 1000)
# plot(y,col='grey80',pch='.',cex=2)
# lines(mu,col='grey60')
# lines(fit_ebtf$posterior$mean,col=2)
# lines(fit_tf,col=3)
# lines(apply(fit_btf$f,2,mean),col=4)
# legend('topright',c('ebtf','cv.tf','btf','true signal','samples'),lty=c(1,2,3,4,NA),col=c(2,3,4,'grey60','grey80'),pch=c(NA,NA,NA,NA,20))
par(mfrow=c(1,2))
runtime_mean = apply(blocks_res$runtimes,2,mean)[c(1,4,5)]
rmse_mean = apply(blocks_res$rmses,2,mean)[c(1,4,5)]
plot(log2(runtime_mean),rmse_mean,xlab='run time(log2, seconds)',ylab='RMSE',pch=1:3,col=2:4)
legend('bottomright',c('ebtf','cv.tf','btf'),pch=c(1,2,3),col=c(2,3,4))
rr = blocks_res$rmses[,c(1,4,5)]
colnames(rr) = c('ebtf','cv.tf','btf')
boxplot(rr,ylab='RMSE')

runtime_mean = apply(angle_res$runtimes,2,mean)[c(3,4,5)]
rmse_mean = apply(angle_res$rmses,2,mean)[c(3,4,5)]
plot(log2(runtime_mean),rmse_mean,xlab='run time(log2)',ylab='RMSE',pch=1:3,col=2:4)
boxplot(angle_res$rmses[,c(3,4,5)],ylab='RMSE')


###############################
########## run time comparisons

run_time_comp = function(n_list=2^(c(7:17)),snr=3,seed=12345){
  set.seed(seed)
  runtimes = c()

  n = length(mu)
  for (n in n_list) {
    print(paste('Running ',n))
    xx = DJ.EX(n,signal = sqrt(snr))
    mu = xx$blocks
    sigma2 = var(mu)/snr
    y = mu + rnorm(n,0,sqrt(sigma2))
    fit_ebtf0 = ebtf(y,sigma2,ord = 0)
    runtimes = c(runtimes,fit_ebtf0$run_time)
  }
  return(runtimes)
}

runtimes = run_time_comp()
par(mfrow=c(1,1))
plot(2^(c(7:17)),runtimes,xlab='n',ylab='Run time in seconds',pch=20,col=4)

y = MASS::mcycle$accel
x = MASS::mcycle$times

v = smash.gaus(y,v.est = T,joint = T,ebnm_param = list(prior_family='ash'))
fit_ebtf = ebtf(y,v$var.res,sk=c(1e-4,1e-3,1e-2,1e-1,1,2,3,5),m_init=v$mu.res,ord=2)
lines(fit_ebtf$posterior$mean,lwd=2,col=2)


fit_tf = trendfilter(y,ord=2)
fit_tf_cv = cv.trendfilter(fit_tf)
fit_tf = fit_tf$beta[,fit_tf$lambda==fit_tf_cv$lambda.min]

fit_btf = btf(y,k=2)

plot(y,col='grey60',xlab='',ylab='acceleration')
lines(fit_tf,col=1)
lines(fit_ebtf$posterior$mean,lwd=2,col=2)
legend('bottomright',c('ebtf','cv.tf'),lty=c(1,1),col=c(2,1),lwd=c(2,1))


