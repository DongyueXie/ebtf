library(ebtf)
set.seed(123)
n=2^10
xx = DJ.EX(n,signal = sqrt(3))
mu = xx$heavi

#t = seq(-5,5,length.out = n)
#mu = c(t[t<=0]*-1,t[t>0]*1)
# mu = sin(t*pi)*5

v = rep(1,n)
v[sample(1:n,round(n*0.1))] = 10
v[sample(1:n,round(n*0.1))] = 0.1
x = mu + rnorm(n,0,sqrt(v))

fit = ebtf(x,v,ord=0,printevery = 1)
fit = ebtf(x,v,ord=1,printevery = 1)
fit = ebtf(x,v,ord=2,printevery = 1)

plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean,col=4)

plot(x,col='grey80')
lines(mu,col='grey50')
lines(ti.thresh(x,sqrt(v)))
lines(smash.gaus(x,sqrt(v)))



