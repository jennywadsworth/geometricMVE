# geometricMVE
## R package for fitting parametric and piecewise-linear models for limit set based multivariate extremal inference in exponential margins
 
## References
- Wadsworth, J. L. and Campbell, R. (2024) Statistical inference for multivariate extremes via a geometric approach JRSSB, 86 (5), 1243-1265 https://academic.oup.com/jrsssb/article/86/5/1243/7637082
- Campbell, R. and Wadsworth, J. L. (2024) Piecewise-linear modeling of multivariate geometric extremes. https://arxiv.org/abs/2412.05195
- Lee, J. and Wadsworth, J. L. (2025) Geometric criteria for identifying extremal dependence and flexible modeling via additive mixtures. https://arxiv.org/abs/2512.24392

## Example code
### Installation
```r
library(remotes)
remotes::install_github("jennywadsworth/geometricMVE")
```
### 2d Example
```r
library(evd)
set.seed(1)
n<-2500
x<-rbvevd(n, dep=0.5, mar1=c(0,1,0))
x<-qexp(exp(-exp(-x)))

# Radial-angular variables

r<-apply(x,1,sum)
w<-x/r

# Establishing a high threshold for R|W

# Empirical:
qr<-fit.thresh(r=r,w=w, method = "empirical")

# KDE:
qr<-fit.thresh(r=r,w=w, method = "KDE")

# Visualization
plotfittedthresh(qr)

# Fitting truncated gamma model, using R exceedances defined through output of fit.thresh

# Fit using different basic parametric gauge functions

fit1<-fit.geometric.par(thresh.fit=qr,model="log") # logistic-type gauge
fit2<-fit.geometric.par(thresh.fit=qr,model="gauss") # gaussian-type gauge
fit3<-fit.geometric.par(thresh.fit=qr,model="invlog") # Inverted logistic-type gauge

# Fit using additive mixtures (agnostic to extremal dependence class)
fit4<-fit.geometric.par(thresh.fit=qr,model="expgauss") # Exponential-Gaussian mix
fit5<-fit.geometric.par(thresh.fit=qr,model="expinvlog") # Exponential-Inverted logistic mix

# Plot fitted gauge functions / limit sets
plotfittedgauge(fit1)
plotfittedgauge(fit4)


# Diagnostics
ppdiag(fit1)
ppdiag(fit1, type="horizontal") # clearer scale
qqdiag(fit1) # PP plot trasnformed on to standard exponential margins

# Fit using piecewise-linear gauge function

fitpwl1<-fit.geometric.pwl(thresh.fit=qr) 
plotfittedgauge(fitpwl1)

# Fix shape

fitpwl2<-fit.geometric.pwl(thresh.fit=qr,fixshape = T) 
plotfittedgauge(fitpwl2)

# Bound fit

fitpwl3<-fit.geometric.pwl(thresh.fit=qr,bound.fit = T) 
plotfittedgauge(fitpwl3)

# Change reference locations

fitpwl4<-fit.geometric.pwl(thresh.fit=qr, locs=seq(0,1,len=5)) 
plotfittedgauge(fitpwl4)


# Diagnostics
ppdiag(fitpwl1)
qqdiag(fitpwl1)


# Simulation of new pseudo-observations using the model structure
newx<-sim.geometric(fit = fit1, nsim=10000)
plot(x,xlim=c(0,12),ylim=c(0,12))
points(newx,col=2,pch=20)

# Simulate above a higher threshold
newx2<-sim.geometric(fit = fit1, k=2, nsim=10000)
plot(x,xlim=c(0,15),ylim=c(0,15))
points(newx2,col=3,pch=20)


# Estimate probability 8<X<12, 4<Y<6

cond.prob<-mean(newx[,1]>8 & newx[,1]<12 & newx[,2]>4 & newx[,2]<6)
prob<-cond.prob*mean(excind)
prob

cond.prob1<-mean(newx2[,1]>8 & newx2[,1]<12 & newx2[,2]>4 & newx2[,2]<6)
cond.prob2<-Rexc.prob.k(k=2, fit=fit1) # Undoes extra conditioning of being above a higher threshold
prob<-cond.prob1*cond.prob2*mean(excind)
prob

```
#### Note on providing threshold exceedances directly
For simplicity, fit.geometric.par and fit.geometric.pwl take the output of fit.thresh and calculate threshold exceedances from that. It is also possible to supply threshold exceedances and values directly, which will be needed if they are not calculated via fit.thresh. Example syntax for this (in this case using fit.thresh for the threshold calculation):

```
excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind,]
r0wexc<-qr$r0w[excind]

fit1<-fit.geometric.par(r=rexc,w=wexc,r0w=r0wexc,model="log") # Same as fit1 above
```

### 3d Example

```
library(mvtnorm)
library(rgl) # for plots

set.seed(2)
n<-2500
sig<-matrix(c(1,0.8,0.5,0.8,1,0.6,0.5,0.6,1),nrow=3,byrow=T)
x<-qexp(pnorm(rmvnorm(n,sigma=sig)))

# Radial-angular variables

r<-apply(x,1,sum)
w<-x/r

# Establishing a high threshold for R|W

# Empirical:
qr<-fit.thresh(r=r,w=w, method = "empirical")

# KDE:
qr<-fit.thresh(r=r,w=w, method = "KDE")

# Visualization
plotfittedthresh(qr)

# Fitting truncated gamma model, using different basic parametric gauge functions

fit1<-fit.geometric.par(thresh.fit=qr,model="log") # logistic-type gauge
fit2<-fit.geometric.par(thresh.fit=qr,model="gauss") # gaussian-type gauge
fit3<-fit.geometric.par(thresh.fit=qr,model="invlog") # Inverted logistic-type gauge

# Plot fitted gauge functions / limit sets
plotfittedgauge(fit2)

# Add scaled points
points3d(x/log(n))

# Diagnostics
ppdiag(fit2)
ppdiag(fit2, type="horizontal")
qqdiag(fit2)

# Fit using piecewise-linear gauge function

# Define reference angles

par.locs = seq(0,1,by=1/6)
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])


fitpwl1<-fit.geometric.pwl(thresh.fit=qr,locs=par.locs) 

# Plot fitted gauge function
plotfittedgauge(fitpwl1)

# Add scaled points
points3d(x/log(n))

# Diagnostics
ppdiag(fitpwl1)
qqdiag(fitpwl2)


# Simulation of new pseudo-observations using the model structure
newx<-sim.geometric(fit = fit2, nsim=10000)
plot3d(x,xlim=c(0,12),ylim=c(0,12),zlim=c(0,12))
points3d(newx,col=2,pch=20)
```
