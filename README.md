# geometricMVE
## R package for fitting parametric and piecewise-linear models for the gauge function of a multivariate sample in exponential margins

## References
- Wadsworth, J. L. and Campbell, R. (2024) Statistical inference for multivariate extremes via a geometric approach JRSSB, 86 (5), 1243-1265
- Campbell, R. and Wadsworth, J. L. (2024) Piecewise-linear modeling of multivariate geometric extremes. https://arxiv.org/abs/2412.05195

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

# Fitting truncated gamma model, using R exceedances defined through rolling-windows quantiles

excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind,]

# Threshold value corresponding to each w:

r0w<-qr$r0w[excind]

# Fit using different basic parametric gauge functions

fit1<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="log") # logistic-type gauge
fit2<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="gauss") # gaussian-type gauge
fit3<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="invlog") # Inverted logistic-type gauge

# Plot fitted gauge functions / limit sets
plotfittedgauge(fit1)

# Diagnostics
ppdiag(fit1)
qqdiag(fit1)

# Fit using piecewise-linear gauge function

fitpwl1<-fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w) 
plotfittedgauge(fitpwl1)

# Fix shape

fitpwl2<-fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,fixshape = T) 
plotfittedgauge(fitpwl2)

# Bound fit

fitpwl3<-fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,bound.fit = T) 
plotfittedgauge(fitpwl3)

# Change reference locations

fitpwl4<-fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w, locs=seq(0,1,len=5)) 
plotfittedgauge(fitpwl4)


# Diagnostics
ppdiag(fitpwl1)
qqdiag(fitpwl1)
```
### 3d Example

```
library(mvtnorm)
library(rgl) # for plots

set.seed(1)
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

# Fitting truncated gamma model, using R exceedances defined through rolling-windows quantiles

excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind,]

# Threshold value corresponding to each w:

r0w<-qr$r0w[excind]

# Fit using different basic parametric gauge functions

fit1<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="log") # logistic-type gauge
fit2<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="gauss") # gaussian-type gauge
fit3<-fit.geometric.par(r=rexc,w=wexc,r0w=r0w,model="invlog") # Inverted logistic-type gauge

# Plot fitted gauge functions / limit sets
plotfittedgauge(fit2)

# Add scaled points
points3d(x/log(n))

# Diagnostics
ppdiag(fit2)
qqdiag(fit2)

# Fit using piecewise-linear gauge function

# Define reference angles

par.locs = seq(0,1,by=1/6)
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])


fitpwl1<-fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,locs=par.locs) 

# Plot fitted gauge function
plotfittedgauge(fitpwl1)

# Add scaled points
points3d(x/log(n))

# Diagnostics
ppdiag(fitpwl1)
qqdiag(fitpwl2)
```
