# geometricMVE
## R package for fitting parametric and piecewise-linear models for the gauge function of a multivariate sample in exponential margins

## References
Wadsworth, J. L. and Campbell, R. (2024) Statistical inference for multivariate extremes via a geometric approach JRSSB, 86 (5), 1243-1265
Campbell, R. and Wadsworth, J. L. (2024) Piecewise-linear modeling of multivariate geometric extremes. https://arxiv.org/abs/2412.05195

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

r<-x[,1]+x[,2]
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
qqfiag(fit1)

```
