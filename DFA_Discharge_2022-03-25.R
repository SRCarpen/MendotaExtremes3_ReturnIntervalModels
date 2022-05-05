# Explore DFA (detrended fluctuation analysis) using Yahara data
# SRC 2019-01-01

rm(list = ls())
graphics.off()

library('nonlinearTseries')
library('pracma')

source('Organize_Precip+Gages_2019-01-01.R')

# INPUT ----------------------------------------------------------------------------------------------------
# read data
load(file='Discharge_PB+YW_1990-2021.Rdata')
# PB+YW discharge
dat1 = na.omit(PBYW2)  # rename for compatibility with existing code
print('Pheasant Branch and Yahara-Windsor P load dataset',quote=F)
print('PB is .x and YW is .y',quote=F)
print(dat1[1,]) # show first line with column names

# Select variable & Standardize =====================================================================================================
T = dat1$dindex.x
X = dat1$Totdis 
X = log(dat1$Totdis) # log transform; optional 

windows()
plot(T,X,type='l',xlab='time index',ylab='discharge')

# Make Z score
mu = mean(X,na.rm=T)
sig2 = var(X,na.rm=T)
sig = sqrt(sig2)
Xtr = (X - mu)/sig 

X.acf = acf(Xtr,lag.max=1,type='correlation',plot=F,na.rm=T)
print(c('lag 1 acf = ',X.acf$acf[2]),quote=F)
print('',quote=F)
print('Is Hurst exponent > 0.5 ?',quote=F)
hurstexp(Xtr)
print('',quote=F)
# Check trends
print('Kendalls tau test for trends',quote=F)
kstat = cor.test(T,Xtr,use='complete.obs',method='kendall')
print(kstat)

# Detrend if needed 
lm0 = lm(Xtr ~ T)
Xdt = lm0$residuals
print('',quote=F)
print('detrending model',quote=F)
print(summary(lm0))

# Detrending slightly improves Hurst exponent
print('',quote=F)
print('Is Hurst exponent > 0.5 ?',quote=F)
H.estimate = hurstexp(Xdt)
print(H.estimate$Hal)
print('',quote=F)
print('2H = 2 - gamma: Ghil et al 2011',quote=F)
print('gamma = 2 - 2H',quote=F)
print(c('gamma from Hurst ',2-2*H.estimate$Hal),quote=F)
print(c('Days in sample ',length(X)),quote=F)

# NOTE:  DETRENDING HAS NO EFFECT ON ALPHA OR GAMMA ESTIMATES FROM DFA

# test DFA
X.dfa = dfa(Xtr,window.size.range = c(10,100),npoints=50,do.plot=F)
dfa.est = estimate(X.dfa,regression.range=NULL,do.plot=F,add.legend=F)

# calculate alpha
xvar = log(X.dfa$window.sizes)
yvar = log(X.dfa$fluctuation.function)

windows()
plot(xvar,yvar,type='p',pch=1,col='blue',xlab='window size',
     ylab='fluctuation function',main='quadratic')

# Quadratic fit to dfa result
print('',quote=F)
xvar2 = xvar^2
mod.dfa = lm(yvar ~ xvar + xvar2)
print('fit of window size X fluctuation function',quote=F)
print(summary(mod.dfa))
print(c('Quadratic AIC = ',AIC(mod.dfa)),quote=F)
# add model estimates to plot
points(xvar,mod.dfa$fitted.values,type='l',col='blue')

# plot alpha versus xvar with curvature
bhat = mod.dfa$coefficients
# linear
#slope.dfa = bhat[2] #+ (2*bhat[3]*xvar)
#gamma.dfa = 2*(1 - slope.dfa)
# quadratic
slope.dfa = bhat[2] + (2*bhat[3]*xvar)
gamma.dfa = 2*(1 - slope.dfa)
windows(width=8,height=5)
par(mfrow=c(1,2),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(xvar,slope.dfa,type='l',col='blue',xlab='window size',ylab='slope, alpha')
plot(xvar,gamma.dfa,type='l',col='blue',xlab='window size',ylab='gamma')

# Calculate gamma
alpha = mod.dfa$coefficients[2]
gamma = 2*(1-alpha)

print('',quote=F)
print('Quadratic')
print(c('autocorrelation lag 1 = ',X.acf$acf[2]),quote=F)
print(c('DFA self-affinity alpha (variance slope) = ',alpha),quote=F)
print(c('DFA Correlation exponent gamma = ',gamma),quote=F)

# Linear fit to dfa result
print('',quote=F)
mod.dfa = lm(yvar ~ xvar )
print('fit of window size X fluctuation function',quote=F)
print(summary(mod.dfa))
print(c('Linear AIC = ',AIC(mod.dfa)),quote=F)
# add model estimates to plot
points(xvar,mod.dfa$fitted.values,type='l',col='blue')

windows()
plot(xvar,yvar,type='p',pch=1,col='blue',xlab='window size',
     ylab='fluctuation function',main='linear')
# add model estimates to plot
points(xvar,mod.dfa$fitted.values,type='l',col='blue')

# Calculate gamma
alpha = mod.dfa$coefficients[2]
gamma = 2*(1-alpha)

# error of gamma: var(gamma) = 4*var(alpha)
p1 = summary.lm(mod.dfa)
p2 = coef(p1)
se.alpha = p2[2,2]
var.alpha = se.alpha^2
var.gamma = 4*var.alpha  # 2* d(mod.dfa)/dxvar
se.gamma = sqrt(var.gamma)

print('',quote=F)
print(c('autocorrelation lag 1 = ',X.acf$acf[2]),quote=F)
print(c('DFA self-affinity alpha (variance slope) = ',alpha),quote=F)
print(c('alpha s.e. = ',se.alpha),quote=F)
print(c('DFA Correlation exponent gamma = ',gamma,' s.e. = ',se.gamma),quote=F)


