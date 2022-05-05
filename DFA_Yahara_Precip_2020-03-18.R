# Explore DFA (detrended fluctuation analysis) using Yahara data
# SRC 2019-01-01

rm(list = ls())
graphics.off()

library('nonlinearTseries')
library('pracma')

source('Organize_Precip+Gages_2019-01-01.R')

# INPUT ----------------------------------------------------------------------------------------------------
# read Precip data
load(file='DCRA_precip_1940-2021.Rdata')
dat0 = datall # rename to fit the code below
rm(datall)
print('Madison DCRA precip dataset',quote=F)
print(dat0[1,]) # show first line with column names

windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(dat0$dindex,dat0$PRCP,type='l',col='blue',xlab='year',ylab='precip')

# read P load data
load('AnnLoads_PB+YP_1995-2021.Rdata')
dat1 = Pload  # rename for compatibility with existing code
dat1$TotPload = dat1$PB.kg.d + dat1$YW.kg.d   # PB + YW combined loads
rm(Pload)
# add day index for plotting
dat1$dindex = dat1$year + (dat1$doy/367)
print('Yahara-Windsor P load dataset',quote=F)
print(dat1[1,]) # show first line with column names

windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(dat1$dindex,dat1$TotPload,type='l',col='blue',xlab='year',ylab='P Load')

# Load pigments
#save(BGAdark,Chldark,file='BGA+Chl_dark_centered_v2_2008-2021.Rdata')
# clBGA is centered log BGA, z
load(file='BGA+Chl_dark_centered_Z_2008-2021.Rdata') # z scores remove differences of sensors
print('Mendota Phycocyanin daily dark hours',quote=F)
print(BGAdark[1,])

windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(BGAdark$dindex,BGAdark$zlBGA,type='l',col='forestgreen',xlab='year',ylab='log phyco')

# Select variable & Standardize =====================================================================================================
dat0 = na.omit(dat0)
T = dat0$dindex
X = dat0$PRCP   # precip
#T = dat1$dindex
#X = dat1$TotPload  # P load
#T = BGAdark$dindex
#X = BGAdark$zlBGA  # phycocyanin

# Make Z score
mu = mean(X)
sig2 = var(X)
sig = sqrt(sig2)
Xtr = (X - mu)/sig 
X.acf = acf(Xtr,lag.max=1,type='correlation',plot=F,na.rm=T)
print(c('lag 1 acf = ',X.acf$acf[2]),quote=F)
print('',quote=F)
print('Is Hurst exponent > 0.5 ?',quote=F)
H.estimate = hurstexp(Xtr)
print(H.estimate$Hal)
print('',quote=F)
print('2H = 2 - gamma: Ghil et al 2011',quote=F)
print('gamma = 2 - 2H',quote=F)
print(c('gamma from Hurst ',2-2*H.estimate$Hal),quote=F)
print(c('Days in sample ',length(X)),quote=F)

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
X.dfa = dfa(Xdt,window.size.range = c(10,100),npoints=50,do.plot=F)
dfa.est = estimate(X.dfa,regression.range=NULL,do.plot=F,add.legend=F)

# calculate alpha
xvar = log10(X.dfa$window.sizes)
yvar = log10(X.dfa$fluctuation.function)

windows()
plot(xvar,yvar,type='p',pch=1,col='blue',xlab='window size',ylab='fluctuation function')

print('',quote=F)
xvar2 = xvar^2
mod.dfa = lm(yvar ~ xvar )
print('fit of window size X fluctuation function',quote=F)
print(summary(mod.dfa))
print(c('AIC = ',AIC(mod.dfa)),quote=F)
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




