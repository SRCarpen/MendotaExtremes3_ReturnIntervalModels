# DFA on merged Precip, Dischg, Pload, BGA 
# Data sets have near-identical N

rm(list = ls())
graphics.off()

library('nonlinearTseries')
library('pracma')

# Functions ===========================================================

# Linear Detrender 
detrend = function(tvec,xvec) {
  lm0 = lm(xvec ~ tvec)
  Xdt = lm0$residuals
  #print('',quote=F)
  #print('detrending model',quote=F)
  #print(summary(lm0))
  return(Xdt)
}

# Detrend by year
dt_by_year = function(Tindex,X) {
  year = trunc(Tindex)
  dtemp = as.data.frame(cbind(year,Tindex,X))
  nt = length(Tindex)
  uy = unique(year)
  nuy = length(uy)
  Xdt = rep(0,nt) # vector to hold detrend by year
  istart = 1 # first index for Xdt
  for(i in 1:nuy) {
    dsub = subset(dtemp,subset=(year == uy[i]))
    Ttemp = dsub$Tindex
    Xtemp = dsub$X
    ntemp = length(Ttemp)
    istop = istart + ntemp - 1
    yeardt = detrend(Ttemp,Xtemp)
    Xdt[istart:istop] = yeardt
    istart = istop+1
  }
  return(Xdt) 
}

# DFA function;  variate name, detrended variate 
dfa.fun = function(vname,Xdt) {  

X.dfa = dfa(Xdt,window.size.range = c(10,100),npoints=50,do.plot=F)
dfa.est = estimate(X.dfa,regression.range=NULL,do.plot=F,add.legend=F)

# calculate alpha
xvar = log10(X.dfa$window.sizes)
yvar = log10(X.dfa$fluctuation.function)

windows()
plot(xvar,yvar,type='p',pch=1,col='blue',xlab='window size',
     ylab='fluctuation function',main=vname)

print('',quote=F)
print(c('DFA for ',vname),quote=F)
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
print(c('DFA self-affinity alpha (variance slope) = ',alpha),quote=F)
print(c('alpha s.e. = ',se.alpha),quote=F)
print(c('DFA Correlation exponent gamma = ',gamma,' s.e. = ',se.gamma),quote=F)

outvec=c(vname,gamma,se.gamma)
return(outvec)
}

# END FUNCTIONS =============================================================

# Merged data from Merge_Precip_Dischg_Pload_BGA)_for_Plots_2022-04-24.R

#save(dat3,file='PPT_Disch_Pload_BGAdark_2008-2021.Rdata')
load(file='PPT_Disch_Pload_BGAdark_2008-2021.Rdata')

print('dimensions after na.omit',quote=F)
print(dim(dat3))

print('all variates',quote=F)
print(dat3[1,])

#-----------------------------------------------------------------------
print('',quote=F)
print('DFA for precipitation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',quote=F)
vname = c('precip')

# Extract variate
Xobs = dat3$PRCP
T = dat3$dindex

# Detrend if needed 
lm0 = lm(Xobs ~ T)
Xdt = lm0$residuals
print('',quote=F)
print('detrending model',quote=F)
print(summary(lm0))

# run DFA
precip.gamma = dfa.fun(vname,Xdt)
print(precip.gamma,quote=F)

#-----------------------------------------------------------------------
print('',quote=F)
print('DFA for discharge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',quote=F)
vname = c('dischg')

# Extract variate
Xobs = dat3$Totdis
T = dat3$dindex

# Detrend if needed 
#lm0 = lm(Xobs ~ T)
#Xdt = lm0$residuals
#print('',quote=F)
#print('detrending model',quote=F)
#print(summary(lm0))

# run DFA
dischg.gamma = dfa.fun(vname,Xobs)
print(dischg.gamma,quote=F)

#-----------------------------------------------------------------------
print('',quote=F)
print('DFA for P load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',quote=F)
vname = c('Pload')

# Extract variate
Xobs = dat3$TotPload
T = dat3$dindex

# Detrend if needed 
#lm0 = lm(Xobs ~ T)
#Xdt = lm0$residuals
#print('',quote=F)
#print('detrending model',quote=F)
#print(summary(lm0))

# run DFA
Pload.gamma = dfa.fun(vname,Xobs)
print(Pload.gamma,quote=F)

#-----------------------------------------------------------------------
print('',quote=F)
print('DFA for BGA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',quote=F)
vname = c('BGA')

# Extract variate
Xobs = dat3$zlBGA
Tindex = dat3$dindex
#T = Tindex


# NOTE:  DETRENDING HAS NO EFFECT ON ALPHA OR GAMMA ESTIMATES FROM DFA
Xdt = dt_by_year(Tindex,Xobs)

# run DFA
BGA.gamma = dfa.fun(vname,Xdt)
print(BGA.gamma,quote=F)

# Input DFA gamma using all data
variate = c('Precip.','Disch.','P Load','Phyco.')
gamma = c(0.939,0.203,0.0769,0.114)
se = c(0.00465,0.0118,0.0130,0.0036)
se2 = 2*se

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.6,cex.main=1.6)
DFAplot = barplot(gamma,names.arg=variate,ylim=c(0,1),col='lightsteelblue',
                  ylab='gamma',xlab='Time Series',
                  main='Correlation Exponent from DFA')

arrows(x0 = DFAplot,                           # Add error bars
       y0 = gamma + se2, #data_summary$mean + data_summary$se,
       y1 = gamma - se2, #data_summary$mean - data_summary$se,
       angle = 90,
       code = 3,
       length = 0.3,lwd=2)

# Gamma vector & S.E. for matching days
gammam = as.numeric(c(precip.gamma[2],dischg.gamma[2],
                      Pload.gamma[2],BGA.gamma[2]))
sem = as.numeric(c(precip.gamma[3],dischg.gamma[3],
                   Pload.gamma[3],BGA.gamma[3]))

sem2 = 2*sem

# Calculate lag versus effect (decay curves tau^(-gamma)  )
tauvec = c(1:30)
cmat = matrix(0,nr=30,nc=4)
for(i in 1:4)  {
  cmat[,i] = tauvec^(-gammam[i])
}

windows(width=6,height=10)
par(mfrow=c(2,1),mar=c(3, 4, 1, 2) + 0.1, 
    cex.axis=1.5,cex.lab=1.6,cex.main=1.4)
cvec = c('skyblue','darkgoldenrod3','magenta','limegreen')
DFAplot = barplot(gammam,names.arg=variate,ylim=c(0,1),col=cvec, #'lightsteelblue',
                  ylab='gamma',xlab=' ',
                  main='A. Correlation Exponent')
arrows(x0 = DFAplot,                           # Add error bars
       y0 = gammam + sem2, #data_summary$mean + data_summary$se,
       y1 = gammam - sem2, #data_summary$mean - data_summary$se,
       angle = 90,
       code = 3,
       length = 0.3,lwd=2)
box(which='plot',lty='solid',lwd=1,col='limegreen')
#
par( mar=c(4, 4, 1, 2) + 0.1 )
plot(tauvec,cmat[,1],type='l',lwd=3,col='skyblue',ylim=c(0,1),xlab='Lag (days)',
     ylab='Weight',main='B. Decay of Autocorrelation')
points(tauvec,cmat[,2],type='l',lwd=3,col='darkgoldenrod3')
points(tauvec,cmat[,3],type='l',lwd=3,col='magenta')
points(tauvec,cmat[,4],type='l',lwd=3,col='limegreen')

