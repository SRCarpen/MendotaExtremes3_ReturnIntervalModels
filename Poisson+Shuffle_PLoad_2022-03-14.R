# Program to count threshold exceedances and measure gaps between them
# SRC 2019-01-09

rm(list = ls())
graphics.off()

library(extRemes)
library(quantreg)
library(reReg)
library(survival)
library(permute)

# Functions

# Determine if a year is a leap year and set the day number for 30 Sept
isleap = function(yr) {
  rem = (2048-yr)%%4
  leapyr = ifelse(rem==0,1,0)
  Sept30 = 273+leapyr
  return(Sept30)
}

# Make a matrix with columns of water year, row DOY, and some variable X between years Y1 and YN
# Water year is year t-1 day 274 (1 Oct) to year t day 273 (30 Sept)
WatYear = function(year,doy,X,Y1,YN) {
  NY = YN-Y1+1 # of years
  # extract first year
  sep30.0 = isleap(Y1-1)
  sep30.1 = isleap(Y1)
  # make x 
  x0 = subset(X,subset=(year==Y1-1 & doy>sep30.0))
  x1 = subset(X,subset=(year==Y1 & doy<=sep30.1))
  x.wyr = c(x0,x1)
  # make water year
  n.wyr = length(x.wyr)
  y.wyr = rep(Y1,n.wyr)
  # make water doy
  d.wyr = (1:n.wyr)
  # save calendar years
  y0 = subset(year,subset=(year==Y1-1 & doy>sep30.0))
  y1 = subset(year,subset=(year==Y1 & doy<=sep30.1))
  y.cal = c(y0,y1)
  # save calendar doys
  d0 = subset(doy,subset=(year==Y1-1 & doy>sep30.0))
  d1 = subset(doy,subset=(year==Y1 & doy<=sep30.1))
  d.cal = c(d0,d1)
  #print(c(sep30.0,sep30.1))
  #print(c(NY,length(x.wyr),length(y.wyr),length(d.wyr),length(y.cal),length(d.cal)))
  # loop over subsequent years
  for(iy in (Y1+1):YN) {
    # extract first year
    sep30.0 = isleap(iy-1)
    sep30.1 = isleap(iy)
    # make x 
    x0 = subset(X,subset=(year==iy-1 & doy>sep30.0))
    x1 = subset(X,subset=(year==iy & doy<=sep30.1))
    x.wyr = c(x.wyr,x0,x1)
    # make water year
    n.wyr = length(x0) + length(x1)
    y.iy = rep(iy,n.wyr)
    y.wyr = c(y.wyr,y.iy)
    # make water doy
    d.iy = (1:n.wyr)
    d.wyr = c(d.wyr,d.iy)
    # save calendar years
    y0 = subset(year,subset=(year==iy-1 & doy>sep30.0))
    y1 = subset(year,subset=(year==iy & doy<=sep30.1))
    y.cal = c(y.cal,y0,y1)
    # save calendar DOYs
    d0 = subset(doy,subset=(year==iy-1 & doy>sep30.0))
    d1 = subset(doy,subset=(year==iy & doy<=sep30.1))
    d.cal = c(d.cal,d0,d1)
    #print(iy)
    #print(c(sep30.0,sep30.1))
    #print(c(length(x.wyr),length(y.wyr),length(d.wyr),length(y.cal),length(d.cal)))
  }
  print('columns: water year, water doy, x, calendar year, calendar doy',quote=F)
  mat.wyr = matrix(c(y.wyr,d.wyr,x.wyr,y.cal,d.cal),nr=length(y.wyr),nc=5,byrow=F)
  return(mat.wyr)
}

# INPUT ----------------------------------------------------------------------------------------------------
# read data
#dat1 = read.csv("Yahara_Windsor_Flow+P.csv", stringsAsFactors=FALSE)
load('AnnLoads_PB+YP_1995-2021.Rdata')
dat1 = Pload  # rename for compatibility with existing code
dat1$TotPload = dat1$PB.kg.d + dat1$YW.kg.d   # PB + YW combined loads
rm(Pload)
print('Yahara-Windsor P load dataset',quote=F)
print(dat1[1,]) # show first line with column names

# Density plot
dkern = density(dat1$TotPload,kernel='epanechnikov',n=128)
windows()
plot(dkern$x,dkern$y,type='l',lwd=2,col='blue')

# Make water year matrices  ------------------------------------------------------------------
print('Water year matrix P load',quote=F)

# make water year matrix for load
Pload = WatYear(dat1$year,dat1$doy,dat1$TotPload,1995,2021)

# test plot
doy.dec = Pload[,1] + (Pload[,2]/366)
windows(width=15,height=6)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doy.dec,Pload[,3],type='l',col='red',ylog='T',xlab='year',ylab='P load')

print('Summarize Pload data',quote=F)
pvec = c(0.5,0.75,0.8,0.9,0.95)
print(summary(Pload[,3]))
print(quantile( Pload[,3],probs = pvec) )
print('Summarize Pload data',quote=F)
print(summary(Pload[,3]))

# FUNCTION for exceedance statistics ----------------------------------
ReturnTimes = function(THR,DOY,X)  {
  # calculate return time statistics for a given
  # threshold (THR), time index (DOY), and variate (X)
  
  # Get indices of observations that exceed the threshold
  ixc = which(X > THR)
  nxc = length(ixc)
  doyxc = DOY[ixc]
  Xixc = X[ixc]
  dixc = diff(ixc)
  
  # summarize gaps
  print('',quote=F)
  print('Summary of gaps between exceedances',quote=F)
  print(summary(dixc))
  print('Variance/Mean',quote=F)
  print(var(dixc)/mean(dixc))
  print('',quote=F)
  
  outlist = list(ixc,doyxc,Xixc,dixc)
  
  print('output list contains',quote=F)
  print('ixc, the indices of X values exceeding the threshold',quote=F)
  print('doyxc, the day numbers of exceedances',quote=F)
  print('Xixc, the vector of X values exceeding the threshold',quote=F)
  print('dixc, the vector of return times',quote=F)
  print('',quote=F)
  
  return(outlist)
  
}   # END FUNCTION -----------------------------------------------------

# Pick an exceedance threshold 
THR = 35 # mm, 25 mm ~ 1 inch; 20 mm was fitted threshold for Pareto fit
abline(h=THR,lty=2,lwd=2)

# Get indices of observations that exceed the threshold
gaplist = ReturnTimes(THR,doy.dec,Pload[,3])
ixc = gaplist[[1]]  # indices of exceedances
nxc = length(ixc)
doyxc = gaplist[[2]]  # day numbers of exceedances
yearxc = trunc(doyxc)  # years in which exceedance occurred
Xixc = gaplist[[3]]   # size of individual exceedances
dixc = gaplist[[4]]   # gap lengths or return intervals

# fit model to gaps
Ngap = length(dixc)
status.gap = rep(1,Ngap)
status.gap[Ngap] = 0
Surv.Mod = Surv(dixc, status.gap)
print('summary of survival model',quote=F)
print(summary(Surv.Mod))
Cox1 = coxph(Surv.Mod ~ doyxc[1:Ngap])
#Cox1 = coxph(dixc ~ doyxc[1:Ngap])
print('summary of zero-order Cox model',quote=F)
print(summary(Cox1))

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 0.5, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doyxc[1:Ngap],dixc,type='p',pch=19,cex=0.7,col='red', #log='y',
     xlab='start day',ylab='gap length')

# calculate stats
N = length(dixc)
mu = mean(dixc)
sig2 = var(dixc)
ratio = sig2/mu

# get range of interval size
rixc = range(ixc)
# scale by mean interval size mu
rR = dixc/mu

# Rank rR
rRank = rank(rR)
murP = 1 - (rRank/N)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(rR,murP,type='p',pch=19,col='blue',cex=0.7,xlog='T',ylog='T',xlab='r/mu',ylab='P(r)')

# Shuffled cumulative ranks -----------------------------------------------------
Xoriginal = Pload[,3]
Nshuf = length(Xoriginal)
ishuf = shuffle(Nshuf)
Xshuffle = Xoriginal[ishuf]

gapshuf = ReturnTimes(THR,doy.dec,Xshuffle)
dixc.shuf = gapshuf[[4]]   # gap lengths or return intervals

# calculate stats
N.shuf = length(dixc.shuf)
mu.shuf = mean(dixc.shuf)

# get range of interval size
rixc = range(dixc.shuf)
# scale by mean interval size mu
rR.shuf = dixc.shuf/mu.shuf

# Rank rR
rRank.shuf = rank(rR.shuf)
murP.shuf = 1 - (rRank.shuf/N.shuf)
# add to previous plot
points(rR.shuf,murP.shuf,type='p',pch=17,col='magenta',cex=0.7)

# Poisson cumulative ranks ---------------------------------------------------
xrange = range(rR)
xvals = seq(xrange[1],xrange[2],length.out=30)
yvals = (exp(-xvals)) 

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(xvals,yvals,type='l',lwd=2,col='red',xlog='T',ylog='T',xlab='r/mu',ylab='P(r)')
# add shuffle to previous plot
points(rR.shuf,murP.shuf,type='p',pch=20,col='mediumorchid',cex=1)

# Overplot
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(xvals,yvals,type='l',lwd=3,col='red',ylog='T',xlab='gap/mean',ylab='mean*p(gap)',
     main='')
points(rR,murP,type='p',pch=19,col='blue',cex=0.7,ylog='T')#,xlab='gap/mean',ylab='P(gap)')
grid()
legend('topright',legend=c('Expected if gaps are independent','Observed'),
       lwd=c(3,NA),pch=c(NA,19),col=c('red','blue'),cex=1.5,bty='n')

windows(width=12,height=8)
par(mfrow=c(2,1),mar=c(4, 4.2, 0.5, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doy.dec,Pload[,3],type='l',col='forestgreen',ylog='T',xlab='year',ylab='Precip, mm/d')
abline(h=THR,lty=2,lwd=2)
plot(xvals,yvals,type='l',lwd=2,col='red',
     log='y',xlab='Interval Between Events / Mean Interval',ylab='Log10 Rank / N')
points(rR,murP,type='p',pch=19,col='forestgreen',cex=0.7,xlab='gap/mean',ylab='P(gap)')

# Calculate dispersion index by year
start.doy = doyxc[1:N]
start.y = yearxc[1:N]
uy = unique(start.y)
nuy = length(uy)
#
disp.mat = matrix(0,nr=nuy,nc=5) 
#
for(i in 1:nuy) {
  disp.mat[i,1] = uy[i]  # year
  dat.xc = subset(dixc,subset=(yearxc==uy[i]))
  disp.mat[i,2] = length(dat.xc)  #N
  disp.mat[i,3] = mean(dat.xc,na.rm=T)    # mean
  disp.mat[i,4] = var(dat.xc,na.rm=T)     # var
  disp.mat[i,5] = ifelse(disp.mat[i,2]==1,0,disp.mat[i,4]/disp.mat[i,3])
}

#print(disp.mat)

print('',quote=F)
print('Kendall correlation year x N/year',quote=F)
print(cor.test(disp.mat[,1],disp.mat[,2],alternative='two.sided',method='kendall'))
print('Kendall correlation year x dispersion',quote=F)
print(cor.test(disp.mat[,1],disp.mat[,5],alternative='two.sided',method='kendall'))

# 4 panel figure
windows()
par(mfrow=c(4,1),mar=c(4, 4.5, 1.5, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
plot(doy.dec,Pload[,3],type='l',col='forestgreen',ylog='T',xlab='',ylab='Precip., mm/d')
abline(h=THR,lty=2,lwd=2)
#
plot(disp.mat[,1],disp.mat[,2],type='b',pch=19,lty=2,cex=1.5,lwd=2,col='blue',xlab='',
     ylab='Extremes/Year')
grid()
plot(disp.mat[,1],disp.mat[,5],type='b',log='y',pch=19,lty=2,cex=1.5,lwd=2,col='blue',
     xlab='Year',ylab='Dispersion')
grid()
plot(xvals,yvals,type='l',lwd=2,col='red',
     log='y',xlab='Interval Between Events / Mean Interval',ylab='Rank / N')
points(rR,murP,type='p',pch=19,col='blue',cex=0.9,xlab='gap/mean')
grid()

# Try Layout
m.plot = rbind(c(1,4),c(2,4),c(3,4))
print(m.plot)
windows(width=12,height=8)
layout(m.plot)
par(mar=c(2, 4.5, 1, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
plot(doy.dec,Pload[,3],type='l',col='forestgreen',ylog='T',xlab='',ylab='Precip., mm/d')
abline(h=THR,lty=2,lwd=2)
#
plot(disp.mat[,1],disp.mat[,2],type='b',pch=19,lty=2,cex=1.5,lwd=2,col='blue',xlab='',
     ylab='Extremes/Year')
grid()
par(mar=c(4, 4.5, 1, 1) + 0.1)
plot(disp.mat[,1],disp.mat[,5],type='b',log='y',pch=19,lty=2,cex=1.5,lwd=2,col='blue',
     xlab='Year',ylab='Dispersion')
grid()
plot(xvals,yvals,type='l',lwd=2,col='red',ylim=c(0.001,1),xlim=c(0,8),
     log='y',xlab='Interval Between Events / Mean Interval',ylab='Rank / N')
points(rR,murP,type='p',pch=19,col='blue',cex=1.1,xlab='gap/mean')
grid()

# Layout omitting daily precip
m.plot = rbind(c(1,3),c(2,3))
print(m.plot)
windows(width=12,height=8)
layout(m.plot)
par(mar=c(4, 4.5, 1.5, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
plot(disp.mat[,1],disp.mat[,2],type='b',pch=19,lty=2,cex=1.5,lwd=2,col='blue',xlab='',
     ylab='Extremes/Year')
grid()
text(1996,72,'A',cex=1.8)
plot(disp.mat[,1],disp.mat[,5],type='b',log='y',pch=19,lty=2,cex=1.5,lwd=2,col='blue',
     xlab='Year',ylab='Dispersion')
grid()
text(2020,100,'B',cex=1.8)
plot(xvals,yvals,type='l',lwd=2,col='red',ylim=c(0.002,1),xlim=c(0,9),
     log='y',xlab='Interval Between Events / Mean Interval',ylab='Rank / N')
points(rR,murP,type='p',pch=19,col='blue',cex=1.1,xlab='gap/mean')
# add shuffle to previous plot
points(rR.shuf,murP.shuf,type='p',pch=17,col='mediumorchid',cex=1.1)
points(xvals,yvals,type='l',lwd=2,col='red')
grid()
legend('bottomleft',legend=c('Poisson','Data','Shuffled Data'),
       lwd=c(3,NA,NA),pch=c(NA,19,17),col=c('red','blue','mediumorchid'),cex=1.5,bty='n')
text(7.6,0.8,'C',cex=1.8)

# Same as above with no log scale of Poisson plot
m.plot = rbind(c(1,3),c(2,3))
print(m.plot)
windows(width=12,height=8)
layout(m.plot)
par(mar=c(4, 4.5, 1.5, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
plot(disp.mat[,1],disp.mat[,2],type='b',pch=19,lty=2,cex=1.5,lwd=2,col='blue',xlab='',
     ylab='Extremes/Year')
grid()
text(1996,72,'A',cex=1.8)
plot(disp.mat[,1],disp.mat[,5],type='b',log='y',pch=19,lty=2,cex=1.5,lwd=2,col='blue',
     xlab='Year',ylab='Dispersion')
grid()
text(2020,100,'B',cex=1.8)
plot(xvals,yvals,type='l',lwd=2,col='black',ylim=c(0.001,1),xlim=c(0,9),#log='y',
     xlab='Interval Between Events / Mean Interval',ylab='Rank / N')
points(rR,murP,type='p',pch=19,col='blue',cex=1.1,xlab='gap/mean')
# add shuffle to previous plot
points(rR.shuf,murP.shuf,type='p',pch=17,col='red',cex=1.1)
points(xvals,yvals,type='l',lwd=2,col='black')
grid()
legend(x=4,y=0.8,legend=c('Poisson','Data','Shuffled Data'),
       lwd=c(3,NA,NA),pch=c(NA,19,17),col=c('black','blue','red'),cex=1.5,bty='n')
text(8,1,'C',cex=1.8)

# check correlations
cmat = as.data.frame(cbind(yvals,murP,murP.shuf))
print('correlations',quote=F)
print(cor(cmat,use='pairwise.complete.obs'))
print(' ',quote=F)
print('2 sided K-S test, Poisson vs observed',quote=F)
ks1 = ks.test(yvals,murP,alternative='two.sided')
print(ks1)

print('2 sided K-S test, Poisson vs shuffled',quote=F)
ks1 = ks.test(yvals,murP.shuf,alternative='two.sided')
print(ks1)

print('2 sided K-S test, data vs shuffled',quote=F)
ks1 = ks.test(murP,murP.shuf,alternative='two.sided')
print(ks1)

print('Euclidean distances',quote=F)
dtest = dist(x = t(cmat), method = "euclidean", diag = F, upper = F, 
             p = 2)
print(dtest)
