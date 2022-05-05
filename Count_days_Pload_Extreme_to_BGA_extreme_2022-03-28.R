# Count days from P load peak to next Phycocyanin peak
# SRC 2022-03-28

rm(list = ls())
graphics.off()

# load merged data from Merge_Precip_PLoad_BGA_for_CCF_2022-02-24.R
# save line -->
#save(dat3,file='PPT_Pload_BGAdark_2008-2021.Rdata')
load(file='PPT_Pload_BGAdark_2008-2021.Rdata')
print(dat3[1, ])

# Thresholds
THR.P = 35 # 90th pctile is about 35, used this for Poisson, Pareto threshold is 50
THR.BGA = 1 # 90%ile is 0.855 & used this for Poisson, 95%ile is 1.1277. 
# Used 1 for Pareto. Note unit is s.d.

# select necessary columns
dat4 = subset(dat3,select=c('year','doy','zlBGA','PRCP','TotPload'))

# process by years
uy = unique(dat4$year)
nuy = length(uy)

gaplist = c(365) # dummy value for the first gap
for(i in 1:nuy) {
  ydat = subset(dat4,subset=(year == uy[i]))
  ny = length(ydat$TotPload)
  # find doy with extremes for Pload
  iext = which(ydat$TotPload > THR.P)
  pldoy = dat4$doy[iext]
  # find doy with extremes for BGA
  iext = which(ydat$zlBGA > THR.BGA)
  bgdoy = dat4$doy[iext]
  # find next bga extreme for each Pload extreme
  npl = length(pldoy)
  dlay = rep(0,npl)
  for(ip in 1:npl) {
    delt0 = bgdoy - pldoy[ip]
    delt1 = ifelse(delt0 >= 1,delt0,365) # if delt0 is negative set to high value
    mindelt1 = ifelse(min(delt1)>0,min(delt1),365)
    dlay[ip] = mindelt1
  }
  gaplist = c(gaplist,dlay)
}

# remove the first 'starter' entry
ngaps = length(gaplist)
gaplist = gaplist[2:ngaps]

# remove NA (from when there is no minimum)
gaps = ifelse(is.na(gaplist)==TRUE,365,gaplist)

# summary statistics for the gaps
nPextreme = length(gaps)
# count cases of 'no subsequent bloom that year'
dum = ifelse(gaps == 365,1,0)
nobloom = sum(dum)
print(c('number of extremes ',nPextreme),quote=F)
print(c('extremes with no subsequent bloom that year',nobloom),quote=F)
gaps1 = subset(gaps,subset=(gaps < 365))
print('summary of time interval from Pload extreme to BGA extreme',quote=F)
print(summary(gaps1))

windows()
par(mar=c(4, 4.5, 1.5, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
hist(gaps1,breaks=20,col='palegreen',xlab='Days After Extreme Load',
     ylab='Frequency',main='Extreme P Loads and Blooms')

# lettered version
windows()
par(mar=c(4, 4.5, 1.5, 1) + 0.1, cex.axis=1.8,cex.lab=1.8)
hist(gaps1,breaks=20,col='palegreen',xlab='Days After Extreme Load',
     ylab='Frequency',main='Extreme P Loads and Blooms')
text(x=35,y=30,
     paste(nPextreme,' extreme load events;'),cex=1.8,col='blue')
text(x=35,y=28,
     paste(nobloom,' were not followed by a bloom'),cex=1.8,col='blue')



