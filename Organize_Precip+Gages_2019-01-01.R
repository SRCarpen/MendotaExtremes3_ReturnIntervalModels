# Functions to organize and conform Yahara data from DCRA and streams
# SRC 2019-01-01

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
