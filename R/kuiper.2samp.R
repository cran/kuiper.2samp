#'  2-sample Kuiper Test Function: performs Kuiper Test for two sets samples of observations

#' @param x an array of sample observations
#' @param y the other array of sample observations
#' 
#' @return Kuiper test statistic and p-value
#' 
#' @examples 
#' kuiper.2samp(rnorm(1e3),rnorm(1e3))

#' @note 
#'  The computation of the p-value takes references from Paltani(2004)
#'  which states that the functions (in the set of four formulas) never underestimates the false positive probability
#'  however it can be a bit high when the sample size in the range of 40 to 50
#'  a factor of 1.5 is quoted at the 1e-7 level

#' @references 
#'  Kuiper, N. H. (1960). "Tests concerning random points on a circle". Proceedings of the Koninklijke Nederlandse Akademie van Wetenschappen, Series A. 63: 38-47.
#'  Paltani, S., "Searching for periods in X-ray observations using Kuiper's test. Application to the ROSAT PSPC archive", Astronomy and Astrophysics, v.240, p.789-790, 2004.

#' @export 
kuiper.2samp<-function(x,y)
{
#identify and remove na in x data
    x<-x[!is.na(x)]

#calculate length of x data
  nx<-length(x)
  if(nx<1)
    stop("Not enough 'x' data")
pval<-NULL

#identify and remove na in y data
y<-y[!is.na(y)]

#calculate length of y data
ny<-length(y)
if(ny<1)
stop("Not Enough'y' data")

#obtain the rank of x and y combine series
r<-rank(c(x,y),ties.method="min")

#create CDF of x and y respectively
c<-c(x,y)
csort<-sort(c)
c.min<-min(c(x,y))
c.max<-max(c(x,y))

#The CDF of x and y seperately with first row - the order number, and second row - the frequency (%)
x1<-rbind(sort(x),seq(1/nx,1,by=1/nx))
y1<-rbind(sort(y),seq(1/ny,1,by=1/ny))

#combine CDF_x and CDF_y to the combined scale c(x,y)
temp<-matrix(0,2,nx+ny)

for (i in (1:(nx+ny))){
    temp[1,i]<-ifelse(sum(csort[i]==sort(x))>0,x1[2,which.max(csort[i]==sort(x))],max(temp[1,1:i]))
    temp[2,i]<-ifelse(sum(csort[i]==sort(y))>0,y1[2,which.max(csort[i]==sort(y))],max(temp[2,1:i]))
}

#Kuiper statistics calculates the sum of D+ and D- in absoulute value
kp<-abs(max(temp[1,]-temp[2,]))+abs(max(temp[2,]-temp[1,]))


#Calculation of p-value according to S.Paltani (2004)
n<-nx*ny/(nx+ny)
if (kp<0|kp>2){
  stop("The test statistic much be in the range of (0,2) by definition of the Kuiper Test")
} else if (kp<2/n){ 

#(4) in S.Paltani (2004)
  PVAL<-1-factorial(n)*(kp-1/n)^(n-1)
} else if (kp<3/n){

#(5) in S.Paltani (2004)
  k<--(n*kp-1)/2
  r<-sqrt(k^2-(n*kp-2)/2)
  a<--k+r
  b<--k-r
  PVAL<-1-factorial(n-1)*(b^(n-1)*(1-a)-a^(n-1)*(1-b))/(n^(n-2)*(b-a))
} else if ((kp>0.5&&n%%2==0)|(kp>(n-1)/(2*n)&&n%%2==1)){
  
#(6) in S.Paltani (2004)
  temp<-as.integer(floor(n*(1-kp)))+1
  c<-matrix(0,ncol=1,nrow=temp)
  #Loop for calculating the sum (any other structure better than loop in here?)
  for (t in 0:temp){
  y<-kp+t/n
  tt<-y^(t-3)*(y^3*n-y^2*t*(3-2/n)/n-t*(t-1)*(t-2)/n^2)
  p_temp<-choose(n,t)*(1-kp-t/n)^(n-t-1)*tt
  c[t]<-p_temp
  }
  PVAL<-sum(c)
} else {
  
#(3) in S.Paltani (2004)
  z<-kp*sqrt(n)
  s1<-0
  term<-1e-12
  abs<-1e-100
  for (m in 1:1e+08){
    t1<-2*(4*m^2*z^2-1)*exp(-2*m^2*z^2)
    s<-s1
    s1<-s1+t1
    if ((abs(s1-s)/(abs(s1)+abs(s))<term)|(abs(s1-s)<abs)) break
  }
  s2<-0
  for (m in 1:1e+08){
    t2<-m^2*(4*m^2*z^2-3)*exp(-2*m^2*z^2)
    s<-s2
    s2<-s2+t2
    if ((abs(s2-s)/(abs(s2)+abs(s)))<term|(abs(s1-s)<abs)) break
    PVAL<-s1-8*kp/(3*sqrt(n))*s2
  }
}

#Kuiper statistic and p-value as the outputs
output<-list(Kuiper.statistic=kp, p.value=PVAL)
return(output)
}