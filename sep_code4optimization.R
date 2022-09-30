

library(nloptr)
################################################################################
#####One year optimiztion#######################################################
################################################################################


p=4
n=10000
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))#theta for each county each quarter
ii=2
nn=3
cii=array(NA, dim = c(100,n))
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  
  
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  
  
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  
  mm<-seq(101,501,by=100)#
  
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }

  
  postnext<-cbind( postXmm[,ii:nn],rep(7 ,100))
  cii[,ll]<-exp(rep(postbeta0,100)+postnext%*%postbeta[2:p]+rep(s2,1))*thetasave[,6,ll]
  
  
}

ci<-rowMeans(cii)
postbeta<-colMeans(sampleused[[1]][,which(names(sampleused[[1]][1,])=='beta[1]'):(which(names(sampleused[[1]][1,])=='beta[1]')+3)])


lb <- rep(((0)-mu[1])/sd[1],100)
ub <- rep((max(Bupdiff$Bupdiff[101:600])-(mean(Bupdiff$Bupdiff[101:600])))/
            sd(Bupdiff$Bupdiff[101:600]),100)
postbeta<-colMeans(sampleused[[1]][,which(names(sampleused[[1]][1,])=='beta[1]'):(which(names(sampleused[[1]][1,])=='beta[1]')+3)])
beta1hat<-postbeta[1]

x0 <- c(rep(0,100))
x00 <-rep(((0)-mu[1])/sd[1],100)
fn1 <- function(x) #one year num
{
  obj = 0
  for (i in 1:100) {
    a=ci[i]*exp(beta1hat*x[i])
    obj = obj + a 
    #ci includes other terms
  }
  
  return(obj)
}
fn1(rep((0-mu)/sd,100))

fn2 <- function(x)  #one year rate
{
  obj = 0
  for (i in 1:100) {
    a=(ci[i]*exp(beta1hat*x[i]))/Exp[i,6]
    obj = obj + a 
  }
  
  return(obj/100)
}
#fn2(rep((0-mu)/sd,100))

eqn2 <-function(x)#6743
{
  sumt = 0
  for (i in 1:100) {
    c=x[i]
    sumt = sumt + c
  }
  return(sumt-37.18)#mean(c(5669,8566,5996)), (6743.667-100*mu)/sd
}

resnum_6743<-auglag(x0,
                    fn1,
                    heq = eqn2,
                    lower          = lb,
                    upper          = ub,
                    localsolver = "mma",
                    localtol=1e-05)

BUPincr<-resnum_6743$par
result4num<-data.frame(cbind(subregion=final.data$Place[1:100]),num=round(BUPincr*sd[1]+mu[1]))
result4num$num <- as.numeric(as.character(result4num$num))
result4num[with(result4num, order(-result4num[,2])),]
result4plot_num<-result4num[with(result4num, order(-result4num[,2])),]



resrate_6743<-auglag(x0,
                    fn2,
                    heq = eqn2,
                    lower          = lb,
                    upper          = ub,
                    localsolver = "mma",
                    localtol=1e-05)

BUPincr<-resrate_6743$par
result4rate<-data.frame(cbind(subregion=final.data$Place[1:100]),num=round(BUPincr*sd[1]+mu[1]))
result4rate$num <- as.numeric(as.character(result4rate$num))
result4rate[with(result4rate, order(-result4rate[,2])),]
result4plot_rate<-result4rate[with(result4rate, order(-result4rate[,2])),]


rate<-summary(result4plot_rate$num)
num<-summary(result4plot_num$num)

rbind(rate,num)
dim(result4num[result4num$num==0,])
#77  2
dim(result4num[result4rate$num==0,])
#19  2

library(xtable)

rateTABLE<- data.frame(cbind(result4plot_rate[1:50,],result4plot_rate[51:100,]))
colnames(rateTABLE) <- c('subregion','value','subregion','value')
xtable(rateTABLE, type = "latex", file = "result4plot.tex")


numTABLE<- data.frame(cbind(result4plot_num[1:50,],result4plot_num[51:100,]))
colnames(numTABLE) <- c('subregion','value','subregion','value')
xtable(numTABLE, type = "latex", file = "result4plot.tex")

rateTABLE<- data.frame(rbind(rate,num))
colnames(rateTABLE) <- c("Min.", "1st Qu." ,"Median" , "Mean" ,"3rd Qu.", "Max.")
xtable(rateTABLE, type = "latex", file = "result4plot.tex")

#line plot
p=4
n=10000
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))
#theta for each county each quarter,
#thetasave[,7:9,] are baseline 10:12 are prediction with optimal solution

ii=2
nn=3

# one year num 
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  
  
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  
  
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  
  mm<-seq(101,501,by=100)#101-3701 q2-q38
  
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  
  
  postX3900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX4000<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39<-rep(postbeta0,100)+postX3900[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40<-rep(postbeta0,100)+postX4000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41<-rep(postbeta0,100)+postX4100[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta39)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta40)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta41)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  
  postX3900x<-cbind(resnum_6743$par, postXmm[,ii:nn],rep(7 ,100))
  postX4000x<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100x<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39x<-rep(postbeta0,100)+postX3900x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40x<-rep(postbeta0,100)+postX4000x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41x<-rep(postbeta0,100)+postX4100x[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta39x)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta40x)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta41x)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  
  
  
}



#probability for 2022
statesum<-apply(thetasave[,10,],2, sum)
statebase<-apply(thetasave[,7,],2,sum)
length(which((statebase-statesum)>0))
#probability for 2022, 23, 24
statesum<-apply(thetasave[,10,],2, sum)+apply(thetasave[,11,],2, sum)+apply(thetasave[,12,],2, sum)
statebase<-apply(thetasave[,7,],2,sum)+apply(thetasave[,8,],2,sum)+apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))
#probability for 2024
statesum<-apply(thetasave[,12,],2, sum)
statebase<-apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))


statesum<-apply(thetasave[,12,],2, sum)
statebase<-0.8*apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))


statesum<-apply(thetasave[,12,],2, sum)
statebase<-0.8*apply(thetasave[,6,],2,sum)
length(which((statebase-statesum)>0))



# one year rate 
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  
  mm<-seq(101,501,by=100)#101-3701 q2-q38
  
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  
  
  postX3900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX4000<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39<-rep(postbeta0,100)+postX3900[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40<-rep(postbeta0,100)+postX4000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41<-rep(postbeta0,100)+postX4100[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta39)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta40)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta41)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  
  postX3900x<-cbind(resrate_6743$par, postXmm[,ii:nn],rep(7 ,100))
  postX4000x<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100x<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39x<-rep(postbeta0,100)+postX3900x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40x<-rep(postbeta0,100)+postX4000x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41x<-rep(postbeta0,100)+postX4100x[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta39x)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta40x)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta41x)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  
  
  
}

#probability for 2022
statesumR<-apply(apply(thetasave[,10,], 2, '/',Exp[,6]),2, sum)
statebaseR<-apply(apply(thetasave[,7,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))

#probability for 2022, 23, 24
statesum<-apply(apply(thetasave[,10,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,11,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,12,], 2, '/',Exp[,6]),2, sum)

statebase<-apply(apply(thetasave[,7,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,8,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))



low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)

#statesum<-apply(thetasave[,10:12,], 3, sum)
#hist(statesum)
#abline(v=10182.6,col="red")
#length(statesum[which(statesum<10182.6)])/10000


thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(thetasave[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(thetasave[i,j,])[1]
    high[i,j]= hdi(thetasave[i,j,])[2]
    mean[i,j]=mean(thetasave[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(thetasave[i,(9+j),])[1]
    prelow[i,j]= hdi(thetasave[i,(9+j),])[2]
    Predicted[i,j]=mean(thetasave[i,(9+j),],na.rm = T)
  }
}




y=matrix(as.numeric(final.data$death),nrow=100,byrow=F)




estima2<-data.frame(time=c(2016:2024),Estimated=colSums(mean), 'Upper Bound'=colSums(high),
                    'Lower Bound'=colSums(low), Observed=c(colSums(y[,1:6],na.rm = T),rep(NA,3)),
                    Predicted=c(rep(NA,6),colSums(Predicted)),
                    prehigh=c(rep(NA,6),colSums(prehigh)), 
                    prelow=c(rep(NA,6),colSums(prelow)))



newexp<-cbind(Exp,Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6])
# use one year rate, then repeat line 379-413
estima_r<-data.frame(time=c(2016:2024),Estimated=colMeans(mean/newexp[,1:9]), 'Upper Bound'=colMeans(high/newexp[,1:9]),
                     'Lower Bound'=colMeans(low/newexp[,1:9]), Observed=c(colMeans(y[,1:6]/Exp,na.rm = T),rep(NA,3)),
                     Predicted=c(rep(NA,6),colMeans(Predicted/newexp[,10:12])),
                     prehigh=c(rep(NA,6),colMeans(prehigh/newexp[,10:12])), 
                     prelow=c(rep(NA,6),colMeans(prelow/newexp[,10:12])))



COL=c("Estimated" = "black",
      'Credible Bound' = "steelblue",
      'Credible Bound' = "steelblue",
      "Observed" = "red",
      "Predicted" = "green",
      "Predicted Bound" = "darkgreen",
      "Predicted Bound" = "darkgreen" 
)
oneyearnum<-ggplot(data=estima2,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Number")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())


COL=c("Estimated" = "black",
      'Credible Bound' = "steelblue",
      'Credible Bound' = "steelblue",
      "Observed" = "red",
      "Predicted" = "green",
      "Predicted Bound" = "darkgreen",
      "Predicted Bound" = "darkgreen" 
)

oneyearrate<-ggplot(data=estima_r,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Rate")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())




#map 
usa_counties = map_data("county")
nc = subset(usa_counties, region == "north carolina") 

library(usa)
zcs <- usa::zipcodes
head(zcs)
NCcity <- zcs[zcs[ ,3] == "NC", ]

Charlotte<- c(city= "Charlotte", lat=mean(NCcity[NCcity$city== "Charlotte",]$lat),long=mean(NCcity[NCcity$city== "Charlotte",]$long))
Raleigh<- c(city= "Raleigh", lat=mean(NCcity[NCcity$city== "Raleigh",]$lat),long=mean(NCcity[NCcity$city== "Raleigh",]$long))
Greensboro<- c(city= "GS", lat=mean(NCcity[NCcity$city== "Greensboro",]$lat),long=max(NCcity[NCcity$city== "Greensboro",]$long))
Durham<- c(city= "Durham", lat=mean(NCcity[NCcity$city== "Durham",]$lat),long=mean(NCcity[NCcity$city== "Durham",]$long))
WinstonSalem<- c(city= "WS", lat=max(NCcity[NCcity$city== "Winston Salem",]$lat),long=mean(NCcity[NCcity$city== "Winston Salem",]$long))
Fayetteville<- c(city= "Fayetteville", lat=mean(NCcity[NCcity$city== "Fayetteville",]$lat),long=mean(NCcity[NCcity$city== "Fayetteville",]$long))
Cary<- c(city= "Cary", lat=mean(NCcity[NCcity$city== "Cary",]$lat),long=mean(NCcity[NCcity$city== "Cary",]$long))
Wilmington<- c(city= "Wilmington", lat=mean(NCcity[NCcity$city== "Wilmington",]$lat),long=mean(NCcity[NCcity$city== "Wilmington",]$long))
HighPoint<- c(city= "High Point", lat=mean(NCcity[NCcity$city== "High Point",]$lat),long=mean(NCcity[NCcity$city== "High Point",]$long))
Concord<- c(city= "Concord", lat=mean(NCcity[NCcity$city== "Concord",]$lat),long=mean(NCcity[NCcity$city== "Concord",]$long))
Asheville<- c(city= "Asheville", lat=mean(NCcity[NCcity$city== "Asheville",]$lat),long=mean(NCcity[NCcity$city== "Asheville",]$long))
Monroe<-c(city= "Monroe", lat=mean(NCcity[NCcity$city== "Monroe",]$lat),long=mean(NCcity[NCcity$city== "Monroe",]$long))


cityct<-as.data.frame(rbind(Charlotte, Raleigh,Greensboro,Durham,WinstonSalem,Fayetteville,Asheville))

cityct$lat=as.numeric(cityct$lat)
cityct$long=as.numeric(cityct$long)
cityct$cityPopulation=c(9.25290, 4.88334,3.03286 ,2.92301,2.52175, 2.14384,
                        0.92)
cityct$County=c("Mecklenburg", "Wake","Guilford" ,"Durham","Forsyth", "Cumberland"
                ,"Buncombe"	)

cityct$Population=c(1110356, 1111761,537174 ,321488,382295,335509,
                    261191)
cityct$Population=cityct$Population/1000000


BUPincr<-resrate_6743$par
BUPincr<-resnum_6743$par
result4plotactnum<-data.frame(cbind(num=BUPincr*sd[1]+mu[1],subregion=final.data$Place[1:100]))
BUPincr1<-as.numeric(result4plotactnum$num)

result4plot<-data.frame(cbind(BUPincr1,subregion=final.data$Place[1:100]))
result4plot$subregion<-tolower(result4plot$subregion)
map.data = merge(nc,result4plot,by.x="subregion",by.y="subregion")
map.data$BUPincr1<-as.numeric(map.data$BUPincr1)
ggplot()+
  geom_polygon(data=map.data,aes(x=long,y=lat,group=subregion,fill=BUPincr1),color='black',alpha=.8,size=.3)+
  scale_fill_gradient2(name="Increase",low='white',high='darkblue')+
  geom_point(data=cityct,aes(x=long, y=lat,size=Population),shape = 17,color="pink") +
  coord_map()+
  theme_void()+
  theme(legend.text=element_text(size=12))+
  theme(legend.title =element_text(size=14))






################################################################################
#####Three year optimiztion#####################################################
################################################################################


p=4
n=10000
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))#theta for each county each quarter
ii=2
nn=3
cii=array(NA, dim = c(100,n))
dii=array(NA, dim = c(100,n))
eii=array(NA, dim = c(100,n))
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  mm<-seq(101,501,by=100)#
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  #other terms than bupe
  postnext<-cbind( postXmm[,ii:nn],rep(7 ,100))
  cii[,ll]<-exp(rep(postbeta0,100)+postnext%*%postbeta[2:p]+rep(s2,1))*thetasave[,6,ll]
  postnext1<-cbind( postXmm[,ii:nn],rep(8 ,100))
  dii[,ll]<-exp(rep(postbeta0,100)+postnext1%*%postbeta[2:p]+rep(s2,1))
  postnext2<-cbind( postXmm[,ii:nn],rep(9 ,100))
  eii[,ll]<-exp(rep(postbeta0,100)+postnext2%*%postbeta[2:p]+rep(s2,1))

  
  
}

ci<-rowMeans(cii)
di0<-rowMeans(dii)
ei0<-rowMeans(eii)

postbeta<-colMeans(sampleused[[1]][,which(names(sampleused[[1]][1,])=='beta[1]'):(which(names(sampleused[[1]][1,])=='beta[1]')+3)])
beta1hat<-postbeta[1]

fn4 <- function(x) 
{
  
  
  obj = 0
  for (i in 1:100) {
      a= ci[i]*exp(beta1hat*x[i])*(
        1+(di0[i])*exp(beta1hat*x[i+100])+
         (di0[i])*(ei0[i])*exp(beta1hat*x[i+100])*exp(beta1hat*x[i+200])) 
    obj = obj + a
  }
  
  return(obj)
}

eqn4 <-function(x)
{
  sumt = 0
  for (i in 1:300) {
    c=x[i]
    sumt = sumt + c
  }
  return(sumt-111.54) #64.63122
  
}

lb <- rep(((0)-mu[1])/sd[1],300)
ub <- rep((max(Bupdiff$Bupdiff[101:600])-(mean(Bupdiff$Bupdiff[101:600])))/
            sd(Bupdiff$Bupdiff[101:600]),300)

x0 <- c(rep(0,300))
fn4(c(rep((0-mu)/sd,300)))

res3yr_num<-auglag(x0,
                   fn4,
                   heq = eqn4,
                   lower          = lb,
                   upper          = ub,
                   localsolver = "mma",
                   localtol=1e-05)


fn4r <- function(x) 
{
  
  
  obj = 0
  #  grad=0
  for (i in 1:100) {
    
    a= ci[i]*exp(beta1hat*x[i])*(
      1+(di0[i])*exp(beta1hat*x[i+100])+
        (di0[i])*(ei0[i])*exp(beta1hat*x[i+100])*exp(beta1hat*x[i+200])) 

    obj = obj +(a)/Exp[i,6]    
  }
  
  return(obj/100)
}
eqn4r <-function(x)
{
  sumt = 0
  for (i in 1:300) {
    c=x[i]
    sumt = sumt + c
  }
  return(sumt-111.54) #64.63122
  
}


fn4r(c(rep((0-mu)/sd,300)))

#three year rate version
res3yr_rate<-auglag(x0,
                      fn4r,
                      heq = eqn4r,
                      lower          = lb,
                      upper          = ub,
                      localsolver = "mma",
                      localtol=1e-05)

p=4
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))#theta for each county each quarter
ii=2
nn=3

# three years
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  
  
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  
  
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  
  mm<-seq(101,501,by=100)#101-3701 q2-q38
  
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  
  
  postX3900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX4000<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39<-rep(postbeta0,100)+postX3900[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40<-rep(postbeta0,100)+postX4000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41<-rep(postbeta0,100)+postX4100[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta39)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta40)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta41)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  

  
  postX3900x<-cbind(res3yr_num$par[1:100], postXmm[,ii:nn],rep(7 ,100))
  postX4000x<-cbind(res3yr_num$par[101:200], postXmm[,ii:nn],rep(8 ,100))
  postX4100x<-cbind(res3yr_num$par[201:300], postXmm[,ii:nn],rep(9 ,100))
  xbeta39x<-rep(postbeta0,100)+postX3900x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40x<-rep(postbeta0,100)+postX4000x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41x<-rep(postbeta0,100)+postX4100x[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta39x)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta40x)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta41x)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  
  
  
}


#probability for 2022
statesum<-apply(thetasave[,10,],2, sum)
statebase<-apply(thetasave[,7,],2,sum)
length(which((statebase-statesum)>0))
#probability for 2022, 23, 24
statesum<-apply(thetasave[,10,],2, sum)+apply(thetasave[,11,],2, sum)+apply(thetasave[,12,],2, sum)
statebase<-apply(thetasave[,7,],2,sum)+apply(thetasave[,8,],2,sum)+apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))
#probability for 2024
statesum<-apply(thetasave[,12,],2, sum)
statebase<-apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))


statesum<-apply(thetasave[,12,],2, sum)
statebase<-0.8*apply(thetasave[,9,],2,sum)
length(which((statebase-statesum)>0))


statesum<-apply(thetasave[,12,],2, sum)
statebase<-0.8*apply(thetasave[,6,],2,sum)
length(which((statebase-statesum)>0))


# three year rate 
for (ll in 1:n) { 
  
  l=1*ll
  epsilon<- sampleused[[1]][l,which(names(sampleused[[1]][1,])=='epsilon[1, 1]')
                            :which(names(sampleused[[1]][1,])=='epsilon[100, 1]')]
  s1<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s1[1]')
                      :(which(names(sampleused[[1]][1,])=='s1[100]'))]
  s2<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='s2[1]')
                      :(which(names(sampleused[[1]][1,])=='s2[100]'))]
  
  postbeta0<-sampleused[[1]][l,'beta0']
  
  
  postbeta<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='beta[1]')
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+3)]
  
  
  postalpha<-sampleused[[1]][l,which(names(sampleused[[1]][1,])=='alpha[1]')
                             :(which(names(sampleused[[1]][1,])=='alpha[1]')+2)]
  postalpha0<-sampleused[[1]][l,'alpha0']
  
  postX<-matrix(sampleused[[1]][l,1:(p*600)],600, )
  
  #q1
  postX100<-postX[1:100,]
  xalpha<-rep(postalpha0,100)+postX100[,1:3]%*%postalpha+rep(s1,1)+epsilon
  thetasave[,1,ll]<-exp(xalpha)*Exp[,1]
  theta_last<-exp(xalpha)*Exp[,1]
  
  mm<-seq(101,501,by=100)#101-3701 q2-q38
  
  for (m in mm) {
    postXmm<-postX[m:(m+99),]
    xbeta<-rep(postbeta0,100)+postXmm[,1:4]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  
  
  postX3900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX4000<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX4100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta39<-rep(postbeta0,100)+postX3900[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40<-rep(postbeta0,100)+postX4000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41<-rep(postbeta0,100)+postX4100[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta39)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta40)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta41)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  #postX3900x<-cbind(resrate_1$par, postXmm[,ii:nn],rep(7 ,100))
  
  # postX3900x<-cbind(bun_rate2000$par, postXmm[,ii:nn],rep(7 ,100))
  # postX3900x<-cbind(resnum_2000$par, postXmm[,ii:nn],rep(7 ,100))
  
  
  postX3900x<-cbind(res3yr_rate$par[1:100], postXmm[,ii:nn],rep(7 ,100))
  postX4000x<-cbind(res3yr_rate$par[101:200], postXmm[,ii:nn],rep(8 ,100))
  postX4100x<-cbind(res3yr_rate$par[201:300], postXmm[,ii:nn],rep(9 ,100))
  xbeta39x<-rep(postbeta0,100)+postX3900x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta40x<-rep(postbeta0,100)+postX4000x[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta41x<-rep(postbeta0,100)+postX4100x[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta39x)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta40x)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta41x)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  
  
  
}


#probability for 2022
statesumR<-apply(apply(thetasave[,10,], 2, '/',Exp[,6]),2, sum)
statebaseR<-apply(apply(thetasave[,7,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))

#probability for 2022, 23, 24
statesum<-apply(apply(thetasave[,10,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,11,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,12,], 2, '/',Exp[,6]),2, sum)
statebase<-apply(apply(thetasave[,7,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,8,], 2, '/',Exp[,6]),2, sum)+apply(apply(thetasave[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))


#probability for 2024
statesumR<-apply(apply(thetasave[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-0.8*apply(apply(thetasave[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))


#probability for 2024
statesumR<-apply(apply(thetasave[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-0.8*apply(apply(thetasave[,6,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))






low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)

thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(thetasave[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(thetasave[i,j,])[1]
    high[i,j]= hdi(thetasave[i,j,])[2]
    mean[i,j]=mean(thetasave[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(thetasave[i,(9+j),])[1]
    prelow[i,j]= hdi(thetasave[i,(9+j),])[2]
    
    Predicted[i,j]=mean(thetasave[i,(9+j),],na.rm = T)
  }
}





estimathreeyear<-data.frame(time=c(2016:2024),Estimated=colSums(mean), 'Upper Bound'=colSums(high),
                            'Lower Bound'=colSums(low), Observed=c(colSums(y[,1:6],na.rm = T),rep(NA,3)),
                            Predicted=c(rep(NA,6),colSums(Predicted)),
                            prehigh=c(rep(NA,6),colSums(prehigh)),
                            prelow=c(rep(NA,6),colSums(prelow)))
newexp<-cbind(Exp,Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6])

estimathreeyearrate<-data.frame(time=c(2016:2024),Estimated=colMeans(mean/newexp[,1:9]), 'Upper Bound'=colMeans(high/newexp[,1:9]),
                                'Lower Bound'=colMeans(low/newexp[,1:9]), Observed=c(colMeans(y[,1:6]/Exp,na.rm = T),rep(NA,3)),
                                Predicted=c(rep(NA,6),colMeans(Predicted/newexp[,10:12])),
                                prehigh=c(rep(NA,6),colMeans(prehigh/newexp[,10:12])),
                                prelow=c(rep(NA,6),colMeans(prelow/newexp[,10:12])))




ggplot(data=estimathreeyear,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Number")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())

ggplot(data=estimathreeyearrate,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Rate")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())

BUPincr<-res3yr_num$par[1:100]
BUPincr<-res3yr_num$par[101:200]
BUPincr<-res3yr_num$par[201:300]
BUPincr<-res3yr_rate$par[1:100]
BUPincr<-res3yr_rate$par[101:200]
BUPincr<-res3yr_rate$par[201:300]



#map 
usa_counties = map_data("county")
nc = subset(usa_counties, region == "north carolina") 
library(ggpubr)


result4plotactnum<-data.frame(cbind(num=BUPincr*sd[1]+mu[1],subregion=final.data$Place[1:100]))
BUPincr1<-as.numeric(result4plotactnum$num)

result4plot<-data.frame(cbind(BUPincr1,subregion=final.data$Place[1:100]))
result4plot$subregion<-tolower(result4plot$subregion)
map.data = merge(nc,result4plot,by.x="subregion",by.y="subregion")
map.data$BUPincr1<-as.numeric(map.data$BUPincr1)
allocation2024<-ggplot()+
  geom_polygon(data=map.data,aes(x=long,y=lat,group=subregion,fill=BUPincr1),color='black',alpha=.8,size=.3)+
  scale_fill_gradient2(name="Increase",low='white',high='darkblue')+
  geom_point(data=cityct,aes(x=long, y=lat,size=Population),shape = 17,color="pink") +
   coord_map()+
  theme_void()
ggpubr::ggarrange(allocation2022,allocation2023,allocation2024,
                  ncol = 1,common.legend=T,legend ="right")


ggpubr::ggarrange(Rallocation2022,Rallocation2023,Rallocation2024,
                  ncol=1,common.legend=T,legend ="right")


