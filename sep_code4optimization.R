library(nloptr)
library(xtable)

rm(list = ls())
set.seed(12345)




################################################################################
#####data preparation###########################################################
################################################################################

setwd("/Users/qianyu/Documents/bayes/data")
final.data<-read.csv(file = "final.data.csv")

##################################
#Exp#Exp#Exp#Exp#Exp##Exp#Exp#Exp
Exp_yearly<-final.data$Rate.Denom/1000
Exp =  matrix(Exp_yearly,nrow = 100,byrow = F)
###################################
#Bup#Bup#Bup#Bup#Bup##Bup#Bup#Bup#
Bupdiff <-final.data%>%
  dplyr:: arrange(Year) %>%
  dplyr:: group_by(Place) %>%
  dplyr::mutate(
    Bupdiff =  Bup.patient-c(0, dplyr::lag(Bup.patient)[2:6])
  )%>%
  ungroup()
Bupdiff<-as.data.frame(Bupdiff)

result_sep172<-readRDS(file = "result_sep172.rda")
sampleused<-result_sep172

###################################
#map#map#map#map#map#map#map#map#ma
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




################################################################################
#####One year optimiztion#######################################################
################################################################################

###################################
#One year number ##################
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
#fn1(rep((0-mu)/sd,100)) 
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


numTABLE<- data.frame(cbind(result4plot_num[1:50,],result4plot_num[51:100,]))
colnames(numTABLE) <- c('subregion','value','subregion','value')
xtable(numTABLE, type = "latex", file = "result4plot.tex")



#line plot
p=4
n=10000
ii=2
nn=3
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))
#theta(poisson mean) for each county each quarter,
#thetasave[,7:9,] are baseline 10:12 are prediction with optimal solution
s.y.new=array(NA, dim = c(100,12,n))# the same for s.y.new

set.seed(12345)
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
  postX700<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX800<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta7<-rep(postbeta0,100)+postX700[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta8<-rep(postbeta0,100)+postX800[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta9<-rep(postbeta0,100)+postX900[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta7)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta8)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta9)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  

  postX1000<-cbind(resnum_6743$par, postXmm[,ii:nn],rep(7 ,100))
  postX1100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX1200<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta10<-rep(postbeta0,100)+postX1000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta11<-rep(postbeta0,100)+postX1100[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta12<-rep(postbeta0,100)+postX1200[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta10)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta11)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta12)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  
  # post.theta<-c(thetasave[,1,ll],thetasave[,2,ll],thetasave[,3,ll],
  #               thetasave[,4,ll],thetasave[,5,ll],thetasave[,6,ll],
  #               thetasave[,7,ll],thetasave[,8,ll],thetasave[,9,ll],
  #               thetasave[,10,ll],thetasave[,11,ll],thetasave[,12,ll])
  
  for (k in 1:12) {
    s.y.new[,k,ll]<- rpois(length(thetasave[,k,ll]),thetasave[,k,ll])
  }

  

}



# #probability for 2022
# statesum<-apply(thetasave[,10,],2, sum)
# statebase<-apply(thetasave[,7,],2,sum)
# length(which((statebase-statesum)>0))
# #probability for 2022, 23, 24
# statesum<-apply(thetasave[,10,],2, sum)+apply(thetasave[,11,],2, sum)+apply(thetasave[,12,],2, sum)
# statebase<-apply(thetasave[,7,],2,sum)+apply(thetasave[,8,],2,sum)+apply(thetasave[,9,],2,sum)
# length(which((statebase-statesum)>0))
# 


#2024
statesum<-apply(s.y.new[,12,],2, sum)
statebase<-apply(s.y.new[,9,],2,sum)
length(which((statebase-statesum)>0))

#2024 0.8
statesum<-apply(s.y.new[,12,],2, sum)
statebase<-0.8*apply(s.y.new[,9,],2,sum)
length(which((statebase-statesum)>0))





low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)


thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(s.y.new[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(s.y.new[i,j,])[1]
    high[i,j]= hdi(s.y.new[i,j,])[2]
    mean[i,j]=mean(s.y.new[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(s.y.new[i,(9+j),])[1]
    prelow[i,j]= hdi(s.y.new[i,(9+j),])[2]
    Predicted[i,j]=mean(s.y.new[i,(9+j),],na.rm = T)
  }
}

# for (i in 1:100) {
#   for (j in 1:3) {
#     prehigh[i,j]=hdi(s.y.new[i,j,])[1]
#     prelow[i,j]= hdi(s.y.new[i,j,])[2]
#     Predicted[i,j]=mean(s.y.new[i,j,],na.rm = T)
#   }
# }




y=matrix(as.numeric(final.data$death),nrow=100,byrow=F)

estima2<-data.frame(time=c(2016:2024),Estimated=colSums(mean), 'Upper Bound'=colSums(high),
                    'Lower Bound'=colSums(low), Observed=c(colSums(y[,1:6],na.rm = T),rep(NA,3)),
                    Predicted=c(rep(NA,6),colSums(Predicted)),
                    prehigh=c(rep(NA,6),colSums(prehigh)), 
                    prelow=c(rep(NA,6),colSums(prelow)))



COL=c("Estimated" = "black",
     # 'Credible Bound' = "steelblue",
    #  'Credible Bound' = "steelblue",
      "Observed" = "red",
      "Predicted" = "darkgreen"
     # "Predicted Bound" = "darkgreen",
    #  "Predicted Bound" = "darkgreen" 
)
oneyearnum<-ggplot(data=estima2,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_ribbon(aes(ymin=Lower.Bound, ymax=Upper.Bound),fill="grey70", alpha=0.2)+

  #geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  #geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_ribbon(aes(ymin=prelow, ymax=prehigh),fill="green", alpha=0.2)+
  
  # geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  # geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Number")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())




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




###################################
#One year rate ##################
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


fn2 <- function(x)  #one year rate
{
  obj = 0
  for (i in 1:100) {
    a=(ci[i]*exp(beta1hat*x[i]))/Exp[i,6]
    obj = obj + a 
  }
  
  return(obj/100)
}
fn2(rep((0-mu)/sd,100))

eqn2 <-function(x)#6743
{
  sumt = 0
  for (i in 1:100) {
    c=x[i]
    sumt = sumt + c
  }
  return(sumt-37.18)#mean(c(5669,8566,5996)), (6743.667-100*mu)/sd
}
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



rateTABLE<- data.frame(cbind(result4plot_rate[1:50,],result4plot_rate[51:100,]))
colnames(rateTABLE) <- c('subregion','value','subregion','value')
xtable(rateTABLE, type = "latex", file = "result4plot.tex")




p=4
n=10000
ii=2
nn=3
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))
s.y.new=array(NA, dim = c(100,12,n))
set.seed(12345)

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
  
  
  
  postX700<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX800<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta7<-rep(postbeta0,100)+postX700[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta8<-rep(postbeta0,100)+postX800[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta9<-rep(postbeta0,100)+postX900[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta7)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta8)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta9)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  
  
  postX1000<-cbind(resrate_6743$par, postXmm[,ii:nn],rep(7 ,100))
  postX1100<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX1200<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta10<-rep(postbeta0,100)+postX1000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta11<-rep(postbeta0,100)+postX1100[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta12<-rep(postbeta0,100)+postX1200[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta10)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta11)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta12)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
 # s.y.new[,1,ll]<- rpois(length(thetasave[,10,ll]),thetasave[,10,ll])
  #s.y.new[,2,ll]<- rpois(length(thetasave[,11,ll]),thetasave[,11,ll])
  #s.y.new[,3,ll]<- rpois(length(thetasave[,12,ll]), thetasave[,12,ll])
  
  for (k in 1:12) {
    s.y.new[,k,ll]<- rpois(length(thetasave[,k,ll]),thetasave[,k,ll])
  }
  
  
  
}


# #2021
# statesumR<-apply(apply(s.y.new[,3,], 2, '/',Exp[,6]),2, sum)
# statebaseR<-apply(apply(thetasave[,6,], 2, '/',Exp[,6]),2, sum)
# length(which((statebaseR-statesumR)>0))
# 
# 
# statesumR<-apply(apply(s.y.new[,3,], 2, '/',Exp[,6]),2, sum)
# statebaseR<-0.8*apply(apply(thetasave[,6,], 2, '/',Exp[,6]),2, sum)
# length(which((statebaseR-statesumR)>0))

#2024
statesumR<-apply(apply(s.y.new[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-apply(apply(s.y.new[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))


statesumR<-apply(apply(s.y.new[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-0.8*apply(apply(s.y.new[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))





low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)


thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(s.y.new[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(s.y.new[i,j,])[1]
    high[i,j]= hdi(s.y.new[i,j,])[2]
    mean[i,j]=mean(s.y.new[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(s.y.new[i,(9+j),])[1]
    prelow[i,j]= hdi(s.y.new[i,(9+j),])[2]
    Predicted[i,j]=mean(s.y.new[i,(9+j),],na.rm = T)
  }
}

# for (i in 1:100) {
#   for (j in 1:3) {
#     prehigh[i,j]=hdi(s.y.new[i,j,])[1]
#     prelow[i,j]= hdi(s.y.new[i,j,])[2]
#     Predicted[i,j]=mean(s.y.new[i,j,],na.rm = T)
#   }
# }




y=matrix(as.numeric(final.data$death),nrow=100,byrow=F)
newexp<-cbind(Exp,Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6])
# use one year rate, then repeat line 379-413
estima_r<-data.frame(time=c(2016:2024),Estimated=colMeans(mean/newexp[,1:9]), 'Upper Bound'=colMeans(high/newexp[,1:9]),
                     'Lower Bound'=colMeans(low/newexp[,1:9]), Observed=c(colMeans(y[,1:6]/Exp,na.rm = T),rep(NA,3)),
                     Predicted=c(rep(NA,6),colMeans(Predicted/newexp[,10:12])),
                     prehigh=c(rep(NA,6),colMeans(prehigh/newexp[,10:12])), 
                     prelow=c(rep(NA,6),colMeans(prelow/newexp[,10:12])))




COL=c("Estimated" = "black",
      #'Credible Bound' = "steelblue",
      #'Credible Bound' = "steelblue",
      "Observed" = "red",
      "Predicted" = "darkgreen"
      #,
      #"Predicted Bound" = "darkgreen",
     # "Predicted Bound" = "darkgreen" 
)

oneyearrate<-ggplot(data=estima_r,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_ribbon(aes(ymin=Lower.Bound, ymax=Upper.Bound),fill="grey70", alpha=0.2)+
  
  #geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  #geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_ribbon(aes(ymin=prelow, ymax=prehigh),fill="green", alpha=0.2)+
  
  # geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  # geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Rate")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())



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


#number of counties receiving 0
dim(result4num[result4num$num==0,])
#77  2
dim(result4num[result4rate$num==0,])
#19  2

#five number sum of results
rate<-summary(result4plot_rate$num)
num<-summary(result4plot_num$num)
rbind(rate,num)

rateTABLE<- data.frame(rbind(rate,num))
colnames(rateTABLE) <- c("Min.", "1st Qu." ,"Median" , "Mean" ,"3rd Qu.", "Max.")
xtable(rateTABLE, type = "latex", file = "result4plot.tex")


################################################################################
#####Three year optimiztion#####################################################
################################################################################

###############################
#####Three year number#########
###############################
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


p=4
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
thetasave=array(NA, dim = c(100,12,n))#theta for each county each quarter
s.y.new=array(NA, dim = c(100,12,n))
ii=2
nn=3

set.seed(12345)
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
  
  
  
  postX700<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX800<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta7<-rep(postbeta0,100)+postX700[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta8<-rep(postbeta0,100)+postX800[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta9<-rep(postbeta0,100)+postX900[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta7)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta8)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta9)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  
  
  postX1000<-cbind(res3yr_num$par[1:100], postXmm[,ii:nn],rep(7 ,100))
  postX1100<-cbind(res3yr_num$par[101:200], postXmm[,ii:nn],rep(8 ,100))
  postX1200<-cbind(res3yr_num$par[201:300], postXmm[,ii:nn],rep(9 ,100))
  xbeta10<-rep(postbeta0,100)+postX1000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta11<-rep(postbeta0,100)+postX1100[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta12<-rep(postbeta0,100)+postX1200[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta10)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta11)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta12)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
 # s.y.new[,1,ll]<- rpois(length(thetasave[,10,ll]),thetasave[,10,ll])
  #s.y.new[,2,ll]<- rpois(length(thetasave[,11,ll]),thetasave[,11,ll])
  #s.y.new[,3,ll]<- rpois(length(thetasave[,12,ll]), thetasave[,12,ll])
  
  for (k in 1:12) {
    s.y.new[,k,ll]<- rpois(length(thetasave[,k,ll]),thetasave[,k,ll])
  }

}


# #probability for 2022
# statesum<-apply(thetasave[,10,],2, sum)
# statebase<-apply(thetasave[,7,],2,sum)
# length(which((statebase-statesum)>0))
# #probability for 2022, 23, 24
# statesum<-apply(thetasave[,10,],2, sum)+apply(thetasave[,11,],2, sum)+apply(thetasave[,12,],2, sum)
# statebase<-apply(thetasave[,7,],2,sum)+apply(thetasave[,8,],2,sum)+apply(thetasave[,9,],2,sum)
# length(which((statebase-statesum)>0))
# #probability for 2024
# statesum<-apply(thetasave[,12,],2, sum)
# statebase<-apply(thetasave[,9,],2,sum)
# length(which((statebase-statesum)>0))
# 


#2024
statesum<-apply(s.y.new[,12,],2, sum)
statebase<-apply(s.y.new[,9,],2,sum)
length(which((statebase-statesum)>0))

#2024 0.8
statesum<-apply(s.y.new[,12,],2, sum)
statebase<-0.8*apply(s.y.new[,9,],2,sum)
length(which((statebase-statesum)>0))







low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)

thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(s.y.new[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(s.y.new[i,j,])[1]
    high[i,j]= hdi(s.y.new[i,j,])[2]
    mean[i,j]=mean(s.y.new[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(s.y.new[i,(9+j),])[1]
    prelow[i,j]= hdi(s.y.new[i,(9+j),])[2]

    Predicted[i,j]=mean(s.y.new[i,(9+j),],na.rm = T)
  }
}

# 
# for (i in 1:100) {
#   for (j in 1:3) {
#     prehigh[i,j]=hdi(s.y.new[i,(j),])[1]
#     prelow[i,j]= hdi(s.y.new[i,(j),])[2]
#     Predicted[i,j]=mean(s.y.new[i,(j),],na.rm = T)
#   }
# }


estimathreeyear<-data.frame(time=c(2016:2024),Estimated=colSums(mean), 'Upper Bound'=colSums(high),
                            'Lower Bound'=colSums(low), Observed=c(colSums(y[,1:6],na.rm = T),rep(NA,3)),
                            Predicted=c(rep(NA,6),colSums(Predicted)),
                            prehigh=c(rep(NA,6),colSums(prehigh)),
                            prelow=c(rep(NA,6),colSums(prelow)))

estimathreeyearplot<-ggplot(data=estimathreeyear,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_ribbon(aes(ymin=Lower.Bound, ymax=Upper.Bound),fill="grey70", alpha=0.2)+
  #geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  #geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_ribbon(aes(ymin=prelow, ymax=prehigh),fill="green", alpha=0.2)+
  #geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  #geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+
  labs(x = "Year", y="Number")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())



#map 
BUPincr<-res3yr_num$par[1:100]
BUPincr<-res3yr_num$par[101:200]
BUPincr<-res3yr_num$par[201:300]

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



###############################
#####Three year number#########
###############################
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
s.y.new=array(NA, dim = c(100,12,n))
ii=2
nn=3

set.seed(12345)
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
  
  
  
  postX700<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(7 ,100))
  postX800<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(8 ,100))
  postX900<-cbind(rep(((0)-mu[1])/sd[1],100), postXmm[,ii:nn],rep(9 ,100))
  xbeta7<-rep(postbeta0,100)+postX700[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta8<-rep(postbeta0,100)+postX800[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta9<-rep(postbeta0,100)+postX900[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,7,ll]<-exp(xbeta7)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,8,ll]<-exp(xbeta8)*thetasave[,7,ll]*Exp[,6]/Exp[,6]
  thetasave[,9,ll]<-exp(xbeta9)*thetasave[,8,ll]*Exp[,6]/Exp[,6]
  
  
  
  postX1000<-cbind(res3yr_rate$par[1:100], postXmm[,ii:nn],rep(7 ,100))
  postX1100<-cbind(res3yr_rate$par[101:200], postXmm[,ii:nn],rep(8 ,100))
  postX1200<-cbind(res3yr_rate$par[201:300], postXmm[,ii:nn],rep(9 ,100))
  xbeta10<-rep(postbeta0,100)+postX1000[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta11<-rep(postbeta0,100)+postX1100[,1:(p)]%*%postbeta+rep(s2,1)
  xbeta12<-rep(postbeta0,100)+postX1200[,1:(p)]%*%postbeta+rep(s2,1)
  thetasave[,10,ll]<-exp(xbeta10)*thetasave[,6,ll]*Exp[,6]/Exp[,6]
  thetasave[,11,ll]<-exp(xbeta11)*thetasave[,10,ll]*Exp[,6]/Exp[,6]
  thetasave[,12,ll]<-exp(xbeta12)*thetasave[,11,ll]*Exp[,6]/Exp[,6]
  #s.y.new[,1,ll]<- rpois(length(thetasave[,10,ll]),thetasave[,10,ll])
  #s.y.new[,2,ll]<- rpois(length(thetasave[,11,ll]),thetasave[,11,ll])
  #s.y.new[,3,ll]<- rpois(length(thetasave[,12,ll]), thetasave[,12,ll])
  
  for (k in 1:12) {
    s.y.new[,k,ll]<- rpois(length(thetasave[,k,ll]),thetasave[,k,ll])
  }
}


#2024
statesumR<-apply(apply(s.y.new[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-apply(apply(s.y.new[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))


statesumR<-apply(apply(s.y.new[,12,], 2, '/',Exp[,6]),2, sum)
statebaseR<-0.8*apply(apply(s.y.new[,9,], 2, '/',Exp[,6]),2, sum)
length(which((statebaseR-statesumR)>0))




low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)
Predicted=matrix(NA, nrow = 100,ncol = 3)
prehigh=matrix(NA, nrow = 100,ncol = 3)
prelow=matrix(NA, nrow = 100,ncol = 3)

thetamean9<-matrix(0,100,9)
for (i in 1:100) {
  thetamean9[i,]= rowMeans(s.y.new[i,1:9,],na.rm = T)
  
}

for (i in 1:100) {
  for (j in 1:9) {
    low[i,j]= hdi(s.y.new[i,j,])[1]
    high[i,j]= hdi(s.y.new[i,j,])[2]
    mean[i,j]=mean(s.y.new[i,j,],na.rm = T)
  }
}

for (i in 1:100) {
  for (j in 1:3) {
    prehigh[i,j]=hdi(s.y.new[i,(9+j),])[1]
    prelow[i,j]= hdi(s.y.new[i,(9+j),])[2]
    Predicted[i,j]=mean(s.y.new[i,(9+j),],na.rm = T)
  }
}

# 
# for (i in 1:100) {
#   for (j in 1:3) {
#     prehigh[i,j]=hdi(s.y.new[i,(j),])[1]
#     prelow[i,j]= hdi(s.y.new[i,(j),])[2]
#     
#     Predicted[i,j]=mean(s.y.new[i,(j),],na.rm = T)
#   }
# }

newexp<-cbind(Exp,Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6])

estimathreeyearrate<-data.frame(time=c(2016:2024),Estimated=colMeans(mean/newexp[,1:9]), 'Upper Bound'=colMeans(high/newexp[,1:9]),
                                'Lower Bound'=colMeans(low/newexp[,1:9]), Observed=c(colMeans(y[,1:6]/Exp,na.rm = T),rep(NA,3)),
                                Predicted=c(rep(NA,6),colMeans(Predicted/newexp[,10:12])),
                                prehigh=c(rep(NA,6),colMeans(prehigh/newexp[,10:12])),
                                prelow=c(rep(NA,6),colMeans(prelow/newexp[,10:12])))




estimathreeyearratepp<-ggplot(data=estimathreeyearrate,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_ribbon(aes(ymin=Lower.Bound, ymax=Upper.Bound),fill="grey70", alpha=0.2)+
  #geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  #geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  geom_line(aes( y=Predicted ,color ="Predicted" ),size=1)+
  labs(color = '')+
  scale_color_manual(values = COL)+
  
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  geom_point(aes(x=time, y=Predicted ),shape = 24,fill="green",size=2) +
  geom_ribbon(aes(ymin=prelow, ymax=prehigh),fill="green", alpha=0.2)+
  labs(x = "Year", y="Rate")+
  
  #geom_line(aes( y=prehigh,color = "Predicted Bound") ,linetype="twodash" )+
  #geom_line(aes( y=prelow ,color = "Predicted Bound"), linetype="twodash" )+  labs(x = "Year", y="Rate")+
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
library(ggpubr)

BUPincr<-res3yr_rate$par[1:100]
BUPincr<-res3yr_rate$par[101:200]
BUPincr<-res3yr_rate$par[201:300]



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

res1001_1year<-list(resnum_6743,resrate_6743,estima2,estima_r)
saveRDS(res1001_1year,file = "res1001_1year.rda")

res1001_3year<-list(res3yr_num,res3yr_rate,estimathreeyear,estimathreeyearrate)
saveRDS(res1001_3year,file = "res1001_3year.rda")

