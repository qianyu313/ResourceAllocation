
library(spData)
library('spdep')
library('maptools')
library('shapefiles')
library(ggplot2)
library('ggmap')
library('rgdal')
library('Rmisc')#
library('gridExtra')
library('grid')
library(raster)
library(rgeos)
library(nimble)
library(coda)
library(smfsb)#
library(HDInterval)
library(jsonlite)
library(dplyr)
library(tidyr)
library(questionr)
library(MASS)
library(rjags)#
library("matrixStats")
library(nloptr)


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

################################################################################
#####Model #####################################################################
################################################################################
Bupdiff_scale<-c(scale(Bupdiff$Bupdiff[1:100])[1:100],
                 scale(Bupdiff$Bupdiff[101:600])) #1:100 are x_bup in t=1, 101:600 are different from previous year
unemploy_yearly=as.numeric(final.data$unemploy)/(Bupdiff$Rate.Denom/1000)  #Rate.Denom is population 
unemploy_scale<-scale(unemploy_yearly[1:600])[1:600]
otp_yearly_scale<-scale(final.data$Otps[1:600])[1:600]
X<-cbind(Bupdiff_scale,unemploy_scale ,otp_yearly_scale)
year<-rep(c(1:6), each = 100)
X<-cbind(X,year)
p = dim(X)[2]
y=matrix(as.numeric(final.data$death),nrow=100,byrow=F)
#initial values for coefficients in t=1, and t>1 models
m1<-glm(y[,1]~X[1:100,1:3],family="poisson",offset=log(Exp[,1]))
m2<-glm(c(y[,2],y[,3],y[,4],y[,5],y[,6])
        ~X[101:600,1:p],family="poisson",
        offset=log(c(Exp[,2],Exp[,3],Exp[,4],Exp[,5],Exp[,6])))

###adjacency matrix  for icar
nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO) 
W.nb <- poly2nb(nc,row.names =  rownames(nc$NAME))
nbInfo <- nb2WB(W.nb)
nregions=100
p = dim(X)[2]
T = 6
n = 100

code2 <- nimbleCode({
  for (t in 1:1){
    for (i in 1:n) {
      y[i,t] ~ dpois(theta[i,t])
      log(theta[i,t]) <- log(E[i,t]) + alpha0+ inprod(X[((t-1)*n+i),1:(p-1)],alpha[1:(p-1)])+ s1[i]+epsilon[i,t]
      epsilon[i,t] ~ dnorm(0,sig0.inv)
      
    }
  }
  for (t in 2:T){
    for(i in 1:n){
      y[i,t] ~ dpois(theta[i,t])
      log(theta[i,t]) <- log(theta[i,t-1])+log(E[i,t])-log(E[i,t-1])+beta0+inprod(X[((t-1)*n+i),1:p],beta[1:p])+s2[i]
    }
  }
  
  s1[1:n] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n], tau1, zero_mean = 1)
  s2[1:n] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n], tau2, zero_mean = 1)
  
  beta0 ~ dnorm(muInt,sig.beta.inv.Int)
  alpha0 ~ dnorm(mu0Int,sig.alpha.inv.Int)
  
  for(j in 1:(p-1)){
    alpha[j] ~ dnorm(0,1000)
  }
  for(j in 1:p){
    beta[j] ~ dnorm(0,1000)
  }
  
  tau1~  dinvgamma(.1,.1)
  tau2~  dinvgamma(.1,.1)
  muInt ~ dnorm(0,.001)
  mu0Int ~ dnorm(0,.001)
  sig0.inv ~  dinvgamma(.1,.1)
  sig.beta.inv.Int ~ dinvgamma(.1,.1)
  sig.alpha.inv.Int ~ dinvgamma(.1,.1)
  
})
constants <-list(E=Exp,T=T,p=p,
                 #bd=bd,
                 n = nregions, L = length(nbInfo$adj),
                 adj = nbInfo$adj, weights = nbInfo$weights,
                 num = nbInfo$num)
data <- list(y=y,X=X)
inits2<-list(
  alpha=m1$coefficients[2:p],beta=m2$coefficients[2:(p+1)],
  alpha0=m1$coefficients[1],beta0=m2$coefficients[1],
  mu0Int=0,
  sig.alpha.inv.Int=.5,
  muInt=0, 
  sig.beta.inv.Int=.5,
  s1= rnorm(nregions),s2= rnorm(nregions),tau1=0.5,tau2=0.5
  ,sig0.inv=.5
)
model <- nimbleModel(code2, constants = constants, data = data,inits=inits2)
cModel <- compileNimble(model,showCompilerOutput = TRUE)
conf <- configureMCMC(model, monitors = c('alpha','beta',
                                          'alpha0','beta0','epsilon','X',
                                          's1','s2','y','tau1','tau2')
                      ,enableWAIC = T)
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel) 
result_sep172<- runMCMC( cMCMC, niter = 1500000, nburnin = 1000000,summary = F
                         ,samplesAsCodaMCMC = T, nchains = 2,thin=50)

#saveRDS(result_sep172,file = "result_sep172.rda")
#result_sep172<-readRDS(file = "result_sep172.rda")




################################################################################
#####Model result###############################################################
################################################################################

sampleused<-result_sep172

#t>1
gelman.diag(sampleused[,which(names(sampleused[[1]][1,])=='beta[1]'):
                         (which(names(sampleused[[1]][1,])=='beta[1]')+p)])
hdi(sampleused[,which(names(sampleused[[1]][1,])=='beta[1]'):
                 (which(names(sampleused[[1]][1,])=='beta[1]')+p)],credMass = 0.95)
colMeans(sampleused[,which(names(sampleused[[1]][1,])=='beta[1]'):
                      (which(names(sampleused[[1]][1,])=='beta[1]')+p)][[1]])

#proportion of density < 0
a=sampleused[[1]][,which(names(sampleused[[1]][1,])=='beta[1]')]
length(a[a<0])/10000


#t=1
gelman.diag(sampleused[,which(names(sampleused[[1]][1,])=='alpha[1]'):
                         (which(names(sampleused[[1]][1,])=='alpha[1]')+p)])
hdi(sampleused[,which(names(sampleused[[1]][1,])=='alpha[1]'):
                 (which(names(sampleused[[1]][1,])=='alpha[1]')+p)],credMass = 0.95)
colMeans(sampleused[,which(names(sampleused[[1]][1,])=='alpha[1]'):
                      (which(names(sampleused[[1]][1,])=='alpha[1]')+p)][[1]])

#trace plot
plot(sampleused[,which(names(sampleused[[1]][1,])=='beta[1]'):
                  (which(names(sampleused[[1]][1,])=='beta[1]')+p)],density = FALSE)

#map for county level intercept(mean and sd)
mean<-colMeans(sampleused[,which(names(sampleused[[1]][1,])=='s2[1]'):
                            (which(names(sampleused[[1]][1,])=='s2[1]')+99)][[1]])

usa_counties = map_data("county")
nc = subset(usa_counties, region == "north carolina") 
result4plot<-data.frame(cbind(value=mean,subregion=combined_data$Place[1:100]))
result4plot$subregion<-tolower(result4plot$subregion)
map.data<- merge(nc,result4plot,by.x="subregion",by.y="subregion")
map.data$value<-as.numeric(map.data$value)

ggplot()+
  geom_polygon(data=map.data,aes(x=long,y=lat,group=subregion,fill=value),color='black',alpha=.8,size=.3)+
  scale_fill_gradient2(name="Mean",low='blue',high='limegreen', mid = 'white',midpoint = 0)+
  coord_map()+
  theme_void()+
  theme(legend.text=element_text(size=12))+
  theme(legend.title =element_text(size=14))

sd<-colSds(sampleused[,which(names(sampleused[[1]][1,])=='s2[1]'):
                        (which(names(sampleused[[1]][1,])=='s2[1]')+99)][[1]])

result4plot<-data.frame(cbind(value=sd,subregion=combined_data$Place[1:100]))
result4plot$subregion<-tolower(result4plot$subregion)
map.data<- merge(nc,result4plot,by.x="subregion",by.y="subregion")
map.data$value<-as.numeric(map.data$value)

ggplot()+
  geom_polygon(data=map.data,aes(x=long,y=lat,group=subregion,fill=value),color='black',alpha=.8,size=.3)+
  scale_fill_gradient2(name="SD",low='white',high='darkred')+
  coord_map()+
  theme_void()+
  theme(legend.text=element_text(size=12))+
  theme(legend.title =element_text(size=14))



n=10000#posterior sample size
s.y.new=matrix(0,600,n)
s.addone=array(NA, dim = c(100,100,n))
thetasave=array(NA, dim = c(100,9,n))#theta for each county each quarter
#####addone
theta_year_num=array(NA, dim = c(100,4,n))#estimate using number for next year
theta_year_rate=array(NA, dim = c(100,4,n)) #estimate using rate(number/pop) for next year

p = dim(X)[2]
mu=mean(Bupdiff$Bupdiff[101:600])
sd=sd(Bupdiff$Bupdiff[101:600])
ii=2
nn=3

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
                            :(which(names(sampleused[[1]][1,])=='beta[1]')+(p-1))]
  
  
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
    xbeta<-rep(postbeta0,100)+postXmm[,1:p]%*%postbeta+rep(s2,1)
    thetasave[,(match(m,mm)+1),ll]<-exp(xbeta)*thetasave[,match(m,mm),ll]*Exp[,(match(m,mm)+1)]/Exp[,(match(m,mm))]
    
  }
  
  
  post.theta<-c(thetasave[,1,ll],thetasave[,2,ll],thetasave[,3,ll],
                thetasave[,4,ll],thetasave[,5,ll],thetasave[,6,ll])
  # y.new<-matrix(0,300,1)
  for (i in 1:600) {
    if(!is.na(post.theta[i])){
      s.y.new[i,ll]<-rpois(1,post.theta[i])
    } else
      c[i,ll]<-NA
    
   
    
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
  
  
}
low=matrix(NA, nrow = 100,ncol = 9)
high=matrix(NA, nrow = 100,ncol = 9)
mean=matrix(NA, nrow = 100,ncol = 9)

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


#bayesian p-value
ds.sd <- colSds(s.y.new[1:600,])
ds.mean <- colMeans(s.y.new[1:600,])

d0.mean <- mean(y[],na.rm=T)
d0.sd <- sd(y[],na.rm=T)

mean(d0.sd>ds.sd)
mean(d0.mean>ds.mean)



y=matrix(as.numeric(bup_patient$death),nrow=100,byrow=F)
estima2<-data.frame(time=c(2016:2024),Estimated=colSums(mean), 'Upper Bound'=colSums(high),
                    'Lower Bound'=colSums(low), Observed=c(colSums(y[,1:6],na.rm = T),rep(NA,3))
                    #Predicted=c(rep(NA,24),thetamean[1,])
                    )


newexp<-cbind(Exp,Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6],Exp[,6])


estima_r<-data.frame(time=c(2016:2024),Estimated=colMeans(mean/newexp[,1:9]), 'Upper Bound'=colMeans(high/newexp[,1:9]),
                     'Lower Bound'=colMeans(low/newexp[,1:9]), Observed=c(colMeans(y[,1:6]/Exp[,1:6],na.rm = T),rep(NA,3))
)



COL=c("Estimated" = "black",
      'Credible Bound' = "steelblue",
      "Observed" = "red")
basenum<-ggplot(data=estima_r,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  labs(color = '')+
  scale_color_manual(values = COL)+
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  labs(x = "Year", y="R")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())


baserate<-ggplot(data=estima_r,aes(x=time) )+
  geom_line(aes( y=Estimated, color ="Estimated" ),size=1)+
  geom_line(aes( y=Upper.Bound,color ='Credible Bound') ,linetype="twodash" )+
  geom_line(aes( y=Lower.Bound ,color ='Credible Bound'), linetype="twodash" )+
  geom_line(aes( y=Observed ,color ="Observed"))+
  labs(color = '')+
  scale_color_manual(values = COL)+
  geom_point(aes(x=time, y=Observed ),shape = 23,fill="red",size=2) +
  labs(x = "Year", y="Rate")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.8),   
        legend.title = element_blank(), 
        legend.text = element_text(size=16),
        axis.text.y = element_text( size=14),
        axis.text.x = element_text( size=14),
        axis.title = element_text(size=16))+
  theme(axis.title.x = element_blank())

# library(ggpubr)
# ggarrange(basenum,baserate,ncol=2,common.legend=T)

