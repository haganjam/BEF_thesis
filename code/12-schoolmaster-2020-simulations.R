############################################################################
##Authors: Schoolmaster Jr, D.R., Zirbel, C.R & J.P. Cronin               ##                       
##Title: A graphical causal model for resolving species identity effects  ##
##       and biodiversity-ecosystem function correlations                 ##          
##Journal: Ecology                                                        ##
##                                                                        ##
##Simulation methods and code                                             ##
############################################################################

####Overview############################################################
# To test the conclusion of the causal analysis of the DAG in Fig. 3 that 
# species effects and biodiversity correlations result from statistically
# misspecified models, we generated data from simulations of communities.
# We simulated data for a set of 5 species, each with a unique response to
# an environmental gradient plus an additional tendency to covary with one 
# another independent of the environmental responses. Each species was
# assigned 4 functional traits. The traits also have non-zero covariance 
# with one another. Both the abundance and trait correlation distributions
# were sampled from a Beta distribution.
# After the abundance generating distribution and trait generating 
# distributions were initialized, communities of the 5 species are drawn 
#  from the abundance distribution, one each of 100 different levels of the 
# environmental gradient, uniformly distributed. Using these data, two 
# functional trait composites were calculated : (1) a community-weighted mean
# and (2) a community-weighted variance of one of the traits. A value for the
# site's ecosystem function is calculated from these trait composites. 
# The ecosystem function was calculated as a linear combination of the two 
# functional trait composites, the environmental variable and some Gaussian
# distributed error.
# Next, we fit strategically misspecified models to the simulated data. We fit 
# models that left out either the environmental term, one of the functional trait
# composites, or we used a functional trait composite based on the wrong functional
# trait. After fitting the misspecified model, we updated it to include, one at 
# a time, the abundance of each species or a measure of diversity. We recorded the
# t-value associated with the added variable. The whole process, was repeated 1000 
# times for each type of misspecification. As a check, we also recorded the t-value
# the resulted from updating the correctly specified model with each variable. The 
# resulting distributions of t-values are shown in Fig. 5 of the main text. The code 
# for conducting the simulations in included below.
#################################################################################

#grab a function capable of giving multivarite normal random numbers
library(mvtnorm)
#set environmental gradient
E<-rnorm(100,0,5)
range(E)
#set linear environmental response with tradeoff between intercept and slope
#the trait x sets up the envionmental response of the species. Can adjust th
#numbers to gives differnt intercepts and slopes
mu<-function(x)10-0.05*x+E*x/2

#plot example mean responses to ennironmenal gradient
plot(E,mu(8),type='l',ylim=c(0,70),xlab="Environment",ylab="Species Response")
lines(E,mu(4),type='l')
lines(E,mu(0),type='l')
lines(E,mu(-4),type='l')
lines(E,mu(-8),type='l')


####misspecification 1: neglect effect of E####
noE<-matrix(0,nrow=1000,ncol=6)
yesE<-matrix(0,nrow=1000,ncol=6)
for(j in 1:1000){
  
  #select E-independent speices abundance correlation matrix by drawing pairwise
  #covariance from a Beta distribution
  mat<-matrix(0,5,5)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(10,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  abund.cormat<-mat
  
  #select standard devition of species abundance 
  sp.sd<-10
  sig<-c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd)%*%t(c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd))*abund.cormat
  #select trait x, whiich defines the environmental response from a unifom distribution U(-8,8)
  x<-16*runif(5)-8
  #drawn species abundances for 100 sites form multivariate normal dist using 
  spp<-cbind(10*mu(x[1]),10*mu(x[2]),10*mu(x[3]),10*mu(x[4]),10*mu(x[5]))+rmvnorm(100,mean=c(0,0,0,0,0),sigma = sig)
  #set all negative abundances to zero
  spp[spp<0]<-NA
  colnames(spp)<-paste0("sp",1:5)
  #calculate Shannon diversity at each 'site'
  biodiv<-apply(spp,1,function(x)-sum(x/sum(x,na.rm = T)*log(x/sum(x,na.rm = T)),na.rm = T))
  
  #simulate trait matrix
  #use same method as above to set up trait correlations 
  mat<-matrix(0,4,4)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(6,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  trait.cormat<-mat
  
  sig<-c(10,5,2,1)%*%t(c(10,5,2,1))*trait.cormat
  traits<-rmvnorm(5,mean=c(0,0,0,0),sigma = sig)
  rownames(traits)<-c("sp1","sp2","sp3","sp4","sp5")
  colnames(traits)<-c("T1","T2","T3","T4")
  #calcualte relative abundance of species within sites
  relat<-sweep(spp,1,rowSums(spp,na.rm = T),"/")
  relat[which(is.na(relat))]<-0
  #calculate commuinty-weighted means and variances to be used to drive ecosystem functions
  CWM<-as.matrix(relat)%*%as.matrix(traits)
  FD<-cbind(T1=apply(relat,1,function(x)sum(x^2))*var(traits[,1]),
            T2=apply(relat,1,function(x)sum(x^2))*var(traits[,2]),
            T3=apply(relat,1,function(x)sum(x^2))*var(traits[,3]),
            T4=apply(relat,1,function(x)sum(x^2))*var(traits[,4]))
  
  #define ecosystem function as driven by environmental response, E, CWM of trait 1, and FD of trait 1
  #plus some error
  EF<-cbind(CWM[,1],FD[,1],E)%*%c(5,1,1)+rnorm(100,0,5)
  
  #fit initial model (with neglects E)
  f<-lm(EF~CWM[,1]+FD[,1])
  #fit correctly specified model
  f.yes<-lm(EF~CWM[,1]+FD[,1]+E)
  #grab the t-value for the testable implications that spp abund and biodiv are independent of EF given 
  #corrent function trait composities (CWM and FD) and enviornment
  for(k in 1:6)noE[j,k]<-summary(update(f,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][4, "t value"]
  for(k in 1:6)yesE[j,k]<-summary(update(f.yes,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][5, "t value"]
}

#how many of the simulations resulted in at least one "species identity effect" or "diveristy effect"
sum(apply(noE,1,function(x)any(abs(x)>2)))/1000

#check the same thing for the correctly specified model, here the answer should be around
# 0.05x6 = 0.30
sum(apply(yesE,1,function(x)any(abs(x)>2)))/1000

#plot the results
plot(1:6,apply(noE,2,function(x)mean(x)),ylim=c(-7,7),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(noE,2,function(x)mean(x)),1:6,apply(noE,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(noE,2,function(x)mean(x)),1:6,apply(noE,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Diversity"))
points(1:6,apply(noE,2,function(x)mean(x)),pch=21,bg="grey")

###Plot the results for the correctly specified model. Here the 95% CI should run from -2 to 2
plot(1:6,apply(yesE,2,function(x)mean(x)),ylim=c(-7,7),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(yesE,2,function(x)mean(x)),1:6,apply(yesE,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(yesE,2,function(x)mean(x)),1:6,apply(yesE,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Diversity"))
points(1:6,apply(yesE,2,function(x)mean(x)),pch=21,bg="grey")


###the code below follows the same logic, but examines the effect of 
###misspecification by leaving out CWN and FD or including a Community-weighted mean based
###the incorrect functional trait.

####misspecification 2: neglect effect of CWM####
noCWM<-matrix(0,nrow=1000,ncol=6)
yesCWM<-matrix(0,nrow=1000,ncol=6)
for(j in 1:1000){
  mat<-matrix(0,5,5)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(10,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  abund.cormat<-mat
  
  sp.sd<-200
  sig<-c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd)%*%t(c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd))*abund.cormat
  x<-16*runif(5)-8
  spp<-cbind(10*mu(x[1]),10*mu(x[2]),10*mu(x[3]),10*mu(x[4]),10*mu(x[5]))+rmvnorm(100,mean=c(0,0,0,0,0),sigma = sig)
  
  spp[spp<0]<-NA
  colnames(spp)<-paste0("sp",1:5)
  biodiv<-apply(spp,1,function(x)-sum(x/sum(x,na.rm = T)*log(x/sum(x,na.rm = T)),na.rm = T))
  
  #simulate trait matrix
  mat<-matrix(0,4,4)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(6,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  trait.cormat<-mat
  sig<-c(10,5,2,1)%*%t(c(10,5,2,1))*trait.cormat
  traits<-rmvnorm(5,mean=c(0,0,0,0),sigma = sig)
  rownames(traits)<-c("sp1","sp2","sp3","sp4","sp5")
  colnames(traits)<-c("T1","T2","T3","T4")
  relat<-sweep(spp,1,rowSums(spp,na.rm = T),"/")
  relat[which(is.na(relat))]<-0
  CWM<-as.matrix(relat)%*%as.matrix(traits)
  
  FD<-cbind(T1=apply(relat,1,function(x)sum(x^2))*var(traits[,1]),
            T2=apply(relat,1,function(x)sum(x^2))*var(traits[,2]),
            T3=apply(relat,1,function(x)sum(x^2))*var(traits[,3]),
            T4=apply(relat,1,function(x)sum(x^2))*var(traits[,4]))
  
  EF<-cbind(CWM[,1],FD[,1],E)%*%c(5,1,1)+rnorm(100,0,10)
  f<-lm(EF~FD[,1]+E)
  f.yes<-lm(EF~FD[,1]+CWM[,1]+E)
  for(k in 1:6)noCWM[j,k]<-summary(update(f,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][4, "t value"]
  for(k in 1:6)yesCWM[j,k]<-summary(update(f.yes,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][5, "t value"]
  
}

plot(1:6,apply(noCWM,2,function(x)mean(x)),ylim=c(-7,7),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(noCWM,2,function(x)mean(x)),1:6,apply(noCWM,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(noCWM,2,function(x)mean(x)),1:6,apply(noCWM,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Diversity"))
points(1:6,apply(noCWM,2,function(x)mean(x)),pch=21,bg="grey")
##correct model##
plot(1:6,apply(yesCWM,2,function(x)mean(x)),ylim=c(-7,7),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(yesCWM,2,function(x)mean(x)),1:6,apply(yesCWM,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(yesCWM,2,function(x)mean(x)),1:6,apply(yesCWM,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Diversity"))
points(1:6,apply(yesCWM,2,function(x)mean(x)),pch=21,bg="grey")


#####misspecification 3: neglect effect of trait dispersion (FD)####
noFD<-matrix(0,nrow=1000,ncol=6)
yesFD<-matrix(0,nrow=1000,ncol=6)
for(j in 1:1000){
  mat<-matrix(0,5,5)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(10,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  abund.cormat<-mat
  sp.sd<-200
  sig<-c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd)%*%t(c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd))*abund.cormat
  x<-16*runif(5)-8
  spp<-cbind(10*mu(x[1]),10*mu(x[2]),10*mu(x[3]),10*mu(x[4]),10*mu(x[5]))+rmvnorm(100,mean=c(0,0,0,0,0),sigma = sig)
  
  spp[spp<0]<-NA
  colnames(spp)<-paste0("sp",1:5)
  biodiv<-apply(spp,1,function(x)-sum(x/sum(x,na.rm = T)*log(x/sum(x,na.rm = T)),na.rm = T))
  
  #simulate trait matrix
  mat<-matrix(0,4,4)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(6,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  trait.cormat<-mat
  sig<-c(10,5,2,1)%*%t(c(10,5,2,1))*trait.cormat
  traits<-rmvnorm(5,mean=c(0,0,0,0),sigma = sig)
  rownames(traits)<-c("sp1","sp2","sp3","sp4","sp5")
  colnames(traits)<-c("T1","T2","T3","T4")
  relat<-sweep(spp,1,rowSums(spp,na.rm = T),"/")
  relat[which(is.na(relat))]<-0
  sum(as.matrix(relat)[1,]*as.matrix(traits)[,1])
  CWM<-as.matrix(relat)%*%as.matrix(traits)
  
  FD<-cbind(T1=apply(relat,1,function(x)sum(x^2))*var(traits[,1]),
            T2=apply(relat,1,function(x)sum(x^2))*var(traits[,2]),
            T3=apply(relat,1,function(x)sum(x^2))*var(traits[,3]),
            T4=apply(relat,1,function(x)sum(x^2))*var(traits[,4]))
  
  
  EF<-cbind(CWM[,1],FD[,1],E)%*%c(5,1,1)+rnorm(100,0,10)
  f<-lm(EF~CWM[,1]+E)
  f.yes<-lm(EF~CWM[,1]+FD[,1]+E)
  for(k in 1:6)noFD[j,k]<-summary(update(f,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][4, "t value"]
  for(k in 1:6)yesFD[j,k]<-summary(update(f.yes,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][5, "t value"]
}

sum(apply(noFD[,1:5],1,function(x)any(abs(x)>2)))/1000
sum(apply(noE[,1:5],1,function(x)any(abs(x)>2)))/1000
sum(apply(yesE[,1:5],1,function(x)any(abs(x)>2)))/1000
sum(apply(noCWM[,1:5],1,function(x)any(abs(x)>2)))/1000

apply(noFD,2,function(x)sum(abs(x)>2))/1000
apply(noE,2,function(x)sum(abs(x)>2))/1000
apply(noCWM,2,function(x)sum(abs(x)>2))/1000

apply(yesFD,2,function(x)sum(abs(x)>2))/1000
apply(yesE,2,function(x)sum(abs(x)>2))/1000
apply(yesCWM,2,function(x)sum(abs(x)>2))/1000

plot(0:7,c(NA,apply(noFD,2,function(x)mean(x)),NA),ylim=c(-25,25),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(noFD,2,function(x)mean(x)),1:6,apply(noFD,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(noFD,2,function(x)mean(x)),1:6,apply(noFD,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Div."))
points(1:6,apply(noFD,2,function(x)mean(x)),pch=21,bg="grey")
abline(h=c(2,-2))
###correct model
plot(1:6,apply(yesFD,2,function(x)mean(x)),ylim=c(-7,7),xaxt='n',ylab="Residual t-value",xlab="Effect")
arrows(1:6,apply(yesFD,2,function(x)mean(x)),1:6,apply(yesFD,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(yesFD,2,function(x)mean(x)),1:6,apply(yesFD,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
axis(side=1,1:6,c(paste0('sp',1:5),"Div."))
points(1:6,apply(yesFD,2,function(x)mean(x)),pch=21,bg="grey")


#### misspecification 4: Use a CWM of incorrect functional trait #####
Wmet<-matrix(0,nrow=1000,ncol=6)
Cmet<-matrix(0,nrow=1000,ncol=6)
for(j in 1:1000){
  mat<-matrix(0,5,5)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(10,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  abund.cormat<-mat
  sp.sd<-200
  sig<-c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd)%*%t(c(sp.sd,sp.sd,sp.sd,sp.sd,sp.sd))*abund.cormat
  x<-16*runif(5)-8
  spp<-cbind(10*mu(x[1]),10*mu(x[2]),10*mu(x[3]),10*mu(x[4]),10*mu(x[5]))+rmvnorm(100,mean=c(0,0,0,0,0),sigma = sig)
  
  spp[spp<0]<-0
  colnames(spp)<-paste0("sp",1:5)
  biodiv<-apply(spp,1,function(x)-sum(x/sum(x,na.rm = T)*log(x/sum(x,na.rm = T)),na.rm = T))
  
  #simulate trait matrix
  mat<-matrix(0,4,4)
  mat[upper.tri(mat,diag = F)]<-2*rbeta(6,15,15)-1
  mat<-mat+t(mat)
  diag(mat)<-1
  trait.cormat<-mat
  sig<-c(10,5,2,1)%*%t(c(10,5,2,1))*trait.cormat
  traits<-rmvnorm(5,mean=c(0,0,0,0),sigma = sig)
  rownames(traits)<-c("sp1","sp2","sp3","sp4","sp5")
  colnames(traits)<-c("T1","T2","T3","T4")
  relat<-sweep(spp,1,rowSums(spp,na.rm = T),"/")
  relat[which(is.na(relat))]<-0
  
  CWM<-as.matrix(relat)%*%as.matrix(traits)
  CUM<-as.matrix(spp)%*%as.matrix(traits)
  FD<-cbind(T1=apply(relat,1,function(x)sum(x^2))*var(traits[,1]),
            T2=apply(relat,1,function(x)sum(x^2))*var(traits[,2]),
            T3=apply(relat,1,function(x)sum(x^2))*var(traits[,3]),
            T4=apply(relat,1,function(x)sum(x^2))*var(traits[,4]))
  
  
  EF<-cbind(CWM[,1],FD[,1],E)%*%c(5,1,1)+rnorm(100,0,10)
  #this misspecification grabs a random CWM based on the wrong trait
  f<-lm(EF~CWM[,sample(2:4)[1]]+FD[,1]+E)
  f.yes<-lm(EF~CWM[,1]+FD[,1]+E)
  for(k in 1:6)Wmet[j,k]<-summary(update(f,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][4, "t value"]
  for(k in 1:6)Cmet[j,k]<-summary(update(f.yes,~.+cbind(spp,biodiv)[,k]))[["coefficients"]][5, "t value"]
}

plot(0:7,c(NA,apply(Wmet,2,function(x)mean(x)),NA),ylim=c(-8,8),xaxt='n',ylab="",xlab="Effect",cex.lab=1.25,cex.axis=1.25)
arrows(1:6,apply(Wmet,2,function(x)mean(x)),1:6,apply(Wmet,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(Wmet,2,function(x)mean(x)),1:6,apply(Wmet,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
abline(h=c(-2,2),lty=3)
axis(side=1,1:6,c(paste0('Sp',1:5),"Div."))
points(1:6,apply(Wmet,2,function(x)mean(x)),pch=21,bg="grey")
#correct model
plot(0:7,c(NA,apply(Cmet,2,function(x)mean(x)),NA),ylim=c(-8,8),xaxt='n',ylab="",xlab="Effect",cex.lab=1.25,cex.axis=1.25)
arrows(1:6,apply(Cmet,2,function(x)mean(x)),1:6,apply(Cmet,2,function(x)quantile(x,probs = 0.025)),angle=90,length=.05)
arrows(1:6,apply(Cmet,2,function(x)mean(x)),1:6,apply(Cmet,2,function(x)quantile(x,probs = 0.975)),angle=90,length=.05)
abline(h=0,lty=2)
abline(h=c(-2,2),lty=3)
axis(side=1,1:6,c(paste0('Sp',1:5),"Div."))
points(1:6,apply(Cmet,2,function(x)mean(x)),pch=21,bg="grey")


