
############################################################################
##Authors: Schoolmaster Jr, D.R., Zirbel, C.R & J.P. Cronin               ##                       
##Title: A graphical causal model for resolving species identity effects  ##
##       and biodiversity-ecosystem function correlations                 ##          
##Journal: Ecology                                                        ##
##                                                                        ##
##FTE Worked Example Code                                                 ##
############################################################################


# load some useful packages
library(xtable)
library(lmtest)
path<-"path/to/data/"
#read in the data
dat<-read.csv(paste0(path,"FTEData.csv"),header = TRUE)

#define useful functions
logit<-function(x)log(x/(1-x))

#soil quality from pca
soil.pca<-eigen(cor(dat[,c("pH_H2O","N_percent","Humus_g")]))
soil.quality<--as.matrix(scale(dat[,c("pH_H2O","N_percent","Humus_g")]))%*%soil.pca$vectors[,1]
cor(soil.quality,dat[,c("pH_H2O","N_percent","Humus_g")])

#grab basal area data
ba.dat<-dat[,grep("_ba",colnames(dat))]
#calculate presence/absence
pres.dat<-matrix(NA,nrow=dim(ba.dat)[1],ncol=4)
pres.dat<-ba.dat[,1:4]>0
colnames(pres.dat)<-c("Larix_pr","Fagus_pr","Picea_pr","Quercus_pr")
#create data.frame with all this info
ccover.dat<-data.frame(Canopy_cover=dat[,"Canopy_cover"],ba.dat,pres.dat)
#calculate variaton in density
fd.den<-apply(ba.dat[,1:4],1,var)
#calculate relative density based on basal area
rel.ba<-sweep(ba.dat[,-5],1,rowSums(ba.dat[,-5]),"/")

#trait values (see references in main text:
# F. sylvatica: Bartelink 1997; L. decidua: Fellner et al. 2016;
# P. abies: Wyka et al. 2007; Q. petraea: Davi et al. 2008)
sla<-c(Larix=117,Fagus=172,Picea=70.4,Quercus=153.84) #cm^2/g

#calcuate community weighed mean and variation of sla
cwm.sla<-as.matrix(rel.ba)%*%(as.matrix(sla))
fd.sla<-sqrt(1/(apply(rel.ba,1,function(x)sum(x^2)))*apply(as.matrix(rel.ba)*(t(apply(cwm.sla,1,function(x)(x-sla)^2))),1,sum))

#fit initial model
summary(f<-lm(logit(dat$Canopy_cover)~cwm.sla+soil.quality+(dat$dens_ba_site)))
print(xtable(f,digits = 3))       
#test each conditional indepdence statement (many ways to do this, I will use p-values of likelihood ratio tests)
summary(f.L1<-update(f,~.+dat$Larix_ba))
summary(f.F1<-update(f,~.+dat$Fagus_ba))
summary(f.P1<-update(f,~.+dat$Picea_ba))
summary(f.Q1<-update(f,~.+dat$Quercus_ba))
summary(f.Div1<-update(f,~.+dat$Tree_diversity))

#do likelihood ratio tests
p.L1<-lrtest(f,f.L1)
p.F1<-lrtest(f,f.F1)
p.P1<-lrtest(f,f.P1)
p.Q1<-lrtest(f,f.Q1)
p.Div1<-lrtest(f,f.Div1)
#grab the pvals from each
pvals<-c(L1=p.L1$`Pr(>Chisq)`[2],F1=p.F1$`Pr(>Chisq)`[2],P1=p.P1$`Pr(>Chisq)`[2],Q1=p.Q1$`Pr(>Chisq)`[2],Div1=p.Div1$`Pr(>Chisq)`[2])
ord<-order(pvals)
print(pvals[ord])
#do Holm-Bonferroni adjustment with alpha=0.05 (or whatever multiple testing procedure, if any, you prefer)
#this one makes the hypotheses of no species/diversity effects harder to reject (used by Textor 2018)

adj.alpha<-sapply(1:length(pvals),function(x)0.05/(length(pvals)-x+1))
#reject hypo that these are equal to zero
which(pvals[ord]<=adj.alpha)

#this says that we have at least on failed test associated with the effect of Fagus

#alternatively we could use Fisher's combined test (i.e Shipley 2000) to test that  
# at least one of the 'null hypotheses' of theta_i=0 is false.
#this test makes the set of hypotheses of no species/diversity effects harder to reject

1-pchisq(q = sum(-2*log(pvals)),df=2*length(pvals))
#this says that at least one test has failed, but doesn't directly indicate which one. To achieve that we
#use closure methods,i.e. if an effect is sig in every intersection of the others it is signficant

#create function to do closure
closeit<-function(y)sapply(y,function(x)1-pchisq(q = sum(-2*log(pvals[x])),df=2*length(x)))
#get all combination of 4,3,2,1 to build the intersection heirarchy
H4<-combn(1:5,4,simplify = F)
H3<-combn(1:5,3,simplify = F)
H2<-combn(1:5,2,simplify = F)
H1<-combn(1:5,1,simplify = F)

heir<-c('H4','H3','H2','H1')
#go though each effect, stop if you find one intersection that is not <0.05
ans<-NULL
i=1
for(i in 1:5){
  ans[i]<-FALSE
  j=1
  while((!ans[i])&&j<=4){
    ans[i]<-any(closeit(get(heir[j][[1]]))[sapply(get(heir[j][[1]]),function(x)i%in%x)]>0.05)
    j=j+1}
}
#see which specific effects are signifcant under closure of fisher 
pvals[!ans]


#update model to include community weighted variation in sla
summary(f2<-lm(logit(dat$Canopy_cover)~cwm.sla+fd.sla+soil.quality+(dat$dens_ba_site)))
print(xtable(f2,digits=3))
#test each conditional independence statement 
summary(f.L2<-update(f2,~.+dat$Larix_ba))
summary(f.F2<-update(f2,~.+dat$Fagus_ba))
summary(f.P2<-update(f2,~.+dat$Picea_ba))
summary(f.Q2<-update(f2,~.+dat$Quercus_ba))
summary(f.Div2<-update(f2,~.+dat$Tree_diversity))

#do likelihood ratio tests
p.L2<-lrtest(f2,f.L2)
p.F2<-lrtest(f2,f.F2)
p.P2<-lrtest(f2,f.P2)
p.Q2<-lrtest(f2,f.Q2)
p.Div2<-lrtest(f2,f.Div2)

pvals2<-c(L2=p.L2$`Pr(>Chisq)`[2],F2=p.F2$`Pr(>Chisq)`[2],P2=p.P2$`Pr(>Chisq)`[2],Q2=p.Q2$`Pr(>Chisq)`[2],Div2=p.Div2$`Pr(>Chisq)`[2])
#do Holm-Bonferroni adjustment with alpha=0.05 (or whatever multiple testing procedure, if any, you prefer)
ord<-order(pvals2)
adj.alpha2<-sapply(1:length(pvals2),function(x)0.05/(length(pvals2)-x+1))
#reject hypo that these are equal to zero
which(pvals2[ord]<=adj.alpha2)

#this says that we have failed to reject any of the HO of no species effects 

#alternatively we could use Fisher's combined test (i.e Shipley 2000) to test that  
# at least one of the 'null hypotheses' of theta_i=0 is false

1-pchisq(q = sum(-2*log(pvals2)),df=2*length(pvals2))
#this test suggests that at least one of the tested effects is non-zero, suggesting 
#that we might think about which other traits (and how they combine to affect) to additional 
#affect canopy cover.

#Use closure of fisher's C to find which are significant (Henning & Westfall 2019)
#create function to do closure
closeit<-function(y)sapply(y,function(x)1-pchisq(q = sum(-2*log(pvals2[x])),df=2*length(x)))
#get all combination of 4,3,2,1 to build the intersection heirarchy
H4<-combn(1:5,4,simplify = F)
H3<-combn(1:5,3,simplify = F)
H2<-combn(1:5,2,simplify = F)
H1<-combn(1:5,1,simplify = F)

heir<-c('H4','H3','H2','H1')
#go though each effect, stop if you find one intersection that is not <0.05
ans<-NULL
i=1
for(i in 1:5){
  ans[i]<-FALSE
  j=1
  while((!ans[i])&&j<=4){
    ans[i]<-any(closeit(get(heir[j][[1]]))[sapply(get(heir[j][[1]]),function(x)i%in%x)]>0.05)
    j=j+1}
}
#see which specific effects are signifcant under closure of fisher 
pvals2[!ans]
# we get the strange situation where the global test fails while no single effect is
#judged to be signifiant

#In this case, I judge this
#good enough to accept, although there is a hint that another trait 
#associated with Fagus might be contributing to total canopy cover
