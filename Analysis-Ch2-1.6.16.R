rm(list=ls())
library(spatstat)
data<-read.csv("PowdermillData-4.6.15.csv", header=TRUE)
str(data)
source("BagchiFunctions-ReplicatedPPP.R")

##make counter #1
##do it all by size class

hlist=c(1:4)
ByHeight<-lapply(hlist, function(x) data[data$HeightClass==x,])

ByHeight
list2.0=vector(mode="list")
list2=vector(mode="list")

for (i in seq_along(hlist)){
  num<-unique(ByHeight[[i]]$Species)
  datacount<-length(num)
  datacounter<-c(1:datacount)
  list2.0[[i]]=vector(mode="list", length=datacount)
  list2[[i]]=vector(mode="list", length=datacount)
  list2[[i]]<-lapply(num, function(x) ByHeight[[i]][ByHeight[[i]]$Species==x,])
  print(num)
}

list2
list2.0
num<-unique(ByHeight[[1]]$Species)

list2
list3.0=vector(mode="list")
list3=vector(mode="list")

for (i in seq_along(hlist)){
  
  
  datacount=length(list2[[i]])
  print(datacount)
  datacounter=c(1:datacount)
  #list3.[[i]]=vector(mode="list", length=datacount)
  list3[[i]]=vector(mode="list", length=datacount)
  
  for (j in seq_along(datacounter)){
    num2=unique(list2[[i]][[j]]$PlotNum)
    datacount2<-length(num2)
    print(datacount2)
    #list3.0[[i]][[j]]=vector(mode="list", length=datacount2)
    list3[[i]][[j]]=vector(mode="list", length=datacount2)
    list3[[i]][[j]]<-lapply(num2, function(x) list2[[i]][[j]][list2[[i]][[j]]$PlotNum==x,])
    #list3[[i]][[j]]<-lapply(list3.0[[i]][[j]], na.omit)
  }
}

list3
commX=vector(mode="list")
commY=vector(mode="list")
commXrange=vector(mode="list")
commYrange=vector(mode="list")
commspatdata=vector(mode="list")
commwin=vector(mode="list")
commyin=vector(mode="list")
commreal=vector(mode="list")

for (i in seq_along(hlist)){
  
  
  datacount=length(list3[[i]])
  print(datacount)
  datacounter=c(1:datacount)
  commX[[i]]=vector(mode="list", length=datacount)
  commY[[i]]=vector(mode="list", length=datacount)
  commXrange[[i]]=vector(mode="list", length=datacount)
  commYrange[[i]]=vector(mode="list", length=datacount)
  commspatdata[[i]]=vector(mode="list", length=datacount)
  commwin[[i]]=vector(mode="list", length=datacount)
  commyin[[i]]=vector(mode="list", length=datacount)
  commreal[[i]]=vector(mode="list", length=datacount)
  
  for (j in seq_along(datacounter)){
    
    datacount2<-length(list3[[i]][[j]])
    print(datacount2)
    datacounter2=c(1:datacount2)
    commX[[i]][[j]]=vector(mode="list", length=datacount2)
    commY[[i]][[j]]=vector(mode="list", length=datacount2)
    commXrange[[i]][[j]]=vector(mode="list", length=datacount2)
    commYrange[[i]][[j]]=vector(mode="list", length=datacount2)
    commspatdata[[i]][[j]]=vector(mode="list", length=datacount2)
    commwin[[i]][[j]]=vector(mode="list", length=datacount2)
    commyin[[i]][[j]]=vector(mode="list", length=datacount2)
    commreal[[i]][[j]]=vector(mode="list", length=datacount2)
    
    for (k in seq_along(datacounter2)){
      
      
      commX[[i]][[j]][[k]]<-list3[[i]][[j]][[k]]$Xcoord
      commY[[i]][[j]][[k]]<-list3[[i]][[j]][[k]]$Ycoord
      commXrange[[i]][[j]][[k]]=c(min(commX[[i]][[j]][[k]])-1, max(commX[[i]][[j]][[k]])+1)
      commYrange[[i]][[j]][[k]]=c(min(commY[[i]][[j]][[k]])-1, max(commY[[i]][[j]][[k]])+1)
      commspatdata[[i]][[j]][[k]]<-ppp(commX[[i]][[j]][[k]],commY[[i]][[j]][[k]],commXrange[[i]][[j]][[k]],commYrange[[i]][[j]][[k]])
      commwin[[i]][[j]][[k]]<-convexhull(commspatdata[[i]][[j]][[k]])
      commreal[[i]][[j]][[k]]<-ppp(commX[[i]][[j]][[k]],commY[[i]][[j]][[k]],window=commwin[[i]][[j]][[k]])
    }
  }
}



##now make into separate list for understory plants for height classes 


listH1<-matrix()
listH1.2<-matrix()

datacount<-length(list3[[1]])
datacounter=c(1:datacount)
print(datacount)


for (j in seq_along(datacounter)){
  listH1.2<-cbind(commreal[[1]][[j]])
  listH1<-rbind(listH1, listH1.2)      
} 

listH1<-listH1[-1,]

listH2<-matrix()
listH2.2<-matrix()

datacount<-length(list3[[2]])
datacounter=c(1:datacount)
print(datacount)


for (j in seq_along(datacounter)){
  listH2.2<-cbind(commreal[[2]][[j]])
  listH2<-rbind(listH2, listH2.2)      
} 

listH2<-listH2[-1,]

listH3<-matrix()
listH3.2<-matrix()

datacount<-length(list3[[3]])
datacounter=c(1:datacount)
print(datacount)


for (j in seq_along(datacounter)){
  listH3.2<-cbind(commreal[[3]][[j]])
  listH3<-rbind(listH3, listH3.2)      
} 

listH3<-listH3[-1,]

listH4<-matrix()
listH4.2<-matrix()

datacount<-length(list3[[4]])
datacounter=c(1:datacount)
print(datacount)


for (j in seq_along(datacounter)){
  listH4.2<-cbind(commreal[[4]][[j]])
  listH4<-rbind(listH4, listH4.2)      
} 

listH4<-listH4[-1,]


##Get Ready to Graph

xrange=seq(from=0, to=5, by=.1)
require(abind)


listH1
KfuncsH1 <- lapply(listH1, Kest, r=xrange, correction=NULL, ratio=FALSE)
## put in matrix - distances are columns, reps are rows
KfuncsMatH1 <- do.call('rbind', lapply(KfuncsH1, function(x) x$trans))

#remove NA's
KfuncsMatH1
KfuncsMatH1<-KfuncsMatH1[-c(4,7,11,15,16,17,18,19,21,23,24,25,26,27,31,35,37,38,39,40),]
KfuncsMatH1
KfuncsMatH1<-KfuncsMatH1[,-c(17:61)]
KfuncsMatH1<-KfuncsMatH1[-c(19,8),]

## calculate weights
covwtsH1 <-  do.call('rbind', lapply(listH1, kfunc.weights.calc, type='nx',
                                          r=xrange, corr=NULL))

##remove same parts as above
covwtsH1 <- apply(covwtsH1, 2, function(x) {
  x <- x/mean(x)
  return(x)})
covwtsH1<-covwtsH1[-c(4,7,11,15,16,17,18,19,21,23,24,25,26,27,31,35,37,38,39,40),]
covwtsH1<-covwtsH1[,-c(17:61)]
covwtsH1<-covwtsH1[-c(19,8),]

interceptdatH1 <- data.frame(int=rep(1, nrow(KfuncsMatH1)))

## fit the models
lmmodsH1 <- try(lapply(1:ncol(KfuncsMatH1), function(j)
  kfunclm(k=KfuncsMatH1[,j], weights=covwtsH1[,j],
          dat=interceptdatH1, form='1')), silent=TRUE)
## do the bootstrapping for confidence intervals
lc <- matrix(1, ncol=1)
bsciH1 <- lm.t.boot(lmmodsH1, lincomb=lc, nsim=99, alpha=0.05,
                         simple.method=FALSE)

## pull out the parts needed for plotting and add 0s for distances
## where all reps were 0
lmKpredH1 <-bsciH1$lmKpred
lower.CIH1<-bsciH1$lower
upper.CIH1<-bsciH1$upper
## put 0s back in manually where there were only 0s in K -
## models often fail and give us the same result

resultsH1<-list(lmK=lmmodsH1, lmKpred=lmKpredH1,lower=lower.CIH1, upper=upper.CIH1, sps=listH1)

plotdat.kfuncsH1<-do.call('rbind',resultsH1$sps)
plotdat.kfuncsH1 <-do.call('rbind', lapply(resultsH1$sps, function(x){
  ki <- Kest(x, r=seq(from=0,to=6, by=.1), correction='iso')
  L <- K2L(ki$iso, ki$r)
  dat <- data.frame(distance=ki$r, K=ki$iso, L=L)
  return(dat)
}))

xrange2=seq(from=.1, to=5, by=.3)
pdat3<-data.frame(pred=do.call('rbind', resultsH1$lmKpred),
                  lci=do.call('rbind', resultsH1$lower),
                  uci=do.call('rbind', resultsH1$upper))
pdat3$distance<-xrange
pdat3$Ltheo <- Kpoisson(pdat3$distance, L=TRUE)

library(ggplot2)
none=element_blank() 
pl1<-ggplot(data=plotdat.kfuncsH1)+
  geom_smooth(data=plotdat.kfuncsH1, aes(x=distance, y=L),colour="black", linetype=1)+
  theme(panel.grid.major = none, panel.grid.minor = none)+
  theme(panel.background = none)+
  theme(plot.background=none)+
  geom_abline(intercept=0,slope=0, colour="black", linetype=2)
  
pl1

listH2
KfuncsH2 <- lapply(listH2, Kest, r=xrange, correction=NULL, ratio=FALSE)
## put in matrix - distances are columns, reps are rows
KfuncsMatH2 <- do.call('rbind', lapply(KfuncsH2, function(x) x$trans))

#remove NA's
KfuncsMatH2
KfuncsMatH2<-KfuncsMatH2[-c(7:11,13,15,18:20,24,29,32,35,37,40),]
KfuncsMatH2
KfuncsMatH2<-KfuncsMatH2[,-c(12:61)]
KfuncsMatH2
KfuncsMatH2<-KfuncsMatH2[-c(18,5),]

## calculate weights
covwtsH2 <-  do.call('rbind', lapply(listH2, kfunc.weights.calc, type='nx',
                                     r=xrange, corr=NULL))

##remove same parts as above
covwtsH2 <- apply(covwtsH2, 2, function(x) {
  x <- x/mean(x)
  return(x)})
covwtsH2<-covwtsH2[-c(7:11,13,15,18:20,24,29,32,35,37,40),]
covwtsH2<-covwtsH2[,-c(12:61)]
covwtsH2<-covwtsH2[-c(18,5),]


interceptdatH2 <- data.frame(int=rep(1, nrow(KfuncsMatH2)))

## fit the models
lmmodsH2 <- try(lapply(1:ncol(KfuncsMatH2), function(j)
  kfunclm(k=KfuncsMatH2[,j], weights=covwtsH2[,j],
          dat=interceptdatH2, form='1')), silent=TRUE)
## do the bootstrapping for confidence intervals
lc <- matrix(1, ncol=1)
bsciH2 <- lm.t.boot(lmmodsH2, lincomb=lc, nsim=99, alpha=0.05,
                    simple.method=FALSE)

## pull out the parts needed for plotting and add 0s for distances
## where all reps were 0
lmKpredH2 <-bsciH2$lmKpred
lower.CIH2<-bsciH2$lower
upper.CIH2<-bsciH2$upper
## put 0s back in manually where there were only 0s in K -
## models often fail and give us the same result

resultsH2<-list(lmK=lmmodsH2, lmKpred=lmKpredH2,lower=lower.CIH2, upper=upper.CIH2, sps=listH2)

plotdat.kfuncsH2<-do.call('rbind',resultsH2$sps)
plotdat.kfuncsH2 <-do.call('rbind', lapply(resultsH2$sps, function(x){
  ki <- Kest(x, r=seq(from=0,to=6, by=.1), correction='iso')
  L <- K2L(ki$iso, ki$r)
  dat <- data.frame(distance=ki$r, K=ki$iso, L=L)
  return(dat)
}))

pl1<-pl1+
  geom_smooth(data=plotdat.kfuncsH2, aes(x=distance, y=L),colour="black", linetype=3)
pl1



none=element_blank()

listH3
KfuncsH3 <- lapply(listH3, Kest, r=xrange, correction=NULL, ratio=FALSE)
## put in matrix - distances are columns, reps are rows
KfuncsMatH3 <- do.call('rbind', lapply(KfuncsH3, function(x) x$trans))

#remove NA's
KfuncsMatH3
KfuncsMatH3<-KfuncsMatH3[-c(3,5:7, 9:11, 13, 15:18, 22:24, 26, 31:33, 36, 43),]
KfuncsMatH3
KfuncsMatH3<-KfuncsMatH3[,-c(19:61)]
KfuncsMatH3
KfuncsMatH3<-KfuncsMatH3[-c(12,16),]

## calculate weights
covwtsH3 <-  do.call('rbind', lapply(listH3, kfunc.weights.calc, type='nx',
                                     r=xrange, corr=NULL))

##remove same parts as above
covwtsH3 <- apply(covwtsH3, 2, function(x) {
  x <- x/mean(x)
  return(x)})

covwtsH3<-covwtsH3[-c(3,5:7, 9:11, 13, 15:18, 22:24, 26, 31:33, 36, 43),]
covwtsH3<-covwtsH3[,-c(19:61)]
covwtsH3<-covwtsH3[-c(12,16),]


interceptdatH3 <- data.frame(int=rep(1, nrow(KfuncsMatH3)))

## fit the models
lmmodsH3 <- try(lapply(1:ncol(KfuncsMatH3), function(j)
  kfunclm(k=KfuncsMatH3[,j], weights=covwtsH3[,j],
          dat=interceptdatH3, form='1')), silent=TRUE)
## do the bootstrapping for confidence intervals
lc <- matrix(1, ncol=1)
bsciH3 <- lm.t.boot(lmmodsH3, lincomb=lc, nsim=99, alpha=0.05,
                    simple.method=FALSE)

## pull out the parts needed for plotting and add 0s for distances
## where all reps were 0
lmKpredH3 <-bsciH3$lmKpred
lower.CIH3<-bsciH3$lower
upper.CIH3<-bsciH3$upper
## put 0s back in manually where there were only 0s in K -
## models often fail and give us the same result

resultsH3<-list(lmK=lmmodsH3, lmKpred=lmKpredH3,lower=lower.CIH3, upper=upper.CIH3, sps=listH3)

plotdat.kfuncsH3<-do.call('rbind',resultsH3$sps)
plotdat.kfuncsH3 <-do.call('rbind', lapply(resultsH3$sps, function(x){
  ki <- Kest(x, r=seq(from=0,to=6, by=.1), correction='iso')
  L <- K2L(ki$iso, ki$r)
  dat <- data.frame(distance=ki$r, K=ki$iso, L=L)
  return(dat)
}))



pl1<-pl1+
  geom_smooth(data=plotdat.kfuncsH3, aes(x=distance, y=L),colour="black", linetype=4)
pl1

listH4
KfuncsH4 <- lapply(listH4, Kest, r=xrange, correction=NULL, ratio=FALSE)
## put in matrix - distances are columns, reps are rows
KfuncsMatH4 <- do.call('rbind', lapply(KfuncsH4, function(x) x$trans))

#remove NA's
KfuncsMatH4
KfuncsMatH4<-KfuncsMatH4[-c(1,2,5:7,11,12,14:16,20:25, 27, 28),]
KfuncsMatH4
KfuncsMatH4<-KfuncsMatH4[,-c(29:61)]
KfuncsMatH4
KfuncsMatH4<-KfuncsMatH4[-c(5,9),]

## calculate weights
covwtsH4 <-  do.call('rbind', lapply(listH4, kfunc.weights.calc, type='nx',
                                     r=xrange, corr=NULL))

##remove same parts as above
covwtsH4 <- apply(covwtsH4, 2, function(x) {
  x <- x/mean(x)
  return(x)})


covwtsH4<-covwtsH4[-c(1,2,5:7,11,12,14:16,20:25, 27, 28),]
covwtsH4<-covwtsH4[,-c(29:61)]
covwtsH4<-covwtsH4[-c(5,9),]


interceptdatH4 <- data.frame(int=rep(1, nrow(KfuncsMatH4)))

## fit the models
lmmodsH4 <- try(lapply(1:ncol(KfuncsMatH4), function(j)
  kfunclm(k=KfuncsMatH4[,j], weights=covwtsH4[,j],
          dat=interceptdatH4, form='1')), silent=TRUE)
## do the bootstrapping for confidence intervals
lc <- matrix(1, ncol=1)
bsciH4 <- lm.t.boot(lmmodsH4, lincomb=lc, nsim=99, alpha=0.05,
                    simple.method=FALSE)

## pull out the parts needed for plotting and add 0s for distances
## where all reps were 0
lmKpredH4 <-bsciH4$lmKpred
lower.CIH4<-bsciH4$lower
upper.CIH4<-bsciH4$upper
## put 0s back in manually where there were only 0s in K -
## models often fail and give us the same result

resultsH4<-list(lmK=lmmodsH4, lmKpred=lmKpredH4,lower=lower.CIH4, upper=upper.CIH4, sps=listH4)

plotdat.kfuncsH4<-do.call('rbind',resultsH4$sps)
plotdat.kfuncsH4 <-do.call('rbind', lapply(resultsH4$sps, function(x){
  ki <- Kest(x, r=seq(from=0,to=6, by=.1), correction='iso')
  L <- K2L(ki$iso, ki$r)
  dat <- data.frame(distance=ki$r, K=ki$iso, L=L)
  return(dat)
}))

summary(lm(plotdat.kfuncsH1$L~plotdat.kfuncsH1$distance))
summary(lm(plotdat.kfuncsH2$L~plotdat.kfuncsH2$distance))
summary(lm(plotdat.kfuncsH3$L~plotdat.kfuncsH3$distance))
summary(lm(plotdat.kfuncsH4$L~plotdat.kfuncsH4$distance))


library("lme4")


write.csv(plotdat.kfuncsH1, file="H1.csv")
write.csv(plotdat.kfuncsH2, file="H2.csv")
write.csv(plotdat.kfuncsH3, file="H3.csv")
write.csv(plotdat.kfuncsH4, file="H4.csv")

pl1<-pl1+
  geom_smooth(data=plotdat.kfuncsH4, aes(x=distance, y=L),colour="black", linetype=5)+
  coord_cartesian(ylim = c(-13, 13))
  



png(file="Figure2-CommBySize-11.8.15.png",width=1200,height=800, res=300, bg="transparent")
pl1
dev.off()




