rm(list=ls())
library(spatstat)
data<-read.csv("PowdermillData-NoStems-1.20.15-NOHEIGHT.csv", header=TRUE)
str(data)
qlist=c(1:16)
data2<-lapply(qlist, function(x) data[data$PlotNum==x,])
data3<-lapply(data2, na.omit)


##Bagchi functions
K2L <- function(k, r, subtractr=TRUE){
  r[r==0] <- NA
  return( sqrt(k/pi)-r*subtractr)
}

Kpoisson <- function(r, L=FALSE, ...){
  k <- pi*r^2
  if(L)
    k <- K2L(k, r, ...)
  return(k)
}

Kthomas <- function(r, kappa, sigma, L=FALSE, ...){
  k <- pi*r^2 + (1/kappa)*(1-exp(-(r^2/(4*sigma^2))))
  if(L)
    k <- K2L(k, r, ...)
  return(k)
}

Kstrauss <- function(r, gamma, R, L=FALSE, ...){
  sapply(r, function(ri){
    if(ri <= R)
      K <- gamma*pi*ri^2
    else
      K <- pi*ri^2 - (1-gamma)*pi*R^2
    if(L)
      K <- K2L(K, ri, ...)
    return(K)
  })
}

resultCalc <- function(subplots, lincomb, nsim=99, alpha=0.05,
                       weights.type){
  require(abind)
  ##  cat('.') ## causes hangs - useful for testing
  ## calculate K funcsions
  Kfuncs <- lapply(subplots, Kest, r=0:25, correction='border', ratio=TRUE)
  ## put in matrix - distances are columns, reps are rows
  KfuncsMat <- do.call('rbind', lapply(Kfuncs, function(x) x$border))
  
  ## calculate weights
  covwts <-  do.call('rbind', lapply(subplots, kfunc.weights.calc, type=weights.type,
                                     r=0:25, corr='border'))
  
  covwts <- apply(covwts, 2, function(x) {
    x <- x/mean(x)
    return(x)})
  
  ## Some distances may have only 0s - no point modelling
  ## The mean k-function is 0!
  remove <- apply(KfuncsMat, 2, function(x) all(x==0))
  ## Remove these distances from Kfuncs mat and variance covariate
  KfuncsMat <- KfuncsMat[, !remove]; covwts <- covwts[, !remove]
  ## remove replicates where K function couldn't be calculated
  ## even at some distances (as this would cause bias)
  rmrows <- apply(KfuncsMat, 1, function(x) all(!is.na(x)))
  ## Remove the corresponding rows 
  KfuncsMat <- KfuncsMat[rmrows,]; covwts <- covwts[rmrows,]
  ##  Make data frame - in this case just an intercept
  interceptdat <- data.frame(int=rep(1, nrow(KfuncsMat)))
  
  ## fit the models
  lmmods <- try(lapply(1:ncol(KfuncsMat), function(j)
    kfunclm(k=KfuncsMat[,j], weights=covwts[,j],
            dat=interceptdat, form='1')), silent=TRUE)
  ## do the bootstrapping for confidence intervals
  bsci <- lm.t.boot(lmmods, lincomb=lincomb, nsim=nsim, alpha=alpha,
                    simple.method=FALSE)
  
  ## pull out the parts needed for plotting and add 0s for distances
  ## where all reps were 0
  lmKpred <- lower.CI <- upper.CI <- as.list(rep(0, length(remove)))
  lmKpred[!remove] <- bsci$lmKpred
  lower.CI[!remove] <-  bsci$lower; upper.CI[!remove] <- bsci$upper
  ## put 0s back in manually where there were only 0s in K -
  ## models often fail and give us the same result
  
  return(list(lmK=lmmods, lmKpred=lmKpred,
              lower=lower.CI, upper=upper.CI, sps=subplots))}




X0=vector(mode="list", length=16)
Y0=vector(mode="list", length=16)
Xrange=vector(mode="list", length=16)
Yrange=vector(mode="list", length=16)
extra0=vector(mode="list", length=16)
spatdata=vector(mode="list", length=16)
nomarks=vector(mode="list", length=16)
win=vector(mode="list", length=16)
real=vector(mode="list", length=16)
yin=vector(mode="list", length=16)

for (i in qlist){
  X0[[i]]<-data3[[i]]$Xcoord
  Y0[[i]]<-data3[[i]]$Ycoord
  Xrange[[i]]=as.vector(c(min(X0[[i]])-1, max(X0[[i]])+1))
  Yrange[[i]]=as.vector(c(min(Y0[[i]])-1, max(Y0[[i]])+1))
  extra0[[i]]<-data.frame(data3[[i]]$Species, data3[[i]]$SpeciesType,data3[[i]]$SpType2,data3[[i]]$AvgDamage, data3[[i]]$AvgDamageType)
  nomarks[[i]]<-ppp(X0[[i]],Y0[[i]],Xrange[[i]],Yrange[[i]], units="meters") 
  win[[i]]<-convexhull.xy(nomarks[[i]]$x,nomarks[[i]]$y)
  yin[[i]]<-as.polygonal(win[[i]])
  real[[i]]<-ppp(X0[[i]], Y0[[i]], window=yin[[i]])
  spatdata[[i]]<-ppp(X0[[i]],Y0[[i]],window=yin[[i]], marks=extra0[[i]])
}



results.KB.ratioweights <-  lapply(nomarks, ratio.weights.calc)
results.KB.abundanceweights<-lapply(nomarks, abundance.weights.calc)
results.KB.kfuncweights<-lapply(nomarks, kfunc.weights.calc, r=5, type='nx')
results.KB.kfunclm<-lapply(real,kfunclm, dat=real)
lc <- matrix(1, ncol=1)
results.calc.KB<-lapply(spatdata,kfunclm,subplot=spatdata,nsim=99, lincomb=lc, weights.type='nx')



library(GISTools)
library(foreign)
library(maptools)
shapedata<-readShapePoints("Points-KB.shp")
plot(shapedata)
boundary<-readShapePoly("PNRboundary.shp")
plot(boundary)
plot(shapedata, add=TRUE)
##writePolyShape(shapedata,"shapedata.shp")
sub<-subset(shapedata,Plot==5)
plot(boundary)
plot(sub)

pts=coordinates(sub)
pts=rbind(pts,pts[1,])
sp=SpatialPolygons(list(Polygons(list(polygon(pts)),1)))
subby<-summary(sub)

x<-convexhull.xy(nomarks[[1]]$x, nomarks[[1]]$y)
y<-as.polygonal(x)
subby<-ppp(X0[[1]],Y0[[1]],window=y)

x1<-centroid.owin(nomarks[[4]], as.ppp=TRUE)
y1<-disc(radius=40, centre=x1)
subby1<-ppp(X0[[4]],Y0[[4]],window=y1)
plot(subby1)

subby1$x

overstory<-subset(data, SpeciesType==1)
understory<-subset(data, SpeciesType==3)

ov.species.unique<-unique(overstory$Species)
ov.species.num<-length(ov.species.unique)

un.species.unique<-unique(understory$Species)
un.species.num<-length(un.species.unique)

aceru<-subset(overstory, Species=="aceru")
plotnum.unique<-unique(aceru$PlotNum)
plotnum<-length(plotnum.unique)
aceru.list<-lapply(plotnum.unique, function(x) aceru[aceru$PlotNum==x])

linbe<-subset(understory, Species=="linbe")
plotnum.unique.linbe<-unique(linbe$PlotNum)
plotnum.linbe<-length(plotnum.unique.linbe)
linbe.list<-lapply(plotnum.unique.linbe, function(x) linbe[linbe$PlotNum==x,])


understory.list<-lapply(un.species.unique, function(x) understory[understory$Species==x,])
num.und.species<-length(understory.list)
und.list.by.species<-vector(mode="list")
plotnum.by.species.unique<-unique(understory.list[[1]]$PlotNum)
plotnum.species<-length(plotnum.by.species.unique)
for (i in seq_along(num.und.species)){
  plotnum.by.species.unique<-unique(understory.list[[i]]$PlotNum)
  plotnum.species<-length(plotnum.by.species.unique)
  und.list.by.species[[i]]<-vector(mode="list", length=plotnum.by.species)
  und.list.by.species[[i]]<-lapply(plotnum.by.species.unique, function(x) understory.list[[i]][understory.list$PlotNum==x,]) 
}


x2.und<-vector(mode="list", length=un.species.num)
y2.linbe<-vector(mode="list", length=plotnum.linbe)
xrange2.linbe<-vector(mode="list", length=plotnum.linbe)
yrange2.linbe<-vector(mode="list", length=plotnum.linbe)
extra2.linbe<-vector(mode="list", length=plotnum.linbe)
nomarks.linbe<-vector(mode="list", length=plotnum.linbe)
x.linbe<-vector(mode="list", length=plotnum.linbe)
y.linbe<-vector(mode="list", length=plotnum.linbe)
real.linbe<-vector(mode="list", length=plotnum.linbe)
spatdata.linbe<-vector(mode="list", length=plotnum.linbe)
subby.linbe<-vector(mode="list", length=plotnum.linbe)
for (i in seq_along(plotnum.linbe)){
  x2.linbe[[i]]<-linbe.list[[i]]$Xcoord
  y2.linbe[[i]]<-linbe.list[[i]]$Ycoord
  xrange2.linbe[[i]]=c(min(x2.linbe[[i]])-1, max(x2.linbe[[i]])+1)
  yrange2.linbe[[i]]=c(min(y2.linbe[[i]])-1, max(y2.linbe[[i]])+1)
  extra2.linbe[[i]]<-data.frame(linbe.list[[i]]$Species, linbe.list[[i]]$SpeciesType, linbe.list[[i]]$SpType2,linbe.list[[i]]$AvgDamage, linbe.list[[i]]$AvgDamageType)
  nomarks.linbe[[i]]<-ppp(x2.linbe[[i]],y2.linbe[[i]],xrange2.linbe[[i]],yrange2.linbe[[i]])
  x[[i]]<-convexhull.xy(nomarks.linbe[[i]]$x, nomarks.linbe[[i]]$y)
  y[[i]]<-as.polygonal(x[[i]])
  subby[[i]]<-ppp(x2.linbe[[i]],y2.linbe[[i]],window=y[[i]])
}


x2.linbe<-vector(mode="list", length=plotnum.linbe)
y2.linbe<-vector(mode="list", length=plotnum.linbe)
xrange2.linbe<-vector(mode="list", length=plotnum.linbe)
yrange2.linbe<-vector(mode="list", length=plotnum.linbe)
extra2.linbe<-vector(mode="list", length=plotnum.linbe)
nomarks.linbe<-vector(mode="list", length=plotnum.linbe)
x.linbe<-vector(mode="list", length=plotnum.linbe)
y.linbe<-vector(mode="list", length=plotnum.linbe)
real.linbe<-vector(mode="list", length=plotnum.linbe)
spatdata.linbe<-vector(mode="list", length=plotnum.linbe)
subby.linbe<-vector(mode="list", length=plotnum.linbe)
for (i in seq_along(plotnum.linbe)){
  x2.linbe[[i]]<-linbe.list[[i]]$Xcoord
  y2.linbe[[i]]<-linbe.list[[i]]$Ycoord
  xrange2.linbe[[i]]=c(min(x2.linbe[[i]])-1, max(x2.linbe[[i]])+1)
  yrange2.linbe[[i]]=c(min(y2.linbe[[i]])-1, max(y2.linbe[[i]])+1)
  extra2.linbe[[i]]<-data.frame(linbe.list[[i]]$Species, linbe.list[[i]]$SpeciesType, linbe.list[[i]]$SpType2,linbe.list[[i]]$AvgDamage, linbe.list[[i]]$AvgDamageType)
  nomarks.linbe[[i]]<-ppp(x2.linbe[[i]],y2.linbe[[i]],xrange2.linbe[[i]],yrange2.linbe[[i]])
  x[[i]]<-convexhull.xy(nomarks.linbe[[i]]$x, nomarks.linbe[[i]]$y)
  y[[i]]<-as.polygonal(x[[i]])
  subby[[i]]<-ppp(x2.linbe[[i]],y2.linbe[[i]],window=y[[i]])
}

acepe<-subset(understory, Species=="acepe")
plotnum.unique.acepe<-unique(acepe$PlotNum)
acepe.list<-lapply(plotnum.unique.acepe, function(x) acepe[acepe$PlotNum==x])

##make counter #1
plist=c(1,3)

##separate into sublist by Species
data4<-lapply(plist, function(x) data[data$SpeciesType==x,])

##separate sublists into sublists by species
data5.0=vector(mode="list")
data5=vector(mode="list")
for (i in seq_along(plist)){
  num<-unique(data4[[i]]$Species)
  datacount<-length(num)
  data5.0[[i]]=vector(mode="list", length=datacount)
  data5[[i]]=vector(mode="list", length=datacount)
  data5.0[[i]]<-lapply(num, function(x) data4[[i]][data4[[i]]$Species==x,])
  data5[[i]]<-lapply(data5.0[[i]], na.omit)
  print(num)
}

data6.0=vector(mode="list")
data6=vector(mode="list")

for (i in seq_along(plist)){
  
  
  datacount<-length(data5[[i]])
  datacounter=c(1:datacount)
  print(datacount)
  data6.0[[i]]=vector(mode="list", length=datacount)
  data6[[i]]=vector(mode="list", length=datacount)
  
  for (j in seq_along(datacounter)){
    num2=unique(data5[[i]][[j]]$PlotNum)
    datacount2<-length(num2)
    print(datacount2)
    data6.0[[i]][[j]]=vector(mode="list", length=datacount2)
    data6[[i]][[j]]=vector(mode="list", length=datacount2)
    data6.0[[i]][[j]]<-lapply(num2, function(x) data5[[i]][[j]][data5[[i]][[j]]$PlotNum==x,])
    data6[[i]][[j]]<-lapply(data6.0[[i]][[j]], na.omit)
  }
}

X2=vector(mode="list")
Y2=vector(mode="list")
Xrange2=vector(mode="list")
Yrange2=vector(mode="list")
spatdata2=vector(mode="list")
extra2=vector(mode="list")
nomarks2=vector(mode="list")
modeled2=vector(mode="list")
figs2=vector(mode="list")
summ2=vector(mode="list")
win2=vector(mode="list")
yin2=vector(mode="list")
real2=vector(mode="list")

for (i in seq_along(plist)){
  
  
  datacount<-length(data6[[i]])
  datacounter=c(1:datacount)
  print(datacount)
  X2[[i]]=vector(mode="list", length=datacount)
  Y2[[i]]=vector(mode="list", length=datacount)
  Xrange2[[i]]=vector(mode="list", length=datacount)
  Yrange2[[i]]=vector(mode="list", length=datacount)
  spatdata2[[i]]=vector(mode="list", length=datacount)
  extra2[[i]]=vector(mode="list", length=datacount)
  nomarks2[[i]]=vector(mode="list", length=datacount)
  modeled2[[i]]=vector(mode="list", length=datacount)
  figs2[[i]]=vector(mode="list", length=datacount)
  summ2[[i]]=vector(mode="list", length=datacount)
  win2[[i]]=vector(mode="list", length=datacount)
  yin2[[i]]=vector(mode="list", length=datacount)
  real2[[i]]=vector(mode="list", length=datacount)
  
  for (j in seq_along(datacounter)){
    datacount2<-length(data6[[i]][[j]])
    print(datacount2)
    datacounter2<-c(1:datacount2)
    X2[[i]][[j]]=vector(mode="list", length=datacount2)
    Y2[[i]][[j]]=vector(mode="list", length=datacount2)
    Xrange2[[i]][[j]]=vector(mode="list", length=datacount2)
    Yrange2[[i]][[j]]=vector(mode="list", length=datacount2)
    spatdata2[[i]][[j]]=vector(mode="list", length=datacount2)
    extra2[[i]][[j]]=vector(mode="list", length=datacount2)
    nomarks2[[i]][[j]]=vector(mode="list", length=datacount2)
    modeled2[[i]][[j]]=vector(mode="list", length=datacount2)
    figs2[[i]][[j]]=vector(mode="list", length=datacount2)
    summ2[[i]][[j]]=vector(mode="list", length=datacount2)
    win2[[i]][[j]]=vector(mode="list", length=datacount2)
    yin2[[i]][[j]]=vector(mode="list", length=datacount2)
    real2[[i]][[j]]=vector(mode="list", length=datacount2)
    
    for (k in datacounter2){
      datacount3<-length(data6[[i]][[j]][[k]])
      
      if (datacount3>1){
        X2[[i]][[j]][[k]]<-data6[[i]][[j]][[k]]$Xcoord
        Y2[[i]][[j]][[k]]<-data6[[i]][[j]][[k]]$Ycoord
        Xrange2[[i]][[j]][[k]]=c(min(X2[[i]][[j]][[k]])-1, max(X2[[i]][[j]][[k]])+1)
        Yrange2[[i]][[j]][[k]]=c(min(Y2[[i]][[j]][[k]])-1, max(Y2[[i]][[j]][[k]])+1)
        extra2[[i]][[j]][[k]]<-data.frame(data6[[i]][[j]][[k]]$Species, data6[[i]][[j]][[k]]$SpeciesType,data6[[i]][[j]][[k]]$SpType2,data6[[i]][[j]][[k]]$AvgDamage, data6[[i]][[j]][[k]]$AvgDamageType)
        nomarks2[[i]][[j]][[k]]<-ppp(X2[[i]][[j]][[k]],Y2[[i]][[j]][[k]],Xrange2[[i]][[j]][[k]],Yrange2[[i]][[j]][[k]])
        win2[[i]][[j]][[k]]<-vector(mode="list",length=100)
        win2[[i]][[j]][[k]]<-convexhull.xy(nomarks2[[i]][[j]][[k]]$x, nomarks2[[i]][[j]][[k]]$y)
        print(win2[[i]][[j]][[k]])
        yin2[[i]][[j]][[k]]<-as.polygonal(win2[[i]][[j]][[k]])
        real2[[i]][[j]][[k]]<-ppp(X2[[i]][[j]][[k]],Y2[[i]][[j]][[k]],window=yin2[[i]][[j]][[k]])
        spatdata2[[i]][[j]][[k]]<-ppp(X2[[i]][[j]][[k]],Y2[[i]][[j]][[k]],window=yin2[[i]][[j]][[k]], marks=extra2[[i]][[j]][[k]])
        
      }
    } 
  }
}
}

## win2[[i]][[j]][[k]]<-centroid.owin(nomarks2[[i]][[j]][[k]], as.ppp=TRUE)
##yin2[[i]][[j]][[k]]<-disc(radius=100, centre=win2[[i]][[j]][[k]], unitname=c("metres"))

