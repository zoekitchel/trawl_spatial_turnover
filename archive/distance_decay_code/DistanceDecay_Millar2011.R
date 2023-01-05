#Jackknife function
#Inefficient, because it recomputes distances each time (but is still quick)
Jknife=function(distances,similarities,measure="binary",dist="dist") {
  dist=match.fun(dist)
  if(dim(distances)[1]!=nrow(as.matrix(1-beta.pair.abund$beta.bray.bal)))
    stop("\nDifferent number of sites detected in input data")
  nsites=dim(distances)[1]
  Results=matrix(NA,nrow=nsites,ncol=2)
  for(i in 1:nsites) {
    d=as.vector(as.dist(distances.dt[-i]))
    s=1-as.vector(1-dist(counts[-i,],method=measure))
    #log-binomial fit
    g=glm(s~d,family=binomial(link=log))
    Results[i,]=coef(g) }
  #Get parameters of interest
  Results=cbind(Results[,1],-Results[,2],
                exp(Results[,1]),-log(2)/Results[,2])
  colnames(Results)=c("a","beta","s0","halfd")
  Summary=rbind(apply(Results,2,mean),
                sqrt(diag(var(Results)*(nsites-1)*(1-1/nsites))))
  rownames(Summary)=c("Jnife mean","Jknife s.e.")
  return(list(Fits=Results,Summary=Summary))
}
    
#Bootstrap function
#Inefficient, because it recomputes distances each time 
# (and hence not so quick for a very large choice of the number of simulations)
#Use the argument zeros=F to remove like pairs. The bootstrap variance is then
#weighted 
Bstrap=function(gradient,counts,nboots=1000,measure="binary",dist="dist",
                like.pairs=T) {
  dist=match.fun(dist)
  if(length(gradient)!=nrow(as.matrix(counts)))
    stop("\nDifferent number of sites detected in input data")
  nsites=length(gradient)
  Results=matrix(NA,nrow=nboots,ncol=2)
  nPairs=rep(NA,nboots)
  for(i in 1:nboots) {
    indices=sample(1:nsites,replace=T)
    d=as.vector(dist(gradient[indices],method="euclidean"))
    s=as.vector(1-dist(counts[indices,],method=measure))
    df=data.frame(d,s)
    n=nrow(df)
    #remove like site pairs (zero distance) if requested
    if(!like.pairs) df=subset(df,subset=c(d>0))
    n.used=nrow(df)
    if(n.used<2) stop("Less than 2 site-pairs. Can not fit glm")
    #log-binomial fit
    g=glm(s~d,family=binomial(link=log),data=df)
    Results[i,]=coef(g); nPairs[i]=n.used }
  #Get parameters of interest
  Results=cbind(Results[,1],-Results[,2],
              exp(Results[,1]),-log(2)/Results[,2])
  colnames(Results)=c("a","beta","s0","halfd")
  if(!like.pairs) {
   cat("\n",100*(n-mean(nPairs))/n,"% removed due to like pairs")
   cat("\n Bstrap s.e.s have been corrected for the smaller no. of pairs\n\n") }
  Variances=apply(Results,2,wt.var,w=nPairs/n)*mean(nPairs)/n
  Summary=rbind(apply(Results,2,weighted.mean,w=nPairs/n),sqrt(Variances))
  rownames(Summary)=c("Bstrap mean","Bstrap s.e.")
  return(list(Fits=Results,Summary=Summary,
               CtrlList=list(nboots=nboots,like.pairs=like.pairs)))
}

#Test with latitude 33 data
fish_33.df=read.csv("Lat33.csv",header=T)
counts=fish_33.df[,-1] #First column is depth
#counts = dist matrix of species abundances
gradient=fish_33.df$Depth
#gradient= vector of depths

###ACTUAL FIT###
d=as.vector(dist(gradient,method="euclidean"))
#i already have this vector of distances in "distances"
s=as.vector(1-dist(counts,method="binary"))
#I already have similarities in bray_balanced.dt$similarities
glm(s~d,family=binomial(link=log))

###JACKKNIFE###
J=Jknife(gradient,counts) #Takes about 4s
J$Summary
#                      a         beta         s0     halfd
#Jnife mean  -0.80747900 0.0028483647 0.44599706 243.36712
#Jknife s.e.  0.06946021 0.0002014233 0.03107768  17.01097

library(corpcor) #To get weighted variance functions wt.var
###BOOTSTRAP###
B1=Bstrap(gradient,counts,nboots=100) #Does about 1000 per min.
B1$Summary
#From a 10000 bootstrap run
#                      a         beta         s0     halfd
#Bstrap mean -0.69693365 0.0031785115 0.49907891 219.27556
#Bstrap s.e.  0.06241046 0.0002384570 0.03106484  16.13263

###BOOTSTRAP with like site pairs removed###
B2=Bstrap(gradient,counts,nboots=100,like.pairs=F) #Does about 1000 per min.
B2$Summary
#From a 10000 bootstrap run
#                     a         beta         s0     halfd
#Bstrap mean -0.8045717 0.0028810080 0.44838149 241.86885
#Bstrap s.e.  0.0698081 0.0002091287 0.03105285  17.41755
#Note the BIG difference in s0 and halfd from the choice of using
#like pairs in the bootstrap or not.
#With like pairs removed, bootstrap results are very similar to jackknife.


