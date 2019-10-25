library(parallel)
library(ks)
library(devtools)
#install("../CppFunctions")
library(CppFunctions)

source("functions2.r")

simulaciones=function(n=100,M=10,B=150){
  # Kernel smoothing
  require(ks)
  # Implementation of distributions for nonparametric density estimation
  require(benchden)
  res=matrix(0,nrow=6,ncol=7)
  colnames(res)=c("H","FP","Kde", "BagHist", "BagHistFP", "Bagkde", "Rash")
  rownames(res)=c("normal","chi2","mezcla1","mezcla2","bart","triangular")
  for(i in 1:M){
    err=matrix(NA,nrow=6,ncol=7)
    j=1
    for(nummodel in c(1,3,5,7,8,11)){
      samples=gendata(nummodel,n)
      bopt=broptRcpp(samples$train)$opt
      his=hist_rcpp(samples$train,breaks=mybreaks(samples$train,nbr=bopt))
      
      bhist=BagHistfp(xx=samples$train,grille=samples$test, B)
      modelrash=rash(samples$train,grille=samples$test,nbr = bopt, B)
      #modelavshift(samples$train,samples$test,nbr=bopt,M=B)
      modelbagkde <- Bagkde(xx = samples$train, grille = samples$test, B)
      
      err[j,]=c(error(samples$dobs,predict.hist(his,samples$test))[1],
              error(samples$dobs,approxfun(x=his$mids,y=his$density)(samples$test))[1],    
              error(samples$dobs,onekdeucv(samples$train,samples$test))[1],
              error(samples$dobs,bhist$bh)[1],
              error(samples$dobs,bhist$bhfp)[1], 
              #error(samples$dobs,modelavshift)[1]
              error(samples$dobs,modelbagkde)[1],
              error(samples$dobs,modelrash)[1] )
      
      #MISE[j,]=(samples$dobs-bhist$bhfp)^2 +
      j=j+1
    }
    res=res+err
  }
  res=res/M
}

cls <- makeCluster(detectCores() - 1)
vars2export <- c("BagHistfp", "Bagkde", "bropt", "broptfp", "dtriangle", "error", "gendata",         
                 "mel", "melange", "mybreaks", "onekdeucv",
                 "predict.hist", "predict.hist.x",   "rash", 
                 "riskfp", "riskhist", "rtriangle", "simulaciones")   
clusterExport(cls, vars2export)
clusterEvalQ(cls, library(CppFunctions))
#Rprof()
res <- simulaciones(n = 50, M = 100, B = 20)
#Rprof(NULL)
#summaryRprof()
res
#res <- parLapplyLB(cls, rev(c(20, 50, 100, 200, 500, 1000, 2000)),
#              function(n) simulaciones(n = n, M = 100, B = 200))
#)

stopCluster(cls)

save(res, file = "resultados_mise-par2.Rdata")
