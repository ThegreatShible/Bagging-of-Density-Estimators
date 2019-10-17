library(parallel)
library(ks)

source("functions2.r")

simulaciones=function(n=100,M=10,B=150){
  # Kernel smoothing
  require(ks)
  # Implementation of distributions for nonparametric density estimation
  require(benchden)
  AA=matrix(0,nrow=6,ncol=7)
  colnames(AA)=c("H","FP","Kde", "BagHist", "BagHistFP", "Bagkde", "Rash")
  rownames(AA)=c("normal","chi2","mezcla1","mezcla2","bart","triangular")
  for(i in 1:M){
      
    A=matrix(NA,nrow=6,ncol=7)
    j=1
    for(nummodel in c(1,3,5,7,8,11)){
      
      sample=gendata(nummodel,n)
      bopt=bropt(sample$train)$opt
      zz=hist(sample$train,breaks=mybreaks(sample$train,nbr=bopt),plot=F)
      
      #bhist=BagHistfp(xx=sample$train,grille=sample$test,nbr = bopt, B)
      bhist=BagHistfp(xx=sample$train,grille=sample$test, B)
      modelrash=rash(sample$train,grille=sample$test,nbr = bopt, B)
      #modelavshift(sample$train,sample$test,nbr=bopt,M=B)
      modelbagkde <- Bagkde(xx = sample$train, grille = sample$test, B)
      
      A[j,]=c(error(sample$dobs,predict.hist(zz,sample$test))[1],
               error(sample$dobs,approxfun(x=zz$mids,y=zz$density)(sample$test))[1],    
               error(sample$dobs,onekdeucv(sample$train,sample$test))[1],
               error(sample$dobs,bhist$bh)[1],
               error(sample$dobs,bhist$bhfp)[1], 
               #error(sample$dobs,modelavshift)[1]
               error(sample$dobs,modelbagkde)[1],
               error(sample$dobs,modelrash)[1] )
      
      #MISE[j,]=(sample$dobs-bhist$bhfp)^2 +
      j=j+1
    }
    AA=AA+A
  }
  AA=AA/M
}

cls <- makeCluster(detectCores() - 1)
vars2export <- c("BagHistfp", "Bagkde", "bropt", "broptfp", "dtriangle", "error", "gendata",         
                "mel", "melange", "mybreaks", "onekdeucv",
                "predict.hist", "predict.hist.x",   "rash", 
                "riskfp", "riskhist", "rtriangle", "simulaciones")   
clusterExport(cls, vars2export)

res <- simulaciones(n = 20, M = 10, B = 20)
res
#res <- parLapplyLB(cls, rev(c(20, 50, 100, 200, 500, 1000, 2000)),
#              function(n) simulaciones(n = n, M = 100, B = 200))
#)

stopCluster(cls)

save(res, file = "resultados_mise-par2.Rdata")
