
rm(list=ls())
source("functions2.r")
library(devtools)
#install("../CppFunctions")
library(CppFunctions)
library(parallel)
#####################
#Evol. error

evol.error <- function(modele=1,n=100,B=200,M=5){
    
    A=matrix(NA,nrow=M,ncol=B)
    C=matrix(NA,nrow=M,ncol=B)
    D=matrix(NA,nrow=M,ncol=B)
    E=matrix(NA,nrow=M,ncol=B)
    for(kk in 1:M){
      print(kk)
      dd=gendata(modele,n)
      bopt=broptRcpp(dd$train)$opt
      
      bhfp=BagHistfp.err(xx=dd$train,grille=dd$test,B=B,dobs=dd$dobs)
      A[kk,]=bhfp$erreurbh[,1]
      C[kk,]=bhfp$erreurbhfp[,1]
      D[kk,]=BagKDE.err(dd$train, dd$test,B,dd$dobs)$erreur[,1]
      E[kk,]=rash.err(dd$train,dd$test,nbr=bopt, B,dobs=dd$dobs)$erreur[,1]
    }
list(A=A,C=C,D=D,E=E)
    }
vars2export <- c("BagHistfp.err","BagKDE.err","BagHistfp", "kde", "Bagkde","evol.error",
                                 "rash.err", "gendata", "mybreaks", "predict.hist", "predict.hist.x",
                                 "error","melange","mel", "rberdev", "dberdev", "rtriangle", "dtriangle", "ind")   
cls <- makeCluster(detectCores() - 1)
clusterExport(cls, vars2export)
clusterEvalQ(cls,library(CppFunctions))
system.file(
  res <- parLapply(cls, c(1, 3, 5, 8, 11,13,20,21),
                   function(modele) evol.error(modele = modele, n=500, M = 100, B = 200))
)
stopCluster(cls)




save(res, file = "res_evolerror_new.Rdata")

