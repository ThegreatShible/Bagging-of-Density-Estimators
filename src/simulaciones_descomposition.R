
library("foreach")
library("iterators")
library("parallel")
rm(list = ls())

#B nombre d'echantillonage pour le bagging
# M nombre de d?composition du dataset avec chaque partie contenant n elements

innerFunc <- function(n,M,K,B,numModel) {
  res = matrix(0, nrow = 1, ncol = 21)
  for (m in 1:M) {
    distributionData =gendata(numModel,n)
    #optimal number of chunks in the histogram
    bopt=bropt(distributionData$train)$opt
    histogram=hist(distributionData$train,breaks=mybreaks(distributionData$train,nbr=bopt),plot=F)
    #O(distributionData.length * histogram.breaks)
    h=predict.hist(histogram,distributionData$test)
    fp=approxfun(x=histogram$mids,y=histogram$density)(distributionData$test)
    kde=onekdeucv(distributionData$train,distributionData$test)
    bhist = BagHistfp(xx=distributionData$train,grille=distributionData$test, B)
    estim1=bhist$bhfp
    estim2=bhist$bh
    estim3=Bagkde(xx=distributionData$train,grille=distributionData$test,B)
    estim4=rash(xx=distributionData$train,grille=distributionData$test, nbr=bopt,B)
    
    finalh=0
    finalh_b=0
    
    finalfp=0
    finalfp_b=0
    
    finalbh=0
    finalbh_b=0
    
    finalbhfp=0
    finalbhfp_b=0
    
    finalkde=0
    finalkde_b=0
    
    finalbagkde=0
    finalbagkde_b=0
    
    finalrash=0
    finalrash_b=0
    
    for (k in 1:K){
      #if(k%%20 == 0) cat(k,">>")
      sample=onlyTrainGendata(numModel,n)
      bopt=bropt(sample$train)$opt
      # -- Histogram
      histogram=hist(sample$train,breaks=mybreaks(sample$train,nbr=bopt),plot=F)
      h=predict.hist(histogram,distributionData$test)
      # -- FP
      fp=approxfun(x=histogram$mids,y=histogram$density)(distributionData$test)
      # KDE
      kde=onekdeucv(sample$train,distributionData$test)
      
      bagHist <- BagHistfp(xx=sample$train,grille=distributionData$test,B)
      # BagHist
      estim2l=bagHist$bh
      
      # BagFP
      estim1l=bagHist$bhfp
      
      #BagKde
      estim3l=Bagkde(xx=sample$train,grille=distributionData$test,B)
      
      #Rash
      estim4l=rash(xx=sample$train,grille=distributionData$test, nbr=bopt,B)
      
      finalh      = finalh      + (h - distributionData$dobs)^2
      finalfp     = finalfp     + (fp - distributionData$dobs)^2 
      finalkde    = finalkde    + (kde - distributionData$dobs)^2
      
      finalbh   = finalbh   + (estim2l - distributionData$dobs)^2 
      finalbhfp   = finalbhfp   + (estim1l - distributionData$dobs)^2 
      finalbagkde = finalbagkde + (estim3l-distributionData$dobs)^2
      finalrash = finalrash + (estim4l-distributionData$dobs)^2
      
      finalbh_b = finalbh_b + estim2l
      finalbagkde_b=finalbagkde_b + estim3l
      finalbhfp_b = finalbhfp_b + estim1l
      finalrash_b = finalrash_b + estim4l
      finalkde_b  = finalkde_b  + kde
      finalh_b    = finalh_b    + h
      finalfp_b   = finalfp_b   + fp 
      
      #print(finalbhfp)
    }
    
    finalh      = finalh/K
    finalbh    = finalbh/K
    
    finalfp     = finalfp/K
    finalbhfp   = finalbhfp/K
    finalkde    = finalkde/K
    finalbagkde= finalbagkde/K
    finalrash = finalrash/K
    
    finalh_b    = finalh_b/K
    finalbh_b = finalbh_b/K
    finalfp_b   = finalfp_b/K
    finalbhfp_b = finalbhfp_b/K
    finalkde_b  = finalkde_b/K
    finalbagkde_b  = finalbagkde_b/K
    finalrash_b = finalrash_b/K
    
    res = res + c(mean(finalh                   , na.rm = T),
                  mean((finalh_b - distributionData$dobs)^2   , na.rm = T),
                  mean((h-finalh_b)^2           , na.rm = T),
                  
                  
                  mean(finalfp                  , na.rm = T),
                  mean((finalfp_b - distributionData$dobs)^2  , na.rm = T),
                  mean((fp - finalfp_b)^2       , na.rm = T),
                  
                  mean(finalkde                 , na.rm = T),
                  mean((finalkde_b - distributionData$dobs)^2 , na.rm = T),
                  mean((kde - finalkde_b)^2     , na.rm = T),
                  
                  mean(finalbh                   , na.rm = T),
                  mean((finalbh_b - distributionData$dobs)^2   , na.rm = T),
                  mean((estim2-finalbh_b)^2      , na.rm = T),
                  
                  mean(finalbhfp                , na.rm = T),
                  mean((finalbhfp_b - distributionData$dobs)^2, na.rm = T),
                  mean((estim1 - finalbhfp_b)^2  , na.rm = T),
                  
                  mean(finalbagkde                 , na.rm = T),
                  mean((finalbagkde_b - distributionData$dobs)^2 , na.rm = T),
                  mean((estim3 - finalbagkde_b)^2  , na.rm = T),
                  
                  mean(finalrash                 , na.rm = T),
                  mean((finalrash_b - distributionData$dobs)^2 , na.rm = T),
                  mean((estim4 - finalrash_b)^2  , na.rm = T) )
  }
  
  #res0[j,]=res/M
  # write.table(res0,"resbb_par.txt")
  #j=j+1
  res / M

  
}
funcArray <- c(1:7)
innerFunc2 <- function(n,M,K,B, numModel){
  
}

Kfunc <- function(n,K, testArray, numModel){
  #TODO implement funcArray, add rray package
  
  L <- length(testArray)
  samples <- matrix(onlyTrainGendata(numModel,n*K), nrow=K, ncol=n)
  predictions <- array(dim = c(7,K,L))
  tests1 <- array(testArray, dim=c(1,1,L))
  tests <-  rray_broadcast(tests1, c(7,K,L))
  #Todo check if for is the solution
  for (i in (1:length(funcArray))){
    f <- funcArray(i)
    funcRes <- f(samples)
    predictions[i,,] <- funcRes
  }
  
  
  
  
  
  
  
           
  
}
  
  
  
  
  
  
  
  
  
  
  


sesgo_par = function(n = 100, M = 5, K = 10, B = 10){

    res0 <- foreach::foreach(numModel = c(1, 3, 5, 8, 11, 13, 20, 21), .combine = rbind, 
                            .export = c("gendata", "bropt", "riskhist", "mybreaks", 
                                        "broptfp", "riskfp", "ind",
                                        "predict.hist", "predict.hist.x", "onekdeucv",
                                        "kde", "BagHistfp", "melange", "mel", "rberdev",
                                        "dberdev", "rtriangle", "dtriangle", "Bagkde", "rash"),
                           .packages = "ks") %dopar% {
                             innerFunc(n,M,K,B, numModel)
                             
                           }

  colnames(res0) = c("errorH",   "sesgoH",   "varH", 
                     "errorFP",  "sesgoFP",  "varFP",
                     "errorKde",  "sesgoKde",  "varKde",
                     "errorBH",   "sesgoBH",   "varBH",
                     "errorBFP", "sesgoBFP", "varBFP",
                     "errorBKde", "sesgoBKde", "varBKde",
                     "errorRash", "sesgoRash", "varRash")
  rownames(res0) = c("normal", "chi2", "mix1", "bart", "triangular", 
                     "rara1", "rara2", "rar3")
  
  return(res0)
}

library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
system.time(aa <- sesgo_par(n = 500, B = 200, K = 100, M = 100))
#system.time(aa <- sesgo_par(n = 500, B = 10, K = 5, M = 100))

save(aa, file = "resultados_descomposicion.Rdata")

