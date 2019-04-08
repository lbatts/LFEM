 ##-----------------------------
## code for AIC and BIC on basic fun
## Author: LB 
## Date: 23/10/2017
## 
##-----------------------------
 
 
 
 
 
 
 ####set up BIC fucntion
bic<-function(model){
k= 3 + length(unique(model$sigma))+(model$Lambda.params)
n=model$sample.size

k*log(n) - 2*(model$obs.llike[length(model$obs.llike)])

}

     ####set up AIC fucntion
     aic<-function(model){
     k=3 + length(unique(model$sigma))+(model$Lambda.params)
     
     -2*(model$obs.llike[length(model$obs.llike)]) +2*k
     } 
     
     
     #sub AIC
     
     sub_aic<-function(model){
       k=3 + length(unique(model$sigma))+(model$Lambda.params)
       
       -2*(model$sub.Obs) +2*k
     } 
     
     #icl bic
     
     icl_bic <- function(model){
       k= 3 + length(unique(model$sigma))+(model$Lambda.params)
       n=model$sample.size
       EN = model$Entropy
       
     (2*EN) + (k*log(n)) - 2*(model$obs.llike[length(model$obs.llike)])
     
     
       
     }