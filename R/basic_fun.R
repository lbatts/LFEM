##-----------------------------
## Function for EM algo with VB using TMB
## Author: LB 
## Date: 07/12/2017
## This version allows multiple surveys from diffrent times of year
## Note L and l inputs should be approximate for the survey which comes first alphabetically
## 
##-----------------------------
   basic.LFEM<-function(year0,no.years, age1, L, l, k.reparam, sigma.start ,No.comp, Lengths,niter,SD.type,rel.tolerance,dllroot,sub.Obs.lim)    {
     
     ### year 0 should be vector of length(no.surveys)
     ### age zero should be vector of all surveys age zeros, this idicates where to calculate mus from
     
     #####set up Mu's. If no variation between years than only need one set one Mus
     no.surveys<-length(no.years)  
     
     EN <- NA  #set up object for entropy 
     
     if(missing(SD.type)){SD.type <- 4}
     if(SD.type==4 && length(sigma.start)>1){
       stop("MESSAGE:  One sigma.start value needed for constant SD")
     } else if (SD.type==3 && length(sigma.start)==1){
       stop("MESSAGE: A vector of two sigma.start values needed for this SD type")
     }
     
     no.lambda.param <- sum(no.years*(No.comp-1))
     
     
     
     Lengths_matrix<-data.matrix(Lengths)####need this matrix for TMB model...converts survey characters to numbers in alphbetical order     
     Lengths_matrix[,2]<-Lengths_matrix[,2] - (min(year0))    ###sets up "years" column for use within TMB model..i.e. first year == 0 
     Lengths_matrix[,1]<-Lengths_matrix[,1] - 1               ####sets up "surveys" coumn for use within TMB model...i.e. first survey==0
     Lengths_matrix<-Lengths_matrix[order(Lengths_matrix[,1],Lengths_matrix[,2]),] #####order numeric matrix by first column (Survey) the secound column (Year). This means individual obs are the same order as Lengths (which follows)
     Lengths_matrix[is.na(Lengths_matrix)]<-0   ##This just means that if there is only one survey or one year there isn't NAs in matrix
     Lengths<-dlply(Lengths, .(Survey,Year))   ####set up data for function   ####NOTE...surveys set up in alphabetical order
     samplesize <- sum(Lengths_matrix[,4])
     no.years.covered<-  length(levels(factor(Lengths_matrix[,2]))) 
     
     
     #############################################Mus initialised     
     
     L.em<-L 
     l.em<-l    
     k.reparam.em<-k.reparam
     lambda.start <- rep(1 / No.comp, No.comp)
     #lambda.start<-as.brob(lambda.start)    ##no use as cannot convert into array
     lambda.array<-array(lambda.start,dim=c(length(levels(factor(Lengths_matrix[,2]))),No.comp,no.surveys))
     
     lambda.array.plus <- lambda.array
     sigma.em<-sigma.start
     
     
     Components<-(1:No.comp)        ##sets up component vector for TMB model
     surveyyear<-Lengths_matrix[,1:2]      ##sets up matrix for extracting parameter values specific for survey and year in TMB model
     
     #####load the specific TMB file and set up sigma.em
     if(SD.type==3){
       
       
       linf.em<-(L.em-(l.em*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))
       K.em <- -log(k.reparam.em)  
       tzero.em<-  age1[1]-((1/log(k.reparam.em))*log((L.em-l.em)/(L.em-(l.em*k.reparam.em^((No.comp-1)))))) 
       
       mu.em.array<-array(NA,dim=c(1,No.comp,no.surveys))
       
       
       
       mu.em.array[1,1,]<- linf.em*(1-exp(-K.em*(age1-tzero.em)))
       mu.em.array[1,No.comp,]<- linf.em*(1-exp(-K.em*((age1+(No.comp-1))-tzero.em)))
       
       for(i in 2:(No.comp)){
         mu.em.array[,i,] <- mu.em.array[,i-1,] + ((mu.em.array[1,No.comp,] - mu.em.array[1,1,])*(((k.reparam^(i-2))-(k.reparam^(i-1)))/(1-(k.reparam^((No.comp-1))))))
       }      
       
       ##schnute fournier SD formula
       dyn.load(dynlib(paste0(dllroot,"linearSD")))
       s.em<-sigma.em[1]
       S.em<-sigma.em[2]
       sd.em<-NA
       
       sd.em[1]<-s.em
       for(i in 2:No.comp){
         sd.em[i] <- sd.em[i-1] + ((S.em - s.em)*(((k.reparam^(i-2))-(k.reparam^(i-1)))/(1-(k.reparam^(No.comp-1)))))
       }      
       
     }else if(SD.type==4){
       dyn.load(dynlib(paste0(dllroot,"constantSD")))
       
       
       linf.em<-(L.em-(l.em*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))
       K.em <- -log(k.reparam.em)  
       tzero.em<-  age1[1]-((1/log(k.reparam.em))*log((L.em-l.em)/(L.em-(l.em*k.reparam.em^((No.comp-1)))))) 
       
       mu.em.array<-array(NA,dim=c(1,No.comp,no.surveys))
       
       
       
       mu.em.array[1,1,]<- linf.em*(1-exp(-K.em*(age1-tzero.em)))
       mu.em.array[1,No.comp,]<- linf.em*(1-exp(-K.em*((age1+(No.comp-1))-tzero.em)))
       
       for(i in 2:(No.comp)){
         mu.em.array[,i,] <- mu.em.array[,i-1,] + ((mu.em.array[1,No.comp,] - mu.em.array[1,1,])*(((k.reparam^(i-2))-(k.reparam^(i-1)))/(1-(k.reparam^((No.comp-1))))))
       }      
       
       sd.em<-rep(sigma.em,No.comp)
     }
     
     
     obs.llike <- rep(NA, niter)
     K.obs <- rep(NA,niter)
     Linf.obs <- rep(NA,niter)
     k.reparam.obs <- rep(NA,niter)
     L.obs <- rep(NA,niter)
     
     worker.dens<-matrix(NA,ncol=No.comp,nrow=dim(Lengths_matrix)[1])
     tau.mat<-matrix(NA,ncol=No.comp,nrow=dim(Lengths_matrix)[1])
     
     for(k in 1:niter){ 
       
       print(k)
       
       ##--------
       ## E-STEP
       ##--------
       
       
       
       #####################
       for(i in 1:dim(Lengths_matrix)[1]){
         
         worker.dens[i,] <- lambda.array[((Lengths_matrix[i,2])+1),,((Lengths_matrix[i,1])+1)] * dnorm(Lengths_matrix[i,3],mean=mu.em.array[1,,((Lengths_matrix[i,1])+1)],sd=sd.em)
         
       }
       
       
       obs.llike[k] <- sum(Lengths_matrix[,4]*log(rowSums(worker.dens)))
       K.obs[k] <- K.em
       Linf.obs[k] <- linf.em
       k.reparam.obs[k] <- k.reparam.em
       L.obs[k] <- L.em
       
       
       
       tau.mat <- worker.dens[,]/rowSums(worker.dens[,])
       
       
       
       for (i in 1:no.surveys){
         for (j in 1:no.years.covered){
           
           lambda.array.plus[j,,i] <- colSums(Lengths_matrix[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,4]*tau.mat[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,]) / (sum(Lengths_matrix[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,4]*tau.mat[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,])) 
           
         }
       } 
       
       
       EN<- -sum(Lengths_matrix[,4]*(tau.mat*log(ifelse(tau.mat==0,1.1e-323,tau.mat))))
       
       
       if(k > 2 ){
         if(abs(obs.llike[k] - obs.llike[k-1]) <  abs(obs.llike[k-1] * rel.tolerance)){
           
           
           if(SD.type==4) {
             
             dyn.load(dynlib(paste0(dllroot,"constantSD_OBSLL")))
             loglik_nonvar_OBSLL <- MakeADFun(
               data = list(Lengths=Lengths_matrix,surveyyear=surveyyear,lambda=lambda.array,Components=Components,agezero=age1), 
               parameters = list(log_l=log(l.em),log_L=log(L.em),logit_k_reparam=qlogis(k.reparam.em),log_sigma=log(sigma.em)),  
               DLL = "constantSD_OBSLL",silent=T,
               hessian = T)
             
             opt <- nlminb(start=loglik_nonvar_OBSLL$par,objective=loglik_nonvar_OBSLL$fn,gradient=loglik_nonvar_OBSLL$gr,silent=F)
             
             rep <- sdreport(loglik_nonvar_OBSLL)
             srep <- summary(rep)
             
             
           }else if(SD.type==3) {
             
             dyn.load(dynlib(paste0(dllroot,"linearSD_OBSLL")))
             loglik_var_OBSLL <- MakeADFun(
               data = list(Lengths=Lengths_matrix,surveyyear=surveyyear,lambda=lambda.array,Components=Components,agezero=age1), 
               parameters = list(log_l=log(l.em),log_L=log(L.em),logit_k_reparam=qlogis(k.reparam.em),log_s=log(s.em),log_S=log(S.em)),  
               DLL = "linearSD_OBSLL",silent=T,
               hessian = T)
             
             opt <- nlminb(start=loglik_var_OBSLL$par,objective=loglik_var_OBSLL$fn,gradient=loglik_var_OBSLL$gr,silent=F)
             
             rep <- sdreport(loglik_var_OBSLL)
             srep <- summary(rep)
             
             
             
           }
           
           
           obs.llike <- obs.llike[!is.na(obs.llike)]
           sub.Obs<-sum(Lengths_matrix[Lengths_matrix[,3]<=sub.Obs.lim,4]*log(rowSums(worker.dens[Lengths_matrix[,3]<=sub.Obs.lim,])))
           
           
           return(list(obs.llike=obs.llike,k.reparam.obs=k.reparam.obs,L.obs=L.obs,K.obs=K.obs,Linf.obs=Linf.obs,Mu=mu.em.array,Sd=sd.em,Lambda=lambda.array,
                       k.reparam=k.reparam.em,K=K.em,Linf=linf.em,tzero=tzero.em,l=l.em,L=L.em,sigma=sigma.em, Lambda.params=no.lambda.param,sample.size=samplesize,Entropy=EN,age1=age1,srep=srep,sub.Obs=sub.Obs))
           #stop("converged")
         }
       }   
       
       
       
       
       lambda.array<-ifelse(lambda.array == 0,1.1e-323,lambda.array)   ###hack that gets around log(lambda) issues in TMB when Lambda value ==0          
       
       #
       #
       #
       
       
       ######################
       #####################
       ####Maximise the conditional expectation###########
       if(SD.type==3){
         
         loglik_var_schun <- MakeADFun(
           data = list(Lengths=Lengths_matrix,tau=tau.mat,surveyyear=surveyyear,lambda=lambda.array,Components=Components,agezero=age1), 
           parameters = list(log_l=log(l.em),log_L=log(L.em),logit_k_reparam=qlogis(k.reparam.em),log_s=log(s.em),log_S=log(S.em)),  
           DLL = "linearSD",silent=T)
         
         
         opt <- nlminb(start=loglik_var_schun$par,objective=loglik_var_schun$fn,gradient=loglik_var_schun$gr,silent=T)
         
         
         par1<-opt$par
         
         k.reparam.em <- plogis(par1[3])
         l.em <- exp(par1[1])
         L.em<-exp(par1[2])
         s.em<-exp(par1[4]) 
         S.em<-exp(par1[5])           
         
         
         linf.em<-(L.em-(l.em*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))
         K.em <- -log(k.reparam.em)  
         tzero.em<-  age1[1]-((1/log(k.reparam.em))*log((L.em-l.em)/(L.em-(l.em*k.reparam.em^((No.comp-1)))))) 
         
         
         mu.em.array[1,1,]<- linf.em*(1-exp(-K.em*(age1-tzero.em)))
         mu.em.array[1,No.comp,]<- linf.em*(1-exp(-K.em*((age1+(No.comp-1))-tzero.em)))
         
         for(i in 2:(No.comp)){
           mu.em.array[,i,] <- mu.em.array[,i-1,] + ((mu.em.array[1,No.comp,] - mu.em.array[1,1,])*(((k.reparam.em^(i-2))-(k.reparam.em^(i-1)))/(1-(k.reparam.em^((No.comp-1))))))
         }      
         
         
         
         sd.em<-NA
         sd.em[1]<-s.em
         for(i in 2:No.comp){
           sd.em[i] <- sd.em[i-1] + ((S.em - s.em)*(((k.reparam.em^(i-2))-(k.reparam.em^(i-1)))/(1-(k.reparam.em^(No.comp-1)))))
         }    
         
         sigma.em[1] <- s.em ###this just means that s.em and S.em will be in output
         sigma.em[2] <- S.em ###
         
       }else if(SD.type==4) {
         
         loglik_nonvar <- MakeADFun(
           data = list(Lengths=Lengths_matrix,tau=tau.mat,surveyyear=surveyyear,lambda=lambda.array,Components=Components,agezero=age1), 
           parameters = list(log_l=log(l.em),log_L=log(L.em),logit_k_reparam=qlogis(k.reparam.em),log_sigma=log(sigma.em)),  
           DLL = "constantSD",silent=T)
         
         opt <- nlminb(start=loglik_nonvar$par,objective=loglik_nonvar$fn,gradient=loglik_nonvar$gr,silent=T)
         
         par1<-opt$par
         
         k.reparam.em <- plogis(par1[3])
         l.em <- exp(par1[1])
         L.em<-exp(par1[2])
         sigma.em<-exp(par1[4])          
         linf.em<-(L.em-(l.em*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))
         K.em <- -log(k.reparam.em)  
         tzero.em<-  age1[1]-((1/log(k.reparam.em))*log((L.em-l.em)/(L.em-(l.em*k.reparam.em^((No.comp-1)))))) 
         
         
         mu.em.array[1,1,]<- linf.em*(1-exp(-K.em*(age1-tzero.em)))
         mu.em.array[1,No.comp,]<- linf.em*(1-exp(-K.em*((age1+(No.comp-1))-tzero.em)))
         
         for(i in 2:(No.comp)){
           mu.em.array[,i,] <- mu.em.array[,i-1,] + ((mu.em.array[1,No.comp,] - mu.em.array[1,1,])*(((k.reparam.em^(i-2))-(k.reparam.em^(i-1)))/(1-(k.reparam.em^((No.comp-1))))))
         }      
         
         
         
         sd.em <-rep(sigma.em,No.comp)
       }
       
       ##lambda.array.plus calculated earlier for ease but not used in maximisation step now becomes lambda.array
       lambda.array <- lambda.array.plus
       
       
       
       if(k==niter){
         print("Reached max iterations and not converged")
         print(paste("Result is after ",niter," M steps..."))        
         return(list(obs.llike=obs.llike,Mu=mu.em.array,Sd=sd.em,Lambda=ifelse(lambda.array==(1 / No.comp),0,lambda.array),
                     k.reparam=k.reparam.em,K=K.em,Linf=linf.em,tzero=tzero.em,l=l.em,L=L.em,sigma=sigma.em))
       }
       
     }
   }