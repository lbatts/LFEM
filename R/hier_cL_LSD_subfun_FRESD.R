hier_linearSD_L_FRESD<-function(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdL)
  
{
  
  #####set up Mu's. If no variation between years than only need one set one Mus
  no.surveys<-length(no.years)      
  if(length(sigma.start)<2){
    stop("Two sigma.start values needed")
  }
  
  EN <- NA
  Lengths_matrix<-data.matrix(Lengths)####need this matrix for TMB model...converts survey characters to numbers in alphbetical order     
  Lengths_matrix[,2]<-Lengths_matrix[,2] - (min(year0))    ###sets up "years" column for use within TMB model..i.e. first year == 0 
  Lengths_matrix[,1]<-Lengths_matrix[,1] - 1               ####sets up "surveys" coumn for use within TMB model...i.e. first survey==0
  Lengths_matrix<-Lengths_matrix[order(Lengths_matrix[,1],Lengths_matrix[,2]),] #####order numeric matrix by first column (Survey) the secound column (Year). This means individual obs are the same order as Lengths (which follows)
  Lengths_matrix[is.na(Lengths_matrix)]<-0   ##This just means that if there is only one survey or one year there isn't NAs in matrix
  Lengths<-dlply(Lengths, .(Survey,Year))   ####set up data for function   ####NOTE...surveys set up in alphabetical order
  
  no.lambda.param <- sum(no.years*(No.comp-1))
  samplesize <- sum(Lengths_matrix[,4])
  no.years.covered<-  length(levels(factor(Lengths_matrix[,2]))) 
  
  
  max.years <- (length(levels(factor(Lengths_matrix[,2]))))
  max.years.backforfill <- (length(levels(factor(Lengths_matrix[,2])))+(No.comp-1))
  max.years.comp <- (length(levels(factor(Lengths_matrix[,2])))+2*(No.comp-1))
  
  mu.arr.bff<-array(NA,dim=c(max.years.comp,No.comp,no.surveys))
  
  linf.em<-(L-(l*k.reparam^(No.comp-1)))/(1-(k.reparam^(No.comp-1)))## minus 2 as Final component is not restrivted to vb
  K.em <- -log(k.reparam)  
  tzero.em<-  age1[1]-((1/log(k.reparam))*log((L-l)/(L-(l*k.reparam^((No.comp-1)))))) 
  
  l.var<-linf.em*(1-exp(-K.em*(age1-tzero.em)))
  L.var<- linf.em*(1-exp(-K.em*((age1+(No.comp-1))-tzero.em)))
  
  l.par.vec.plus <- rep(l, max.years.backforfill)#temporary for ease of filling starting array, is cut down for algorithm
  
  L.par.vec.plus <- rep(L,max.years.backforfill)#temporary for ease of filling starting array, is cut down for algorithm
  
  
  #compile("hier_model_version_linearSD_adj.cpp")
  #compile("hier_model_version_linearSD_OBSLL_adj.cpp")
  dyn.load(dynlib(paste0(dllroot,"hier_cL_LSD")))
  dyn.load(dynlib(paste0(dllroot,"hier_cL_LSD_OBSLL")))
  
  for(j in 1:no.surveys){
    mu.arr.bff[1:max.years.backforfill,1,j]<- rep(l.var[j], max.years.backforfill)
    mu.arr.bff[No.comp:max.years.comp,No.comp,j]<- rep(L.var[j], max.years.backforfill)
    for(i in 2:max.years.comp){
      
      
      mu.arr.bff[i,2:No.comp,j]<-   sapply(2:No.comp,FUN=function(z){
        
        x <- i-z
        if(x>=0 && x<max.years.backforfill){
          l<-mu.arr.bff[x+1,1,j]
          L <-mu.arr.bff[x+No.comp,No.comp,j] 
        }else { 
          l<-0
          L <-0
        }
        mu.arr.bff[(i-1),(z-1),j] + ((L - l)*(((k.reparam^(z-2))-(k.reparam^(z-1)))/(1-(k.reparam^(No.comp-1)))))
      })
    }
  }
  #mu.arr.bff
  
  mu.arr <- array(data=mu.arr.bff[No.comp:max.years.backforfill,,], dim=c(max.years,No.comp,no.surveys))
  
  mu.em.array <- mu.arr
  
  l.em<-l
  L.em<-L
  
  L.par.vec.em<-L.par.vec.plus[No.comp:max.years.backforfill] 
  l.par.vec.em<-l.par.vec.plus[No.comp:max.years.backforfill]  
  l.mean.em<-l.em
  L.mean.em<-L.em
  
  k.reparam.em<-k.reparam
  lambda.start <- rep(1 / No.comp, No.comp)
  lambda.array<-array(lambda.start,dim=c(length(levels(factor(Lengths_matrix[,2]))),No.comp,no.surveys))
  lambda.array.plus <- lambda.array
  sigma.em<-sigma.start
  
  
  lL.mat.em<-matrix(NA,ncol=2,nrow=max.years)
  lL.mat.em[,1]<-0
  lL.mat.em[,2]<-0
  raw.sd.l.em<-sdl
  raw.sd.L.em<-sdL
  raw.rho.em<-0
  
  s.em<-sigma.em[1]
  S.em<-sigma.em[2]
  sd.em<-NA
  
  sd.bff <- s.em+ (S.em-s.em)*((mu.arr.bff - l.mean.em)/(L.mean.em - l.mean.em)) #dimensions of observed mu array
  sd.em <- array(data=sd.bff[No.comp:max.years.backforfill,,], dim=c(max.years,No.comp,no.surveys))
  
  
  
  obs.llike <- rep(NA, niter)
  worker.dens<-matrix(NA,ncol=No.comp,nrow=dim(Lengths_matrix)[1])
  tau.mat<-matrix(NA,ncol=No.comp,nrow=dim(Lengths_matrix)[1])
  surveyyear<-Lengths_matrix[,1:2]      ##sets up matrix for extracting parameter values specific for survey and year in TMB model
  
  for(k in 1:niter){ 
    
    print(k)
    print(k.reparam.em)
    print(L.par.vec.em)
    
    for(i in 1:dim(Lengths_matrix)[1]){
      
      worker.dens[i,] <- lambda.array[((Lengths_matrix[i,2])+1),,((Lengths_matrix[i,1])+1)] * dnorm(Lengths_matrix[i,3],mean=mu.em.array[((Lengths_matrix[i,2])+1),,((Lengths_matrix[i,1])+1)],sd=sd.em[((Lengths_matrix[i,2])+1),,((Lengths_matrix[i,1])+1)])
      
    }
    
    
    lambda.array<-ifelse(lambda.array == 0,1.1e-323,lambda.array)   ###hack that gets around log(lambda) issues in TMB when Lambda value ==0         
    
    
    
    
    
    obj_LL <- MakeADFun(
      data = list(Lengths=Lengths_matrix,
                  surveyyear=surveyyear,
                  lambda=lambda.array,
                  agezero=age1),
      
      parameters = list(log_l_mu=log(l.mean.em),
                        log_L_mu=log(L.mean.em),
                        lL=(lL.mat.em),
                        raw_sdl=(raw.sd.l.em),
                        raw_sdL=(raw.sd.L.em),
                        raw_rho=(raw.rho.em),
                        logit_k_reparam=qlogis(k.reparam.em),
                        log_s=log(s.em),
                        log_S=log(S.em)),
      
      random=c("lL"),
      
      map = list(raw_sdl=factor(NA),raw_sdL=factor(NA)),#
      
      DLL = "hier_cL_LSD_OBSLL",
      
      silent = TRUE)
    
    
    
    
    
    obs.llike[k] <- -obj_LL$fn(obj_LL$par)[1]
    
    
    if(k > 2 ){
      if(abs(obs.llike[k] - obs.llike[k-1]) <  abs(obs.llike[k-1] * rel.tolerance)){
        
        obs.opt <- nlminb(start=obj_LL$par,objective=obj_LL$fn,gradient=obj_LL$gr,silent=F)
        
        obs.rep <- sdreport(obj_LL)
        obs.srep <- summary(obs.rep)
        
        linf_overall<-(L.mean.em-(l.mean.em*(k.reparam.em^(No.comp-1))))/(1-(k.reparam.em^(No.comp -1))) 
        tzero_overall<- age1[1] - ((1/log(k.reparam.em))*log((L.mean.em - l.mean.em)/(L.mean.em - l.mean.em*k.reparam.em^(No.comp-1))))
        
        linf_cohort<-(mu.arr.bff[No.comp:max.years.comp,No.comp,1]-(mu.arr.bff[1:max.years.backforfill,1,1]*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))
        tzero_cohort<-  age1[1] -((1/log(k.reparam.em))*log((mu.arr.bff[No.comp:max.years.comp,No.comp,1]-mu.arr.bff[1:max.years.backforfill,1,1])/(mu.arr.bff[No.comp:max.years.comp,No.comp,1]-(mu.arr.bff[1:max.years.backforfill,1,1]*k.reparam.em^(No.comp-1)))))
        K<- -log(k.reparam.em)
        
        obs.llike <- obs.llike[!is.na(obs.llike)]
        return(list(obs.llike=obs.llike,Mu.obs.years=mu.em.array,Mu.all.years=mu.arr.bff,Sd.obs.years=sd.em,Sd.all.years=sd.bff,Lambda=ifelse(lambda.array==(1 / No.comp),0,lambda.array),
                    k.reparam=k.reparam.em,l=l.mean.em,L=L.mean.em,sd.l = sd.l.em,sd.L = sd.L.em, Rho= rho.em,RE.mat=lL.mat.em, 
                    l.par.vec = l.par.vec.em, L.par.vec = L.par.vec.em,sigma=sigma.em,K=K,Linf.overall =linf_overall,tzero.overall=tzero_overall, 
                    Linf.cohort = linf_cohort,tzero.cohort= tzero_cohort,Lambda.params=no.lambda.param,
                    sample.size=samplesize,age1=age1,Final.Estimate.Error = obs.srep,Entropy=EN))
        #stop("converged")
      }
    }   
    
    
    tau.mat <- worker.dens[,]/rowSums(worker.dens[,])
    
    
    
    for (i in 1:no.surveys){
      for (j in 1:no.years.covered){
        
        lambda.array.plus[j,,i] <- colSums(Lengths_matrix[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,4]*tau.mat[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,]) / (sum(Lengths_matrix[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,4]*tau.mat[((Lengths_matrix[,1])+1)==i & ((Lengths_matrix[,2])+1)==j,])) 
        
        
      }
    } 
    
    EN<- -sum(tau.mat*log(ifelse(tau.mat==0,1.1e-323,tau.mat)))
    
    obj <- MakeADFun(
      data = list(Lengths=Lengths_matrix,
                  tau=tau.mat,
                  surveyyear=surveyyear,
                  lambda=lambda.array,
                  agezero=age1),
      
      parameters = list(log_l_mu=log(l.mean.em),
                        log_L_mu=log(L.mean.em),
                        lL=(lL.mat.em),
                        raw_sdl=(raw.sd.l.em),
                        raw_sdL=(raw.sd.L.em),
                        raw_rho=(raw.rho.em),
                        logit_k_reparam=qlogis(k.reparam.em),
                        log_s=log(s.em),
                        log_S=log(S.em)),
      
      random=c("lL"),
      
      map = list(raw_sdl=factor(NA),raw_sdL=factor(NA)),#
      
      DLL = "hier_cL_LSD",silent=T)
    
    #UB<-c(Inf,Inf,log(0.05),rep(Inf,5))
    
    opt <- nlminb(start=obj$par,objective=obj$fn,gradient=obj$gr,silent=T)#,upper=UB)
    rep <- sdreport(obj)
    srep <- summary(rep)
    
    par1<-rep$par.fixed
    
    l.mean.em<-exp(par1[1])
    L.mean.em<-exp(par1[2])
    
    #raw.sd.l.em<-par1[3]         
    #raw.sd.L.em<-par1[3]
    sd.l.em<-exp(raw.sd.l.em)         
    sd.L.em<-exp(raw.sd.L.em)
    
    raw.rho.em<-par1[3]
    rho.em<- -1 + 2 * plogis(par1[3])
    
    k.reparam.em <- plogis(par1[4])
    s.em<-exp(par1[5])              
    S.em<-exp(par1[6])
    
    lL.mat.em<- matrix(srep[rownames(srep) == "lL", "Estimate"], ncol = 2)
    
    
    l.par.vec.em<- exp(lL.mat.em[,1])*l.mean.em
    L.par.vec.em<- exp(lL.mat.em[,2])*L.mean.em
    
    mu.arr.bff[1:(No.comp-1),1,1]<-sum(l.par.vec.em)/length(l.par.vec.em)
    mu.arr.bff[No.comp:max.years.backforfill,1,1]<-l.par.vec.em
    
    mu.arr.bff[((2*(No.comp-1))+1):max.years.comp,No.comp,1]<-L.par.vec.em
    mu.arr.bff[No.comp:(2*(No.comp-1)),No.comp,1]<-sum(L.par.vec.em)/length(L.par.vec.em)
    
    
    l.par.vec.plus<-mu.arr.bff[1:max.years.backforfill,1,1]
    L.par.vec.plus<-mu.arr.bff[No.comp:max.years.comp,No.comp,1]
    
    linf.em<-(L.par.vec.plus-(l.par.vec.plus*k.reparam.em^(No.comp-1)))/(1-(k.reparam.em^(No.comp-1)))## minus 2 as Final component is not restrivted to vb
    K.em <- -log(k.reparam.em)  
    tzero.em<-  age1[1]-((1/log(k.reparam.em))*log((L.par.vec.plus-l.par.vec.plus)/(L.par.vec.plus-(l.par.vec.plus*k.reparam.em^((No.comp-1)))))) 
    
    
    
    for(j in 1:no.surveys){
      
      l.var<-linf.em*(1-exp(-K.em*(age1[j]-tzero.em)))
      L.var<- linf.em*(1-exp(-K.em*((age1[j]+(No.comp-1))-tzero.em)))
      mu.arr.bff[1:max.years.backforfill,1,j]<- l.var
      mu.arr.bff[No.comp:max.years.comp,No.comp,j]<- L.var
      for(i in 2:max.years.comp){
        
        
        mu.arr.bff[i,2:No.comp,j]<-   sapply(2:No.comp,FUN=function(z){
          
          x <- i-z
          if(x>=0 && x<max.years.backforfill){
            l<-mu.arr.bff[x+1,1,j]
            L <-mu.arr.bff[x+No.comp,No.comp,j] 
          }else { 
            l<-0
            L <-0
          }
          mu.arr.bff[(i-1),(z-1),j] + ((L - l)*(((k.reparam.em^(z-2))-(k.reparam.em^(z-1)))/(1-(k.reparam.em^(No.comp-1)))))
        })
      }
    }
    
    mu.em.array <- array(data=mu.arr.bff[No.comp:max.years.backforfill,,], dim=c(max.years,No.comp,no.surveys))
    
    
    sd.bff<-s.em+ (S.em-s.em)*((mu.arr.bff - l.mean.em)/(L.mean.em - l.mean.em)) #dimensions of observed mu array
    sd.em <- array(data=sd.bff[No.comp:max.years.backforfill,,], dim=c(max.years,No.comp,no.surveys))
    
    sigma.em[1] <- s.em ###this just means that s.em and S.em will be in output
    sigma.em[2] <- S.em 
    
    lambda.array <- lambda.array.plus
    
    
  }
}



