##-----------------------------
## Warpper function for hierarchical LFEM model functions
## Author: LB 
## Date: 08/04/2019
## This version allows multiple surveys from diffrent times of year
## Note L and l inputs should be approximate for the survey which comes first alphabetically
## 
##-----------------------------
hier.LFEM<-function(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp,SD.type,RE.type, fix.RESD, Lengths,niter,rel.tolerance,dllroot,sdl,sdL,sdk)    {
  
  if(fix.RESD==F){
    
    if(SD.type==3){
    
      if(RE.type==1){ source(paste(dllroot,"R/hier_cL_LSD_subfun.R",sep=""))
      tmp<-hier_linearSD_L(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
      if(RE.type==2){ source(paste(dllroot,"R/hier_ck_LSD_subfun.R",sep=""))
      tmp<-hier_linearSD_k(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
      if(RE.type==3){source(paste(dllroot,"R/hier_yk_LSD_subfun.R",sep=""))
      tmp<-hier_linearSD_yk(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
    
      
      }else{
      
      if(RE.type==1) {source(paste(dllroot,"R/hier_cL_CSD_subfun.R",sep=""))
      tmp<-hier_constantSD_L(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
      if(RE.type==2) {source(paste(dllroot,"R/hier_ck_CSD_subfun.R",sep=""))
      tmp<-hier_constantSD_k(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
      if(RE.type==3){source(paste(dllroot,"R/hier_yk_CSD_subfun.R",sep="")) 
      tmp<-hier_constantSD_yk(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance)
      return(tmp)}
      
    }
  }else{####start of fixRESD
    
    if(SD.type==3){#linear SD
      
      if(RE.type==1){ 
      source(paste(dllroot,"R/hier_cL_LSD_subfun_FRESD.R",sep=""))
        tmp<-hier_linearSD_L_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdL)
        return(tmp)}
      if(RE.type==2){ 
         source(paste(dllroot,"R/hier_ck_LSD_subfun_FRESD.R",sep=""))
        tmp<-hier_linearSD_k_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdk)
        return(tmp)}
      if(RE.type==3){
        source(paste(dllroot,"R/hier_yk_LSD_subfun_FRESD.R",sep=""))
        tmp<-hier_linearSD_yk_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdk)
        return(tmp)}
      
      
    }else{#constant SD
      if(RE.type==1){ 
        source(paste(dllroot,"R/hier_cL_CSD_subfun_FRESD.R",sep=""))
        tmp<-hier_constantSD_L_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdL)
        return(tmp)}
      if(RE.type==2){ 
        source(paste(dllroot,"R/hier_ck_CSD_subfun_FRESD.R",sep=""))
        tmp<-hier_constantSD_k_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdk)
        return(tmp)}
      if(RE.type==3){
        source(paste(dllroot,"R/hier_yk_CSD_subfun_FRESD.R",sep=""))
        tmp<-hier_constantSD_yk_FRESD(year0,no.years, age1, L, l, k.reparam, sigma.start,No.comp, Lengths,niter,rel.tolerance,sdl,sdk)
        return(tmp)}
  }
}
  
  

  
}#end of wrapper function

  