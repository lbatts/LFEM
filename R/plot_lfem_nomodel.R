#=============
#Plot function for a individual model
#
#================

plot.lfem.nomodel<- function(Lengths,Survey.num,xlimit){
  
  
  names(Lengths)<-c("Survey","Year","Length","RF")
  
  Lengths$Year.num<-as.numeric(factor(Lengths$Year))
  Lengths$Survey<-as.numeric(Lengths$Survey)
  
  ###set up for individual survey 
  dat <- subset(Lengths,Lengths$Survey== Survey.num)
  year <- dat$Year
  dat <- dat[order(year),]
  yr.idx <- dat$Year.num
  
  uniq.yr <- unique(dat$Year)
  uniq.yr.num <- unique(dat$Year.num)
  ny <- length(uniq.yr)
  
  
  alpha.num <- Survey.num ###this is the number of the survey if surveys wre arranged in alphbetical order
  
  l.vec <- dat$Length
  count <- dat$RF
  
  
  #m <- length(Sd)  ##no of components
  
  
  
  
  
  ## predicted lengths
  maxlen<-max(Lengths$Length)
  l.pred <- seq(0, maxlen)
  brks <- seq(0, maxlen, by = 1)
  
  par.store <- NA
  
  for (i in 1:ny) {
    
    yr.dat <- subset(dat, Year == uniq.yr[i])
    l.vec <- yr.dat$Length
    l.all <- l.vec[rep(1:length(l.vec), times = yr.dat$RF)]
    par.store[i] <- max(hist(l.all, breaks = brks,plot=F)$density)
    }
  
  
  maxden <- max(par.store)
  
  pfun <- function(){
    
    if(ny<8){
      
      par(mfcol = c(ny, 1))
      
      
    }else if(ny%%2>0){
      par(mfcol = c((ny+1)/2, 2))
    }else{
    par(mfcol = c(ny/2, 2))
    }
    par(cex = 0.6,las=1,xaxs="i")
    par(mar = c(0, 1, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
    xlimcoord<-NA
    
    for (i in 1:ny) {
      
      
      yr.dat <- subset(dat, Year == uniq.yr[i])
      l.vec <- yr.dat$Length
      l.all <- l.vec[rep(1:length(l.vec), times = yr.dat$RF)]
      hist(l.all, breaks = brks, col = "lightgrey", border = "grey", main = "", xlim = c(0, xlimit), probability = TRUE,axes=F,las=1,ylim=c(0,maxden))#ylim = c(0, 0.02)
      
      
      
      xlimcoord[i]<-grconvertX(x=xlimit, "user", "ndc")
      
      mtext(uniq.yr[i], side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
      
     
      
      
      if(ny<8){
        if (i %in% c((ny/2),ny))
          axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
                  
        axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
        
        
      }else if(ny%%2>0){
        if (i %in% c(((ny+1)/2),ny))
          axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
        if (i %in% 1:((ny+1)/2))
          axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
        
      }else{
        
      if (i %in% c((ny/2),ny))
        axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
      if (i %in% 1:(ny/2))
        axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
      }#box(col = "grey60")
    }
    
    
    mtext("Length of fish (cm)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
          col = "grey20")
    mtext("Proportion of total count", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
          col = "grey20",las=0)
    
  
  
  
  }
  return(pfun())
  
  
##end of basic plot
  
  
  
  
  
  }#end of function