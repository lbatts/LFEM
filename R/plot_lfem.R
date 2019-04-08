#=============
#Plot function for a individual model
#
#================

plot.lfem<- function(model,Lengths,Survey.num,xlimit){
  
  
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
  
  
  Mu <- model$Mu[,,,drop=F]
  Lambda <- model$Lambda[,,,drop=F]
  Sd <- model$Sd
  m <- length(Sd)  ##no of components
  
  
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
    par(cex = 0.6)
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
    xcoords <- matrix(NA,ncol=m,nrow=ny)
    ycoords <- matrix(NA,ncol=m,nrow=ny)
    xlimcoord<-NA
    
    for (i in 1:ny) {
      
      
      dens.mat <- sapply(1:m, FUN = function(z){
        Lambda[uniq.yr.num[i], z,Survey.num] * dnorm(l.pred, mean = Mu[1,z,Survey.num], sd = Sd[z])
      })
      
      yr.dat <- subset(dat, Year == uniq.yr[i])
      l.vec <- yr.dat$Length
      l.all <- l.vec[rep(1:length(l.vec), times = yr.dat$RF)]
      hist(l.all, breaks = brks, col = "lightgrey", border = "grey", main = "", xlim = c(0, xlimit), probability = TRUE,axes=F,ylim=c(0,maxden))#ylim = c(0, 0.02)
      matlines(l.pred, dens.mat, lty = 2,col=1, alpha=0.6)
      lines(l.pred, rowSums(dens.mat),lty=2,col = "red")
      
      
      
      xcoords[i,] <- grconvertX(x=Mu[,,1], "user", "ndc")
      ycoords[i,] <- grconvertY(y=0,"user","ndc")
      xlimcoord[i]<-grconvertX(x=xlimit, "user", "ndc")
      
      mtext(uniq.yr[i], side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
      
      if(i==1){
        legend("topright",legend = c("Components", "Cohort progression","Total density"),lty=c(2,1,2),col=c(1,1,2),cex=0.8)
      }
      
      
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
          col = "grey20")
  
  if(ny<8){
    
    pushViewport(viewport())
    for(j in 1:(ny-1)){
      for(i in 1:m-1){
        if(xcoords[j+1,i+1] <= xlimcoord[j]){
          grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
        }
      }
    }
    popViewport()
    
  }else if(ny%%2>0){
    
    pushViewport(viewport())
    for(j in 1:((ny+1)/2)){
      for(i in 1:m-1){
        if(xcoords[j+1,i+1] <= xlimcoord[j]){
          grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
        }
      }
    }
    
    for(j in (((ny+1)/2)+1):(ny-1)){
      for(i in 1:m-1){
        if(xcoords[j+1,i+1] <= xlimcoord[j]){
          grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
        }
      }
    }
    popViewport()
    
  }else
    
    
    pushViewport(viewport())
    for(j in 1:((ny)/2)){
      for(i in 1:m-1){
        if(xcoords[j+1,i+1] <= xlimcoord[j]){
          grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
        }
      }
    }
    
    for(j in (((ny)/2)+1):(ny-1)){
      for(i in 1:m-1){
        if(xcoords[j+1,i+1] <= xlimcoord[j]){
          grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
        }
      }
    }
    popViewport()
  
  
  
  }
  return(pfun())
  
  }#end of function