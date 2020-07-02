#=============
#Plot function for a individual model
#
#================

plot.monthly<- function(model,Lengths,xlimit,td){
  
  
  names(Lengths)<-c("Survey","Year","Length","RF")
  
  Lengths$Year.num<-as.numeric(factor(Lengths$Year))
  Lengths$month.num<-as.numeric(Lengths$Survey)
  
  ###set up for individual survey 
  dat <- Lengths
  year <- dat$Year
  dat <- dat[order(year),]
  yr.idx <- dat$Year.num
  m.idx <- dat$month.num
  
  
  uniq.yr <- unique(dat$Year)
  uniq.yr.num <- unique(dat$Year.num)
  ny <- length(uniq.yr)
  uniq.m <- unique(dat$Survey)
  uniq.m.num <- unique(dat$month.num)
  nm <- length(uniq.m)
  
  
  n.sample<-length(unique(interaction(dat$Year,dat$Survey)))
  
  l.vec <- dat$Length
  count <- dat$RF
  
  if(is.null(model$Mu)){#start of hier model plot fun
    
    Mu <- model$Mu.obs.years[,,,drop=F]
    Lambda <- model$Lambda[,,,drop=F]
    Sd <- model$Sd[,drop=F] #model$Sd.obs.years[,,,drop=F]
    m <- length(Mu[1,,1])  ##no of components
    
    ## predicted lengths
    maxlen<-max(Lengths$Length)
    l.pred <- seq(0, maxlen)
    brks <- seq(0, ceiling(maxlen), by = 1)
    width <- 1
    # par.store <- NA
    # 
    # for (i in 1:ny) {
    #   
    #   yr.dat <- subset(dat, Year == uniq.yr[i])
    #   l.vec <- yr.dat$Length
    #   l.all <- l.vec[rep(1:length(l.vec), times = yr.dat$RF)]
    #   par.store[i] <- max(hist(l.all, breaks = brks,plot=F)$density)
    #   }
    # 
    # 
    # maxden <- max(par.store)
    # 
    pfun <- function(){
      
      
      par(mfcol=c(1,n.sample))
      
      
      par(cex = 0.6,las=1,xaxs="i")
      par(mar = c(0, 1, 0, 0), oma = c(4, 4, 0.5, 0.5))
      par(tcl = -0.25)
      par(mgp = c(2, 0.6, 0))
      
      #xcoords <- matrix(NA,ncol=m,nrow=ny)
      #ycoords <- matrix(NA,ncol=m,nrow=ny)
      #xlimcoord<-NA
      
      for (i in 1:ny) {
        for (j in 1:nm) {
          
          
          dens.mat <- sapply(1:m, FUN = function(z){
            Lambda[uniq.yr.num[i], z,j] * dnorm(l.pred, mean = Mu[i,z,j], sd = Sd[z])
          })
          
          m.dat <- subset(dat, Year == uniq.yr[i] & Survey == uniq.m[j] )
          l.vec <- m.dat$Length
          l.all <- l.vec[rep(1:length(l.vec), times = m.dat$RF)]
          
          a<-hist(l.all, breaks = brks, col = "lightgrey", border = "grey", main = "", xlim = c(0, xlimit), probability = TRUE,axes=F,las=1,plot=F)#ylim = c(0, 0.02)
          
          title<-paste(as.character(uniq.yr[i]),"-",as.character(uniq.m[j]),sep="")
          barplot(a$density, space=0, horiz=TRUE,main=title,yaxs="i")
          if(i==1 && j==1){
            axis(2, at=(pretty(a$breaks) - a$breaks[1])/width,
                 labels=pretty(a$breaks))
          }
          matlines(y=l.pred, x=dens.mat, lty = 2,col=1, alpha=0.6)
          
          if(td==TRUE){
            lines(y=l.pred, x=rowSums(dens.mat),lty=2,col = "red")
          }
          
          #     
          #     xcoords[i,] <- grconvertX(x=Mu[,,Survey.num], "user", "ndc")
          #     ycoords[i,] <- grconvertY(y=0,"user","ndc")
          #     xlimcoord[i]<-grconvertX(x=xlimit, "user", "ndc")
          #     
          #     mtext(uniq.yr[i], side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
          #     
          #     if(i==1){
          #       if(td==T){
          #         
          #         legend("topright",legend = c("Components", "Cohort progression","Total density"),lty=c(2,1,2),col=c(1,1,2),cex=0.8)
          #       }else
          #         legend("topright",legend = c("Components", "Cohort progression"),lty=c(2,1),col=c(1,1),cex=0.8)
          #     }
          #     
          #     
          #     if(ny<8){
          #       if (i %in% c((ny/2),ny))
          #         axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
          #                 
          #       axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
          #       
          #       
          #     }else if(ny%%2>0){
          #       if (i %in% c(((ny+1)/2),ny))
          #         axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
          #       if (i %in% 1:((ny+1)/2))
          #         axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
          #       
          #     }else{
          #       
          #     if (i %in% c((ny/2),ny))
          #       axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
          #     if (i %in% 1:(ny/2))
          #       axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
          #     }#box(col = "grey60")
          #   }
          #   
          #   
          #   
          #   mtext("Length of fish (cm)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
          #         col = "grey20")
          #   mtext("Proportion of total count", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
          #         col = "grey20",las=0)
          # 
          # if(ny<8){
          #   
          #   pushViewport(viewport())
          #   for(j in 1:(ny-1)){
          #     for(i in 1:m-1){
          #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
          #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
          #       }
          #     }
          #   }
          #   popViewport()
          #   
          # }else if(ny%%2>0){
          #   
          #   pushViewport(viewport())
          #   for(j in 1:((ny+1)/2)){
          #     for(i in 1:m-1){
          #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
          #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
          #       }
          #     }
          #   }
          #   
          #   for(j in (((ny+1)/2)+1):(ny-1)){
          #     for(i in 1:m-1){
          #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
          #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
          #       }
          #     }
          #   }
          #   popViewport()
          #   
          # }else
          #   
          #   
          #   pushViewport(viewport())
          #   for(j in 1:((ny)/2)){
          #     for(i in 1:m-1){
          #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
          #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
          #       }
          #     }
          #   }
          #   
          #   for(j in (((ny)/2)+1):(ny-1)){
          #     for(i in 1:m-1){
          #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
          #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
          #       }
          #     }
          #   }
          #   popViewport()
          # 
        }
      }
      
    }
    return(pfun())
    
    
  }else
    
  ##==========================START of basic plot fun
  
  Mu <- model$Mu[,,,drop=F]
  Lambda <- model$Lambda[,,,drop=F]
  Sd <- model$Sd
  m <- length(Sd)  ##no of components
  
  ## predicted lengths
  maxlen<-max(Lengths$Length)
  l.pred <- seq(0, maxlen)
  brks <- seq(0, ceiling(maxlen), by = 1)
  width <- 1
  # par.store <- NA
  # 
  # for (i in 1:ny) {
  #   
  #   yr.dat <- subset(dat, Year == uniq.yr[i])
  #   l.vec <- yr.dat$Length
  #   l.all <- l.vec[rep(1:length(l.vec), times = yr.dat$RF)]
  #   par.store[i] <- max(hist(l.all, breaks = brks,plot=F)$density)
  #   }
  # 
  # 
  # maxden <- max(par.store)
  # 
  pfun <- function(){
    
    
    par(mfcol=c(1,n.sample))
    
    
    par(cex = 0.6,las=1,xaxs="i")
    par(mar = c(0, 1, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
    #xcoords <- matrix(NA,ncol=m,nrow=ny)
    #ycoords <- matrix(NA,ncol=m,nrow=ny)
    #xlimcoord<-NA
    
    for (i in 1:ny) {
      for (j in 1:nm) {
        
      
      dens.mat <- sapply(1:m, FUN = function(z){
        Lambda[uniq.yr.num[i], z,j] * dnorm(l.pred, mean = Mu[1,z,j], sd = Sd[z])
      })
      
      m.dat <- subset(dat, Year == uniq.yr[i] & Survey == uniq.m[j] )
      l.vec <- m.dat$Length
      l.all <- l.vec[rep(1:length(l.vec), times = m.dat$RF)]
      
      a<-hist(l.all, breaks = brks, col = "lightgrey", border = "grey", main = "", xlim = c(0, xlimit), probability = TRUE,axes=F,las=1,plot=F)#ylim = c(0, 0.02)
       
      title<-paste(as.character(uniq.yr[i]),"-",as.character(uniq.m[j]),sep="")
        barplot(a$density, space=0, horiz=TRUE,main=title,yaxs="i")
           if(i==1 && j==1){
           axis(2, at=(pretty(a$breaks) - a$breaks[1])/width,
                labels=pretty(a$breaks))
           }
      matlines(y=l.pred, x=dens.mat, lty = 2,col=1, alpha=0.6)
      
      if(td==TRUE){
        lines(y=l.pred, x=rowSums(dens.mat),lty=2,col = "red")
      }
      
  #     
  #     xcoords[i,] <- grconvertX(x=Mu[,,Survey.num], "user", "ndc")
  #     ycoords[i,] <- grconvertY(y=0,"user","ndc")
  #     xlimcoord[i]<-grconvertX(x=xlimit, "user", "ndc")
  #     
  #     mtext(uniq.yr[i], side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
  #     
  #     if(i==1){
  #       if(td==T){
  #         
  #         legend("topright",legend = c("Components", "Cohort progression","Total density"),lty=c(2,1,2),col=c(1,1,2),cex=0.8)
  #       }else
  #         legend("topright",legend = c("Components", "Cohort progression"),lty=c(2,1),col=c(1,1),cex=0.8)
  #     }
  #     
  #     
  #     if(ny<8){
  #       if (i %in% c((ny/2),ny))
  #         axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
  #                 
  #       axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
  #       
  #       
  #     }else if(ny%%2>0){
  #       if (i %in% c(((ny+1)/2),ny))
  #         axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
  #       if (i %in% 1:((ny+1)/2))
  #         axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
  #       
  #     }else{
  #       
  #     if (i %in% c((ny/2),ny))
  #       axis(1, col = "grey40", col.axis = "grey20", at = seq(0,maxlen,by=10))
  #     if (i %in% 1:(ny/2))
  #       axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(maxden,2),(round(maxden/3,2)) ))
  #     }#box(col = "grey60")
  #   }
  #   
  #   
  #   
  #   mtext("Length of fish (cm)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
  #         col = "grey20")
  #   mtext("Proportion of total count", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
  #         col = "grey20",las=0)
  # 
  # if(ny<8){
  #   
  #   pushViewport(viewport())
  #   for(j in 1:(ny-1)){
  #     for(i in 1:m-1){
  #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
  #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
  #       }
  #     }
  #   }
  #   popViewport()
  #   
  # }else if(ny%%2>0){
  #   
  #   pushViewport(viewport())
  #   for(j in 1:((ny+1)/2)){
  #     for(i in 1:m-1){
  #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
  #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
  #       }
  #     }
  #   }
  #   
  #   for(j in (((ny+1)/2)+1):(ny-1)){
  #     for(i in 1:m-1){
  #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
  #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
  #       }
  #     }
  #   }
  #   popViewport()
  #   
  # }else
  #   
  #   
  #   pushViewport(viewport())
  #   for(j in 1:((ny)/2)){
  #     for(i in 1:m-1){
  #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
  #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
  #       }
  #     }
  #   }
  #   
  #   for(j in (((ny)/2)+1):(ny-1)){
  #     for(i in 1:m-1){
  #       if(xcoords[j+1,i+1] <= xlimcoord[j]){
  #         grid.lines(x = c(xcoords[j,i],xcoords[j+1,i+1]), y = c(ycoords[j,i],ycoords[j+1,i+1]), gp = gpar(col = 1))
  #       }
  #     }
  #   }
  #   popViewport()
  # 
      }
    }
  
  }
  return(pfun())
  
  
##end of basic plot
  
  
  
  
  
  }#end of function