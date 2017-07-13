require(e1071) # for svm()     
require(rgl) # for 3d graphics.     
require(plot3D)
set.seed(12345)   
seed <- .Random.seed   

meanmat_pop = matrix(0,nrow=length(poplist),ncol=nrow(all),
                     dimnames=list(names(poplist),rownames(all)))
sdmat = matrix(0,nrow=length(poplist),ncol=nrow(all),
               dimnames=list(names(poplist),rownames(all)))
for(i in 1:length(poplist)){
  repmat = matrix(0,nrow=length(poplist[[i]]),ncol=nrow(poplist[[i]][[1]]))
  for(j in 1:length(poplist[[i]])){
    repmat[j,] = poplist[[i]][[j]][,ncol(poplist[[i]][[j]])]
  }
  meanmat_pop[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,sd)
}
dires <- capscale(t(meanmat_pop)~1, distance="bray")
dires_sum <- summary(dires)

t <- data.frame(x=dires_sum$sites[,1], y=dires_sum$sites[,2], z=dires_sum$sites[,3], cl=NA)
t$cl <- as.character(pclass[rownames(dires_sum$sites),1])
t$cl = as.factor(t$cl)
#t$cl[which(t$cl=="Crohn's disease")] = -1
#t$cl[which(t$cl=="Healthy")] = 1
#t = t[-which(rownames(t) %in% c("V1.CD.2","V1.UC.7","V1.CD.11")),]

cweigth = c(table(t$cl)[2]/table(t$cl)[1],1)
names(cweigth)=levels(t$cl)

svm_model <- svm(cl~x+y+z, t, type='C-classification', kernel='linear', class.weights=cweigth)
w <- t(svm_model$coefs) %*% svm_model$SV

detalization <- 100                                                                                                                                                                 
grid <- expand.grid(seq(from=min(t$x),to=max(t$x),length.out=detalization),                                                                                                         
                    seq(from=min(t$y),to=max(t$y),length.out=detalization))                                                                                                         
z <- (svm_model$rho- w[1,1]*grid[,1] - w[1,2]*grid[,2]) / w[1,3]

plot3d(grid[,1],grid[,2],z,size=1, 
       xlab=paste("PCo1",round(dires_sum$cont$importance[2,1]*100, 2),"%","explained variance"),
       ylab=paste("PCo2",round(dires_sum$cont$importance[2,2]*100, 2),"%","explained variance"), 
       zlab=paste("PCo3",round(dires_sum$cont$importance[2,3]*100, 2),"%","explained variance"))  # this will draw plane.

# adding of points to the graphics.
points3d(t$x[which(t$cl=="Crohn's disease")], t$y[which(t$cl=="Crohn's disease")], t$z[which(t$cl=="Crohn's disease")], col='firebrick2',size=20)
points3d(t$x[which(t$cl=="Healthy")], t$y[which(t$cl=="Healthy")], t$z[which(t$cl=="Healthy")], col='royalblue2',size=20)




scatter3D(t$x,t$y,t$z,type="n",bty="u",colvar=NA,shade=3,pch=16,cex=2,phi=10,theta=-120,
          col=c('firebrick2','royalblue2')[as.numeric(t$cl)],
          xlim=c(min(grid[,1]),max(grid[,1])),ylim=c(min(grid[,2]),max(grid[,2])),zlim=c(min(z),max(z)),
          xlab=paste("PCo1",round(dires_sum$cont$importance[2,1]*100, 2),"%","explained variance"),
          ylab=paste("PCo2",round(dires_sum$cont$importance[2,2]*100, 2),"%","explained variance"), 
          zlab=paste("PCo3",round(dires_sum$cont$importance[2,3]*100, 2),"%","explained variance"),
          surf=list(x=matrix(grid[,1],50,50)[1:50,],y=matrix(grid[,2],50,50)[1:50,],z=matrix(z,50,50),
                    shade=0.2, col="black", alpha=0.5, facets = T))








require(e1071) # for svm()     
require(rgl) # for 3d graphics.     
set.seed(12345)   
seed <- .Random.seed   

meanmat_sub = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
                     dimnames=list(names(sublist),rownames(alldiet)))
sdmat = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
               dimnames=list(names(sublist),rownames(alldiet)))
for(i in 1:length(sublist)){
  repmat = matrix(0,nrow=length(sublist[[i]]),ncol=nrow(sublist[[i]][[1]]))
  for(j in 1:length(sublist[[i]])){
    repmat[j,] = sublist[[i]][[j]][,ncol(sublist[[i]][[j]])]
  }
  meanmat_sub[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,sd)
}
meanmat_sub = ifelse(meanmat_sub<0,0,meanmat_sub)
dires <- capscale(meanmat_sub~1, distance="bray")
dires_sum <- summary(dires)

t <- data.frame(x=dires_sum$sites[,1], y=dires_sum$sites[,2], z=dires_sum$sites[,3], cl=NA)
t$cl <- as.character(pclass[rownames(dires_sum$sites),1])
t$cl = as.factor(t$cl)
#t = t[-which(rownames(t) %in% c("V1.CD.2","V1.UC.7","V1.CD.11")),]

cweigth = c(table(t$cl)[2]/table(t$cl)[1],1)
names(cweigth)=levels(t$cl)

svm_model <- svm(cl~x+y+z, t, type='C-classification', kernel='linear', class.weights=cweigth,scale=T)
w <- t(svm_model$coefs) %*% svm_model$SV

detalization <- 50                                                                                                                                                                 
grid <- expand.grid(seq(from=min(t$x),to=max(t$x),length.out=detalization),                                                                                                         
                    seq(from=min(t$y),to=max(t$y),length.out=detalization))                                                                                                         
z <- (svm_model$rho- w[1,1]*grid[,1] - w[1,2]*grid[,2]) / w[1,3]+1.185

scatter3D(t$x,t$y,t$z,type="n",bty="u",colvar=NA,shade=3,pch=16,cex=2,phi=20,theta=-120,
          col=c('firebrick2','royalblue2')[as.numeric(t$cl)],
          #xlim=c(min(grid[,1]),max(grid[,1])),ylim=c(min(grid[,2]),max(grid[,2])),zlim=c(min(z),max(z)),
          xlab=paste("PCo1",round(dires_sum$cont$importance[2,1]*100, 2),"%","explained variance"),
          ylab=paste("PCo2",round(dires_sum$cont$importance[2,2]*100, 2),"%","explained variance"), 
          zlab=paste("PCo3",round(dires_sum$cont$importance[2,3]*100, 2),"%","explained variance"),
          surf=list(x=matrix(grid[,1],50,50)[1:50,],y=matrix(grid[,2],50,50)[1:50,],z=matrix(z,50,50),
                    shade=0.2, col="black", alpha=0.5, facets = T))


plot3d(grid[,1],grid[,2],z,size=0.5, 
       xlab=paste("PCo1",round(dires_sum$cont$importance[2,1]*100, 2),"%","explained variance"),
       ylab=paste("PCo2",round(dires_sum$cont$importance[2,2]*100, 2),"%","explained variance"), 
       zlab=paste("PCo3",round(dires_sum$cont$importance[2,3]*100, 2),"%","explained variance"))  # this will draw plane.
# adding of points to the graphics.
points3d(t$x[which(t$cl=="Crohn's disease")], t$y[which(t$cl=="Crohn's disease")], t$z[which(t$cl=="Crohn's disease")], col='firebrick2',size=20)
points3d(t$x[which(t$cl=="Healthy")], t$y[which(t$cl=="Healthy")], t$z[which(t$cl=="Healthy")], col='royalblue2',size=20)










t <- data.frame(dires_sum$sites[,1:6], cl=NA)
colnames(t) = c("x1","x2","x3","x4","x5","x6","cl")
t$cl <- pclass[rownames(dires_sum$sites),1]
#t$cl[which(t$cl=="Crohn's disease")] = -1
#t$cl[which(t$cl=="Healthy")] = 1
t$cl = as.factor(t$cl)


cweigth = c(4,1)
names(cweigth)=levels(t$cl)

svm_model <- svm(cl~x1+x2, t, type='C-classification', kernel='linear',scale=FALSE,class.weights=cweigth,tolerance=0.00000000001)
w <- t(svm_model$coefs) %*% svm_model$SV

plot(t$x1,t$x2,col=t$cl)
abline(coef=w[1,1:2])



detalization <- 100                                                                                                                                                                 
grid <- expand.grid(seq(from=min(t$x),to=max(t$x),length.out=detalization),                                                                                                         
                    seq(from=min(t$y),to=max(t$y),length.out=detalization))                                                                                                         
z <- (svm_model$rho- w[1,1]*grid[,1] - w[1,2]*grid[,2]) / w[1,3]

plot3d(grid[,1],grid[,2],z)  # this will draw plane.
# adding of points to the graphics.
points3d(t$x[which(t$cl==-1)], t$y[which(t$cl==-1)], t$z[which(t$cl==-1)], col='red')
points3d(t$x[which(t$cl==1)], t$y[which(t$cl==1)], t$z[which(t$cl==1)], col='blue')