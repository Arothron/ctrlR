#outline.fourier.demo
data(pan.section)
section=pan.section
section=cbind(1, section, 1)
section=rbind(1,section, 1)
section_junk<-imagematrix(section, type="grey")
section_junk[which(section_junk<0.95)]<-0
plot(section_junk)
cont<-Conte(round(unlist(locator(1)),0),section_junk)
lines(cont$X, cont$Y,lwd=2,col="red")
plot(NA,xlim=range(cont$X),ylim=range(cont$Y),xlab="",ylab="",bty="n",axes=FALSE)
lines(cont$X, cont$Y,lwd=2)
scalecoord<-locator(2, type="p",pch=4,lwd=2)
scalepixsize<-sqrt(diff(scalecoord$x)^2+diff(scalecoord$y)^2)
X<-cont$X*30/scalepixsize
Y<-cont$Y*30/scalepixsize
Xc <- mean(X)
Yc <- mean(Y)
plot(X,Y,type="l")
#########