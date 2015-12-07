#type.landmark.demo
data("pan.model.mesh")
data("pan.landmark.type")
shade3d(pan.model.mesh,col=3)
spheres3d(pan.landmark.type,radius=3)
texts3d(pan.landmark.type,texts=c("1 type","2 type","3 type"),adj=c(1,-1),cex=2,col="red")
