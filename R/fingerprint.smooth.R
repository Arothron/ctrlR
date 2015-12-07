#' fingerprint.smooth
#' 
#' This function print the algorithm smooth fingerprint on a added mesh noise
#' @param ref.mesh triangular mesh stored as object of class "mesh3d"
#' @param noise sd deviation to defin vertex shifting
#' @param alg algorithm types stored in Morpho::vcgSmooth
#' @param delta setting values for delta (for angweight and fujilaplace)
#' @param lambda setting values for lambda (for taubin algorithm)
#' @param iter number of iteration of smoothing
#' @param range number subcategory for loss and entail of anatomical information
#' @return mesh_n mesh noised
#' @author Antonio Profico
#' @export


fingerprint_smooth=function(ref.mesh,noise,alg,delta,lambda,iter,range){
alg=tolower(alg)
noi.mesh=noise.mesh(ref.mesh, noise = noise, seed = 123) 
dist_smoothed=NULL
mu=-(lambda+0.01)
for(i in 1 :iter){
if(alg=="taubin"){
for(j in 1:length(lambda)){
smoothed=vcgSmooth(noi.mesh,type=alg,it=i,lambda=lambda[j],mu=mu[j])
dist_smoothed=c(dist_smoothed,sum(abs(meshDist(ref.mesh,smoothed,plot=F)$dists)))}}
if(alg%in%c("angweight","fujilaplace")){
for(j in 1:length(delta)){
smoothed=vcgSmooth(noi.mesh,type=alg,it=i,delta=delta[j])
dist_smoothed=c(dist_smoothed,sum(abs(meshDist(ref.mesh,smoothed,plot=F)$dists)))}} 
if(alg%in%c("laplace","hclaplace")){
smoothed=vcgSmooth(noi.mesh,type=alg,it=i)
dist_smoothed=c(dist_smoothed,sum(abs(meshDist(ref.mesh,smoothed,plot=F)$dists)))}}  
thr_eu=sum(abs(meshDist(ref.mesh,noi.mesh,plot=F)$dists))
dist_matrix=NULL
if(alg=="taubin"){
dist_matrix=matrix(dist_smoothed,ncol=length(lambda),byrow=TRUE)
entail_alg=(1 - (dist_matrix/thr_eu)) * 100
rownames(entail_alg)=paste("iter",c(1:iter))
colnames(entail_alg)=paste("lambda",lambda)}
if(alg%in%c("angweight","fujilaplace")){
dist_matrix=matrix(dist_smoothed,nrow=iter,byrow=TRUE)
entail_alg=(1 - (dist_matrix/thr_eu)) * 100
rownames(entail_alg)=paste("iter",c(1:iter))
colnames(entail_alg)=paste("delta",delta)
}
if(alg%in%c("laplace","hclaplace")){
dist_matrix=matrix(dist_smoothed,ncol=1,byrow=TRUE)
entail_alg=(1 - (dist_matrix/thr_eu)) * 100
rownames(entail_alg)=paste("iter",c(1:iter))
colnames(entail_alg)=1}
block_ret=range
block_los=range
retr=seq(range(entail_alg)[2],0,length=block_ret)
loss=seq(0,range(entail_alg)[1],length=block_los)
A=colorpanel(block_ret, "blue", "grey")
B=colorpanel(block_los, "moccasin","brown")
C=c(A,B)
a=levelplot(entail_alg,border="black", 
  scales=list(cex=2,tck=0, x=list(rot=45)), 
  col.regions=colorpanel(block_ret, "blue", "grey"),
  colorkey=list(at=as.numeric(factor(c(retr,loss)[-range-1])),
  labels=as.character((round(c(retr,loss)[-range-1],2)),cex=2),col=C[-range-1]),
  alpha.regions=0.5, at=retr, main=list(label=alg,cex=3),
  xlab=NULL, ylab=NULL)
p=levelplot(entail_alg, border="black",scales=list(tck=0, x=list(rot=45)), 
  at=loss, col.regions=colorpanel(block_los, "moccasin","brown"),xlab=NULL, ylab=NULL)
return(a+p)
}
