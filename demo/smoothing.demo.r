#smoothing.demo
data(exp.SCP1.mesh)
SCP1.noised=noise.mesh(SCP1.mesh)
shade3d(SCP1.noised,col=3)
mesh.dist=meshDist(SCP1.mesh,SCP1.noised)
mesh.threshold=sum(abs(mesh.dist$dists))
#hclaplace
mesh.smoo.hcl=vcgSmooth(SCP1.noised,iter=3,type="hclaplace")
mesh.hcl.dist=meshDist(SCP1.mesh,mesh.smoo.hcl)
hcl.dist=sum(abs(mesh.hcl.dist$dists))
#laplace
mesh.smoo.lap=vcgSmooth(SCP1.noised,iter=3,type="laplace")
mesh.lap.dist=meshDist(SCP1.mesh,mesh.smoo.lap)
lap.dist=sum(abs(mesh.lap.dist$dists))
#angweight
mesh.smoo.ang=vcgSmooth(SCP1.noised,iter=3,type="angweight")
mesh.ang.dist=meshDist(SCP1.mesh,mesh.smoo.ang)
ang.dist=sum(abs(mesh.ang.dist$dists))
#fujilaplace
mesh.smoo.fuj=vcgSmooth(SCP1.noised,iter=3,type="fujilaplace")
mesh.fuj.dist=meshDist(SCP1.mesh,mesh.smoo.fuj)
fuj.dist=sum(abs(mesh.fuj.dist$dists))
#taubin
mesh.smoo.tau=vcgSmooth(SCP1.noised,iter=3,type="taubin")
mesh.tau.dist=meshDist(SCP1.mesh,mesh.smoo.tau)
tau.dist=sum(abs(mesh.tau.dist$dists))
#summary
alg_dist=c(hcl.dist,lap.dist,ang.dist,fuj.dist,tau.dist)
alg_dist_percent=noquote(paste(round((1 - (alg_dist/mesh.threshold)) * 100,2),"%",sep=""))
names(alg_dist_percent)=c("hclaplace","laplace","angweight","fujilaplace","taubin")
alg_dist_percent
###taubin
Scp1_d=SCP1.mesh
scale_factor=c(0.25,0.50,1,2,4)
lambda=round(seq(0.01,0.30,length=10),2)
mu=-(lambda+0.001)
dist_smoothed_tau=NULL
thr=NULL
for (i in 1: length(scale_factor)){{
Scp1_d_s=scalemesh(Scp1_d,size=scale_factor[i])
Scp1_d_s_n=noise.mesh(Scp1_d_s, noise = (0.25/(1/scale_factor[i])), seed = 123) 
thr_eu=sum(abs(meshDist(Scp1_d_s,Scp1_d_s_n,plot=F)$dists))
thr=c(thr,thr_eu)
for (j in 1:length(lambda)){
smoothed=vcgSmooth(Scp1_d_s_n,type="tau",it=3,delta=delta[j],mu=mu[j],lambda=lambda[j])
dist_smoothed_tau=c(dist_smoothed_tau,sum(abs(meshDist(smoothed,Scp1_d_s,plot=F)$dists)))}}}
result_tau=(1-(matrix(dist_smoothed_tau,ncol=10,byrow=T)/thr))*100
rownames(result_tau)=scale_factor
colnames(result_tau)=paste("lambda",lambda)
fix(result_tau)
###fujilaplace
delta=round(seq(0.01,0.30,length=10),2)
dist_smoothed_fuj=NULL
thr=NULL
for (i in 1: length(scale_factor)){{
Scp1_d_s=scalemesh(Scp1_d,size=scale_factor[i])
Scp1_d_s_n=noise.mesh(Scp1_d_s, noise = (0.25/(1/scale_factor[i])), seed = 123) 
thr_eu=sum(abs(meshDist(Scp1_d_s,Scp1_d_s_n,plot=F)$dists))
thr=c(thr,thr_eu)
for (j in 1:length(delta)){
smoothed=vcgSmooth(Scp1_d_s_n,type="fuj",it=3,delta=delta[j])
dist_smoothed_fuj=c(dist_smoothed_fuj,sum(abs(meshDist(smoothed,Scp1_d_s,plot=F)$dists)))}}}
result_fuj=(1-(matrix(dist_smoothed_fuj,ncol=10,byrow=T)/thr))*100
rownames(result_fuj)=scale_factor
colnames(result_fuj)=paste("delta",delta)
fix(result_fuj)
###angweight
delta=round(seq(0.01,0.30,length=10),2)
dist_smoothed_ang=NULL
thr=NULL
for (i in 1: length(scale_factor)){{
Scp1_d_s=scalemesh(Scp1_d,size=scale_factor[i])
Scp1_d_s_n=noise.mesh(Scp1_d_s, noise = (0.25/(1/scale_factor[i])), seed = 123) 
thr_eu=sum(abs(meshDist(Scp1_d_s,Scp1_d_s_n,plot=F)$dists))
thr=c(thr,thr_eu)
for (j in 1:length(delta)){
smoothed=vcgSmooth(Scp1_d_s_n,type="ang",it=3,delta=delta[j])
dist_smoothed_ang=c(dist_smoothed_ang,sum(abs(meshDist(smoothed,Scp1_d_s,plot=F)$dists)))}}}
result_ang=(1-(matrix(dist_smoothed_ang,ncol=10,byrow=T)/thr))*100
rownames(result_ang)=scale_factor
colnames(result_ang)=paste("delta",delta)
fix(result_ang)
###taubin fingerprint user using settings
tau=fingerprint_smooth(ref.mesh=SCP1.mesh,noise=0.25,alg="taubin",
  lambda=round(seq(0.01,0.3,length=10),2),
  delta=NULL,iter=10,range=10)
tau
###fujilaplace fingerprint using user settings
fuj=fingerprint_smooth(ref.mesh=SCP1.mesh,noise=0.25,alg="fujilaplace",
  lambda=NULL,
  delta=round(seq(0.01,0.3,length=10),2),iter=10,range=10)
fuj
###angweight fingerprint using user settings
ang=fingerprint_smooth(ref.mesh=SCP1.mesh,noise=0.25,alg="angweigth",
  lambda=NULL,
  delta=round(seq(0.01,0.3,length=10),2),iter=10,range=10)
ang
###laplace fingerprint 
lap=fingerprint_smooth(ref.mesh=SCP1.mesh,noise=0.25,alg="laplace",
  lambda=NULL,
  delta=NULL,iter=10,range=10)
lap
###hclaplace fingerprint
hcl=fingerprint_smooth(ref.mesh=SCP1.mesh,noise=0.25,alg="hclaplace",
  lambda=NULL,
  delta=NULL,iter=10,range=10)
hcl
###example 1 smoothing tool
data(exp.dog.mesh)
data(exp.dog.SLset)
data(exp.dog.Lset)
example=aro.smooth.tool(model=dog.mesh,SL.set=dog.SLset,L.set=dog.Lset,iter=10,tarface=2000) 
###example 2 smoothing tool
data(exp.teeth.mesh)
data(exp.teeth.SLset)
data(exp.teeth.Lset)
example=aro.smooth.tool(model=teeth.mesh,SL.set=teeth.SLset,L.set=teeth.Lset,iter=10,tarface=10000,lambda.iter = 350,l.lambda.iter=20,deltaFJ.iter = 350,l.deltaFJ.iter=20,deltaAW.iter = 350,l.deltaAW.iter=20) 
###example 3 smoothing tool
data(exp.SCP1.mesh)
data(exp.SCP1.SLset)
data(exp.SCP1.Lset)
example=aro.smooth.tool(model=SCP1.mesh,SL.set=SCP1.SLset,L.set=SCP1.Lset,iter=10,tarface=10000,lambda.iter = 150,l.lambda.iter=20,deltaFJ.iter = 150,l.deltaFJ.iter=20,deltaAW.iter = 150,l.deltaAW.iter=20) 


