#digital.aligning.demo
library("ctrlR.Profico")
data(pan.model.mesh)
data(pan.model.msp)
set=read.amira.set("pan.model.msp.txt",3)
sur=pan.model.mesh
sur_half=cutMeshPlane(sur, set[1,,1], v2 = set[2,,1], v3 = set[3,,1],
normal = NULL,keep.upper = TRUE)
sur_mirr=mirror(sur_half, icpiter=10,subsample = 30)
points=aro.clo.points(t(sur_half$vb)[,-4], set[,,1])
rot_sur_mirr=rotmesh.onto(sur_mirr, t(sur_mirr$vb)[points,-4], 
   t(sur_half$vb)[points,-4])
shade3d(rot_sur_mirr$mesh,col=2)
shade3d(sur_half,col=3)
