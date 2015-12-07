#' one2four
#' 
#' This function divide a mesh in 4 parts
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @param logical: if TRUE plot the final 4 separeted mesh
#' @return storage.mesh list with 4 separated mesh
#' @author Antonio Profico
#' @export

one2four=function(mesh,plot=FALSE){
storage.mesh=list()
bounding=meshcube(mesh)  
centro=c(mean(bounding[,1]),mean(bounding[,2]),mean(bounding[,3]))
centro1=c(mean(bounding[c(1,3,5,7),1]),mean(bounding[c(1,3,5,7),2]),mean(bounding[c(1,3,5,7),3]))
centro2=c(mean(bounding[c(2,4,6,8),1]),mean(bounding[c(2,4,6,8),2]),mean(bounding[c(2,4,6,8),3]))
centro3=c(mean(bounding[c(5,6,7,8),1]),mean(bounding[c(5,6,7,8),2]),mean(bounding[c(5,6,7,8),3]))
centro4=c(mean(bounding[c(1,2,5,6),1]),mean(bounding[c(1,2,5,6),2]),mean(bounding[c(1,2,5,6),3]))
centro5=c(mean(bounding[c(3,4,7,8),1]),mean(bounding[c(3,4,7,8),2]),mean(bounding[c(3,4,7,8),3]))
centro6=c(mean(bounding[c(1,2,3,4),1]),mean(bounding[c(1,2,3,4),2]),mean(bounding[c(1,2,3,4),3]))
newb=rbind(bounding,centro,centro1,centro2,centro3,centro4,centro5,centro6)
orig=newb
uno=cutMeshPlane(mesh,as.numeric(newb[13,]), as.numeric(newb[12,]),
  as.numeric(newb[14,]), normal = NULL,keep.upper = TRUE)  
bounding=meshcube(uno)
centro=c(mean(bounding[,1]),mean(bounding[,2]),mean(bounding[,3]))
centro1=c(mean(bounding[c(1,3,5,7),1]),mean(bounding[c(1,3,5,7),2]),mean(bounding[c(1,3,5,7),3]))
centro2=c(mean(bounding[c(2,4,6,8),1]),mean(bounding[c(2,4,6,8),2]),mean(bounding[c(2,4,6,8),3]))
centro3=c(mean(bounding[c(5,6,7,8),1]),mean(bounding[c(5,6,7,8),2]),mean(bounding[c(5,6,7,8),3]))
centro4=c(mean(bounding[c(1,2,5,6),1]),mean(bounding[c(1,2,5,6),2]),mean(bounding[c(1,2,5,6),3]))
centro5=c(mean(bounding[c(3,4,7,8),1]),mean(bounding[c(3,4,7,8),2]),mean(bounding[c(3,4,7,8),3]))
centro6=c(mean(bounding[c(1,2,3,4),1]),mean(bounding[c(1,2,3,4),2]),mean(bounding[c(1,2,3,4),3]))
newb=rbind(bounding,centro,centro1,centro2,centro3,centro4,centro5,centro6)
due=cutMeshPlane(uno,as.numeric(newb[11,]), as.numeric(newb[13,]),
  as.numeric(newb[10,]), normal = NULL,keep.upper = TRUE)
uno_ok=cutMeshPlane(uno,as.numeric(newb[15,]), as.numeric(newb[10,]),
  as.numeric(newb[11,]), normal = NULL,keep.upper = TRUE)
due_ok=cutMeshPlane(uno,as.numeric(newb[15,]), as.numeric(newb[10,]),
  as.numeric(newb[11,]), normal = NULL,keep.upper = FALSE)  
tre=cutMeshPlane(mesh,as.numeric(orig[13,]), as.numeric(orig[15,]),
  as.numeric(orig[14,]), normal = NULL,keep.upper = TRUE)
bounding=meshcube(tre)
centro=c(mean(bounding[,1]),mean(bounding[,2]),mean(bounding[,3]))
centro1=c(mean(bounding[c(1,3,5,7),1]),mean(bounding[c(1,3,5,7),2]),mean(bounding[c(1,3,5,7),3]))
centro2=c(mean(bounding[c(2,4,6,8),1]),mean(bounding[c(2,4,6,8),2]),mean(bounding[c(2,4,6,8),3]))
centro3=c(mean(bounding[c(5,6,7,8),1]),mean(bounding[c(5,6,7,8),2]),mean(bounding[c(5,6,7,8),3]))
centro4=c(mean(bounding[c(1,2,5,6),1]),mean(bounding[c(1,2,5,6),2]),mean(bounding[c(1,2,5,6),3]))
centro5=c(mean(bounding[c(3,4,7,8),1]),mean(bounding[c(3,4,7,8),2]),mean(bounding[c(3,4,7,8),3]))
centro6=c(mean(bounding[c(1,2,3,4),1]),mean(bounding[c(1,2,3,4),2]),mean(bounding[c(1,2,3,4),3]))
newb=rbind(bounding,centro,centro1,centro2,centro3,centro4,centro5,centro6)
due=cutMeshPlane(uno,as.numeric(newb[11,]), as.numeric(newb[13,]),
  as.numeric(newb[10,]), normal = NULL,keep.upper = TRUE)
tre_ok=cutMeshPlane(tre,as.numeric(newb[15,]), as.numeric(newb[10,]),
  as.numeric(newb[11,]), normal = NULL,keep.upper = TRUE)
qua_ok=cutMeshPlane(tre,as.numeric(newb[15,]), as.numeric(newb[10,]),
  as.numeric(newb[11,]), normal = NULL,keep.upper = FALSE)
  storage.mesh[[1]]=uno_ok
  storage.mesh[[2]]=due_ok
  storage.mesh[[3]]=tre_ok
  storage.mesh[[4]]=qua_ok
  if(plot=="TRUE"){
  open3d()
  for(i in 1:4){
  shade3d(storage.mesh[[i]],col=i+1)}}
  return(storage.mesh)
}