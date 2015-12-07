#' bezier.amira.path
#'
#' This function find npoints evenly spaced from a surface path (ascii format) from Amira
#' @param path.name character: path of surface path .ascii extension file
#' @param npoint numeric: number of evenly spaced points desidered
#' @param mag numeric: how many times will be divided by the number of initial points, if NULL dec.curve will be suppress
#' @return evenly.set numeric: a kxd matrix with evenly-spaced landmark coordinates 
#' @author Antonio Profico, Alessio Veneziano
#' @export
bezier.amira.path=function(path.name,npoint,mag){
mat.bezier=read.path.amira(path.name)
if(is.null(mag)==TRUE){
evenly.set=pointsOnBezier(mat.bezier[,,1],n=npoint)$points  
} 
if(is.numeric(mag)==TRUE){
dec.mat.bezier=dec.curve(mat.bezier,mag,plot=FALSE)
rownames(dec.mat.bezier)=NULL
dec.mat.bezier=as.matrix(dec.mat.bezier)  
evenly.set=pointsOnBezier(dec.mat.bezier,n=npoint)$points}
return(evenly.set)}