#' bezier.amira.set
#'
#' This function find npoints evenly spaced from a landmark file set (.landmarkAscii)
#' @param name.file character: path of landmark file set (.landmarkAscii)
#' @param npoint numeric: number of evenly spaced points desidered
#' @param mag numeric: how many times will be divided by the number of initial points, if NULL dec.curve will be suppress
#' @return evenly.set numeric: matrix with evenly-spaced landmark coordinates (kxd)
#' @author Antonio Profico, Alessio Veneziano
#' @export
bezier.amira.set=function(name.file,nland,npoint,mag=3){
mat.bezier=read.amira.set(name.file,nland)
if(is.null(mag)==TRUE){
evenly.set=pointsOnBezier(mat.bezier[,,1],n=npoint)$points  
} 
if(is.numeric(mag)==TRUE){dec.mat.bezier=dec.curve(mat.bezier[,,1],mag,plot=FALSE)
evenly.set=pointsOnBezier(dec.mat.bezier,n=npoint)$points}
return(evenly.set)}