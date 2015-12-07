#' order.open.outline
#'
#' This function order an outline along the xz plane and eventually convert it in 2 dimensions
#' @param matrix numeric: matrix of outline coordinates (3D)
#' @param to.2D logical: if TRUE the ordered matrix will be converted in 2D
#' @return matt numeric: matrix with ordered points 
#' @author Antonio Profico, Alessio Veneziano
#' @export
order.open.outline=function(matrix,to.2D=FALSE){
matt_junk=matrix
RIF=cbind(seq(0,bezierArcLength(matt_junk)$arc.length,length=dim(matt_junk)[1]),rep(0,dim(matt_junk)[1]),rep(0,dim(matt_junk)[1]))
RIF2=pcAlign(as.matrix(matt_junk),RIF)
matrix_x=(RIF2[,1]-RIF2[1,1])
matrix_y=(RIF2[,2]-RIF2[1,2])
matrix_z=(RIF2[,3]-RIF2[1,3])
matrix_xyz=cbind(matrix_x,matrix_y,matrix_z)
if(matrix_xyz[round(dim(matrix_xyz)[1]/2),1]<0){
matrix_xyz[,1]=matrix_xyz[,1]*-1
matt=matrix_xyz}  
if(matrix_xyz[round(dim(matrix_xyz)[1]/2),3]<0){
matrix_xyz[,3]=matrix_xyz[,3]*-1
matt=matrix_xyz}
if(to.2D==TRUE){
matrix_xyz[,2]=0  
matt=cbind(matrix_xyz[,1],matrix_xyz[,3])  
}
return(matt)}