#' iter.one2four
#' 
#' This repeat iteratively the function one2four
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @param iter numeric; number of iteration
#' @return prova list with 4^(iter+1) separated mesh
#' @author Antonio Profico
#' @export

iter.one2four=function(mesh,iter){
prova=one2four(mesh)  
for( m in 1:iter){  
prova=prova  
listone=list()
for(i in 1:length(prova)){
listone[[i]]=one2four(prova[[i]])
} 
biglist=list()
for(i in 1:length(listone)){
for(j in 1:4){
biglist[[length(biglist)+1]]=listone[[i]][[j]]}}
prova=biglist
}
col=rainbow(length(prova))
for(i in 1:length(prova)){
shade3d(prova[[i]],col=col[i])}  
return(prova) 
}