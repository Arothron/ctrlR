#' list.amira.dir
#'
#' This function create a list formed by files of landmark set of Amira stored in a folder
#' @param path.dir: path directory where Amira landmark set are stored
#' @param nland: number of landmark, if set on "auto" the value will be automatically set
#' @return list.amira numeric: a list with landmark coordinates (kxdxn)
#' @author Antonio Profico, Alessio Venezianp
#' @export

list.amira.dir=function(path.dir,nland){
names=list.files(path.dir)
list.amira=list()
for(i in 1:length(names)){
temp=read.amira.set(paste(path.dir,"/",names[i],sep=""),nland)
dimnames(temp)[[3]]=list(names[i])
list.amira[[i]]=temp
}
if(nland!="auto"){
for(i in 1:length(list.amira)){
dims=dim(list.amira[[i]])[1]
if(dims!=nland){
print(paste(names[i],"has a different landmark number"))}}}
names(list.amira)=names
return(list.amira)}
