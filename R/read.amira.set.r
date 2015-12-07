#' read.amira.set
#'
#' This function convert a file landmark set of Amira in a array
#' @param name.file character: path of landamrkascii file
#' @param nland numeric: number of landmark sampled in Amira, if is set on "auto" will be automatically recognized the landmark number
#' @return array.set numeric: an array with landmark coordinates (kxdxn)
#' @author Antonio Profico, Alessio Veneziano
#' @export
#' 

read.amira.set=function(name.file,nland){
A <- readLines(name.file, n = 100)
end <- which(A == "@1")
end_2 <- which(A == "@2")
if(length(end_2!=0)){
print(paste("file named",name.file, "contains a second matrix of 000"))
B_junk=read.table(name.file,skip=end,nrows=(end_2-end-2))
if(nland!="auto"){
if (dim(B_junk)[1] != nland){
print(paste("nland is different from dim(matrix)[1]: ",paste("nland=",nland,",",sep=""),
  paste("dim(matrix)[1]=", dim(B_junk)[1],sep="")))  
B=matrix(NA,ncol=3,nrow=nland)}
if(dim(B_junk)[1]==nland){
B=B_junk  
}}
if(nland=="auto"){
B=read.table(name.file,skip=end,nrows=(end_2-end-2))
}}
if(length(end_2)==0){
B_junk=read.table(name.file,skip=end)
if(nland!="auto"){
if (dim(B_junk)[1] != nland){
print(paste("nland is different from dim(matrix)[1]: ",paste("nland=",nland,",",sep=""),
paste("dim(matrix)[1]=", dim(B_junk)[1],sep="")))
B=matrix(NA,ncol=3,nrow=nland)}
if(dim(B_junk)[1]==nland){
B=B_junk  
}}
if(nland=="auto"){
B=B_junk}}
array.set=array(as.matrix(B),dim=c(dim(B)[1],3,1))
dimnames(array.set)[[3]]=list(name.file)
return(array.set)}
