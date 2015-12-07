#' export_amira
#'
#' This function export a list of 3D landmark set in separate files (format landmarkAscii)
#' @param lista: list of 3D landmark sets
#' @param path character: path of the folder where saving the Amira landmark sets
#' @author Antonio Profico, Alessio Veneziano
#' @export


export_amira<-function(lista,path=getwd()){
if (is.null(names(lista))==TRUE) {
nomi=paste("set",1:length(lista))
}else{NULL}
  
if (is.null(names(lista))==FALSE){
nomi=names(lista)
}else{NULL}
   
for (i in 1:length(lista)){
cat(paste("# AmiraMesh 3D ASCII 2.0", "\r\n", sep = ""), file =paste(path,"/",nomi[i],".txt",sep=""),append = TRUE, sep = "")
cat(paste("define Markers ",dim(lista[[i]])[1] , "\r\n", sep = ""), file = paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste("Parameters {", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste('NumSets 1,', "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste('ContentType "LandmarkSet",', "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste(" color 0 0 1", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste("}", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste("Markers { float[3] Coordinates } @1", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste("# Data section follows", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
cat(paste("@1", "\r\n", sep = ""), file= paste(path,"/",nomi[i],".txt",sep=""), 
            append = TRUE, sep = "")
  
write.table(format(lista[[i]], scientific = F, trim = T),file= paste(path,"/",nomi[i],".txt",sep=""),
  sep = " ", append = TRUE, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, na = "")  
}
}
