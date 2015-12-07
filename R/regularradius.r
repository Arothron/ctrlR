#' regularradius
#' 
#' Function from the Morphemetric with R
#' @author Jean Claude
#' @export

regularradius<-function(Rx, Ry, n)
{le<-length(Rx)
M<-matrix(c(Rx, Ry), le, 2)
M1<-matrix(c(Rx-mean(Rx), Ry-mean(Ry)), le, 2)
V1<-complex(real=M1[,1], imaginary=M1[,2])
M2<-matrix(c(Arg(V1), Mod(V1)), le, 2)
V2<-NA
for (i in 0:(n-1))
{V2[i+1]<-which.max((cos(M2[,1]-2*i*pi/n)))}
V2<-sort(V2)
list("pixindices"=V2,"radii"=M2[V2,2],"coord"=M1[V2,])}