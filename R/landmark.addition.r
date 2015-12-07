#' landmark.addition
#' 
#' Function from the Morphemetric with R
#' @author Jean Claude
#' @export

landmark.addition<-function(M, n)
{a<-0
while(a<=n)
{p<-dim(M)[1]
k<-dim(M)[2]
N<-matrix (NA,2*p,k)
N[((1:p)*2)-1,]<-M
N[(1:p)*2,]<-(M+(rbind(M[-1,],M[1,])))/2
M<-N}
M}