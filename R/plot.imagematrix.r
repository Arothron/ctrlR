#' plot.imagematrix
#' 
#' Function from the Morphemetric with R
#' @author Jean Claude
#' @export

plot.imagematrix <- function(x, ...) {
  colvec <- switch(attr(x, "type"),
                grey=grey(x),
                rgb=rgb(x[,,1], x[,,2], x[,,3]))
  if (is.null(colvec)) stop("image matrix is broken.")
  colors <- unique(colvec)
  colmat <- array(match(colvec, colors), dim=dim(x)[1:2])
  image(x = 0:(dim(colmat)[2]), y=0:(dim(colmat)[1]),
        z = t(colmat[nrow(colmat):1, ]), col=colors,
        xlab="", ylab="", axes=FALSE, asp=1, ...)
}

imageType <- function(x) {
  attr(x, "type")
}