#' imagematrix
#' 
#' Function from the removed package rimage
#' @author Tomomi TAKASHINA
#' @export

imagematrix <- function(mat, type=NULL, ncol=dim(mat)[1], nrow=dim(mat)[2],
                        noclipping=FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
  if (length(imgdim) == 3 && type == "grey") {
    mat <- rgb2grey(mat)
  }
  if (noclipping == FALSE && ((min(mat) < 0) || (1 < max(mat)))) {
    warning("Pixel values were automatically clipped because of range over.") 
    mat <- clipping(mat)
  }
  mat <- array(mat, dim=imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}

print.imagematrix <- function(x, ...) {
  x.dim <- dim(x)
  cat("size: ", x.dim[1], "x", x.dim[2], "\n")
  cat("type: ", attr(x, "type"), "\n")
}

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

