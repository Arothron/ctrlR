#' aro.clo.points
#'
#' This function finds the vertices numbering nearest a landmark set
#' @param vertex matrix mesh vertex
#' @param landmarks numeric: a kxm matrix of a landmark set
#' @return position numeric: a vector with vertex number nearest the landmark set
#' @author Antonio Profico, Alessio Veneziano, Alessandro Lanteri
#' @examples
#' \dontrun{
#' #load an example: mesh, and L set
#' data(exp.teeth.mesh)
#' data(exp.teeth.Lset)
#' sur=t(vcgQEdecim(teeth.mesh,percent=0.50)$vb[-4, ])
#' set=teeth.Lset
#' example=aro.clo.points(vertex=sur,landmarks=set) 
#' }
#' @export
aro.clo.points=function (vertex, landmarks) 
{
    position = NULL
    for (i in 1:dim(landmarks)[1]) {
        distance = sqrt((vertex[, 1] - landmarks[i, 1])^2 + 
            (vertex[, 2] - landmarks[i, 2])^2 + (vertex[, 
            3] - landmarks[i, 3])^2)
        mindist_selector = which(distance == min(distance))
        position = c(position, mindist_selector)
    }
    return(position)
}