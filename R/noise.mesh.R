#' noise.mesh
#' 
#' This function add noise to a mesh
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @param noise sd deviation to defin vertex shifting
#' @param seed seed for random number generator
#' @return mesh_n mesh noised
#' @author Antonio Profico
#' @export

noise.mesh=function (mesh, noise = 0.025, seed = 123) 
{
    mesh_n = NULL
    set.seed(seed)
    noise = matrix(rnorm(3 * dim(mesh[[1]])[2], 0, sd = noise), 
        nrow = 3, ncol = dim(mesh[[1]])[2])
    mesh_n_vb = rbind(mesh$vb[1:3, ] + noise, rep(1, dim(mesh[[1]])[2]))
    mesh_n = list(vb = mesh_n_vb, it = mesh$it)
    class(mesh_n) <- "mesh3d"
    return(mesh_n)
}
