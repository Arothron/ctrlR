% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/fingerprint.smooth.R
\name{fingerprint_smooth}
\alias{fingerprint_smooth}
\title{fingerprint.smooth}
\usage{
fingerprint_smooth(ref.mesh, noise, alg, delta, lambda, iter, range)
}
\arguments{
\item{ref.mesh}{triangular mesh stored as object of class "mesh3d"}

\item{noise}{sd deviation to defin vertex shifting}

\item{alg}{algorithm types stored in Morpho::vcgSmooth}

\item{delta}{setting values for delta (for angweight and fujilaplace)}

\item{lambda}{setting values for lambda (for taubin algorithm)}

\item{iter}{number of iteration of smoothing}

\item{range}{number subcategory for loss and entail of anatomical information}
}
\value{
mesh_n mesh noised
}
\description{
This function print the algorithm smooth fingerprint on a added mesh noise
}
\author{
Antonio Profico
}

