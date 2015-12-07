#' mesh.smooth.tool
#'
#' This function find the optimal smoothing algorithm setting giving a mesh using the mesh distance as estimator. 
#' @param sur triangular mesh stored as object of class "mesh3d"
#' @param noise numeric: sd deviation to defin vertex shifting 
#' @param algorithms character: algorithm types stored in Morpho::vcgSmooth
#' @param iter numeric: number of smoothing iterations (raccomended no more than 10 iteration)
#' @param tarface numeric: target of triangle number
#' @param sel numeric vector: nteger vector selecting cleaning type
#' @param lambda numeric: setting values for lambda (if NULL automatic estimation well be done) 
#' @param tau numeric: mu value for "taubin" algorithm
#' @param deltaFJ numeric: setting values for deltaFJ (if NULL automatic estimation well be done)
#' @param deltaAW numeric: setting values for deltaAW (if NULL automatic estimation well be done)
#' @param lambda.levels numeric: length range lambda
#' @param deltaFJ.levels numeric: length range deltaFJ
#' @param deltaAW.levels numeric: length range deltaAW
#' @param lambda.start numeric: upper value lambda smoothing range
#' @param deltaFJ.start numeric: upper value deltaFJ smoothing range
#' @param deltaAW.start numeric: upper value deltaAW smoothing range
#' @param lambda.iter numeric: iteration for lambda estimation
#' @param l.lambda.iter numeric: interval for lambda.iter
#' @param deltaFJ.iter numeric: iteration for deltaFJ estimation
#' @param l.deltaFJ.iter numeric: interval for deltaFJ estimation
#' @param deltaAW.iter numeric: iteration for deltaAW estimation
#' @param l.deltaAW.iter numeric: interval for deltaAW estimation
#' @param lambda.f numeric: factor for estimation lambda
#' @param deltaFJ.f numeric: factor for estimation deltaFJ
#' @param deltaAW.f numeric: factor for estimation deltaAW
#' @param iter_scale_tau numeric: iter used in scale factor estimation for taubin algorithm
#' @param iter_scale_fuj numeric: iter used in scale factor estimation for fujilaplace algorithm
#' @param iter_scale_ang numeric: iter used in scale factor estimation for angweight algorithm
#' @return matrix matrix: matrix with all result (absolute sum of mesh distance) for smoothing setting iteration
#' @return noise.dist numeric: absolute sum of mesh distance between the model (or decimated mesh) and its noised version
#' @return tau.par numeric vector: lambda scale factor values
#' @return fuj.par numeric vector: deltaFJ scale factor values
#' @return ang.par numeric vector: deltaAW scale factor values
#' @author Antonio Profico
#' @examples
#' \dontrun{
#' #load the example 1: mesh, and L set
#' data(exp.venus.mesh)
#' run.tool=mesh.smooth.tool(sur=venus.mesh,tarface=NULL,noise=0.075)
#' run.tool.dec= mesh.smooth.tool(sur=venus.mesh,tarface=NULL, noise=0.075,lambda=run.tool$tau.par,delta_AW=run.tool$ang.par)
#' }
#' @export

mesh.smooth.tool=function(sur,
algorithms = c("taubin", "angweight","fujilaplace", "laplace", "hclaplace"),
iter = 10,tarface=10000,deltaFJ = NULL,deltaAW = NULL,lambda = NULL,tau = 0.01,
lambda.levels = 5,deltaFJ.levels = 5,deltaAW.levels = 5,lambda.start = 0.95,
deltaFJ.start = 0.95,deltaAW.start = 0.95,lambda.f = 0.99,deltaFJ.f = 0.99,
deltaAW.f = 0.99,lambda.iter = 70,l.lambda.iter = 7,deltaFJ.iter = 70,
l.deltaFJ.iter = 7,deltaAW.iter = 70,l.deltaAW.iter = 7,iter_scale_tau = 1,
iter_scale_fuj = 1,iter_scale_ang = 1,noise=0.075,sel=c(0,1,3,5,6)){
  
 ########## 
# sur=ply2mesh("willendorf.ply")
# algorithms = c("taubin", "laplace", "hclaplace")
# iter = 10;tarface=5000;deltaFJ = NULL;deltaAW = NULL;lambda = NULL;tau = 0.01;
# lambda.levels = 5;deltaFJ.levels = 5;deltaAW.levels = 5;lambda.start = 0.95;
# deltaFJ.start = 0.95;deltaAW.start = 0.95;lambda.f = 0.99;deltaFJ.f = 0.99;
# deltaAW.f = 0.99;lambda.iter = 70;l.lambda.iter = 7;deltaFJ.iter = 70;
# l.deltaFJ.iter = 7;deltaAW.iter = 70;l.deltaAW.iter = 7;iter_scale_tau = 1;
# iter_scale_fuj = 1;iter_scale_ang = 1;noise=0.075
########## 


model=vcgClean(sur,sel=sel) 
sur_n=noise.mesh(model,noise=noise)
if(is.numeric(tarface)==TRUE){
sur_n=vcgClean(vcgQEdecim(sur_n,tarface=tarface),sel=sel)}
thr_n=sum(abs(meshDist(model,sur_n,plot = F,method = "m")$dists))
if (is.null(lambda) & "taubin" %in% tolower(algorithms)) {
smoothed.t = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 5,6), silent = TRUE), type = "taubin", iter_scale_tau, 
            lambda = lambda.start, mu = (-lambda.start - tau))
s.dist.t = sum(abs(meshDist(model,smoothed.t, plot = F,method = "m")$dists))
t.lambda = lambda.start
t.lambda_junk = NULL
if (s.dist.t > thr_n) {
            p.v.dist = NULL
            iteration = round(seq(1, lambda.iter, length = l.lambda.iter))
            for (j in 1:length(iteration)) {
                t.lambda_junk = t.lambda^(iteration[j]) * lambda.f
                p.mu = -t.lambda_junk - tau
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "taubin", iter_scale_tau, 
                  lambda = t.lambda_junk, mu = p.mu)
                p.m.dist = sum(abs(meshDist(model,p.smoothed, 
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("retrieval lambda range at the first step estimated at", 
                    iteration[j], "iteration\n")) & break
                }
            }
            N1_step = length(p.v.dist)
            if (N1_step == 1) {
                iteration_2 = iteration[N1_step]:iteration[N1_step + 
                  1]
            }
            else if (N1_step == length(iteration)) {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step]
            }
            else {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step + 
                  1]
            }
            t.lambda = lambda.start
            p.mu = -t.lambda - tau
            t.lambda_junk = NULL
            p.m.dist = NULL
            p.v.dist = NULL
            for (i in 1:length(iteration_2)) {
                t.lambda_junk = t.lambda^(iteration_2[i]) * lambda.f
                p.mu = -t.lambda_junk - tau
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "taubin", iteration = iter_scale_tau, 
                  lambda = t.lambda_junk, mu = p.mu)
                p.m.dist = sum(abs(meshDist(model,p.smoothed,  
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("upper lambda parameter range estimated at", 
                    iteration_2[length(p.v.dist)], "iteration_2\n")) & 
                    break
                }
            }
        }
          p.v.dist = NULL

        if (min(p.v.dist) < thr_n) {

            lambda = seq(0, (lambda.start^iteration_2[length(p.v.dist)]) * 
                lambda.f, length = lambda.levels + 1)[-1]
        }
        if (min(p.v.dist) > thr_n) {
            lambda = seq(0, lambda.start, length = lambda.levels + 
                1)[-1]
            printf(verbose, "upper lambda parameter range not estimated, will be used the user setting for up range \n")
        }
        mu = -lambda - tau
}
    mu = -lambda - tau
if (is.null(deltaFJ) & "fujilaplace" %in% tolower(algorithms)) {
        smoothed.t = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 5, 
            6), silent = TRUE), type = "fujilaplace", iteration = iter_scale_fuj, 
            delta = deltaFJ.start)
        s.dist.t = sum(abs(meshDist(model,smoothed.t, plot = F, 
            method = "m")$dists))
        t.deltaFJ = deltaFJ.start
        if (s.dist.t > thr_n) {
            p.v.dist = NULL
            t.deltaFJ_junk = NULL
            iteration = round(seq(1, deltaFJ.iter, length = l.deltaFJ.iter))
            for (j in 1:length(iteration)) {
                t.deltaFJ_junk = t.deltaFJ^(iteration[j]) * deltaFJ.f
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "fujilaplace", 
                  iter_scale_fuj, delta = t.deltaFJ_junk)
                p.m.dist = sum(abs(meshDist(model,p.smoothed,  
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("retrieval deltaFJ estimated at the first step at", 
                    iteration[j], "iteration\n")) & break
                }
            }
            N1_step = length(p.v.dist)
            if (N1_step == 1) {
                iteration_2 = iteration[N1_step]:iteration[N1_step + 
                  1]
            }
            else if (N1_step == length(iteration)) {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step]
            }
            else {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step + 
                  1]
            }
            t.deltaFJ = deltaFJ.start
            t.deltaFJ_junk = NULL
            p.m.dist = NULL
            p.v.dist = NULL
            for (i in 1:length(iteration_2)) {
                t.deltaFJ_junk = t.deltaFJ^(iteration_2[i]) * 
                  deltaFJ.f
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "fujilaplace", 
                  iter_scale_fuj, delta = t.deltaFJ_junk)
                p.m.dist = sum(abs(meshDist(model,p.smoothed,  
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("upper deltaFJ parameter range estimated at", 
                    iteration_2[length(p.v.dist)], "permutation\n")) & 
                    break
                }
            }
        }
        if (min(p.v.dist) < thr_n) {
            deltaFJ = seq(0, (deltaFJ.start^iteration_2[length(p.v.dist)]) * 
                deltaFJ.f, length = deltaFJ.levels + 1)[-1]
        }
        if (min(p.v.dist) > thr_n) {
            deltaFJ = seq(0, deltaFJ.start, length = deltaFJ.levels + 
                1)[-1]
            printf(verbose, "upper delta (fujilaplace) parameter range not estimated, will be used the user setting for up range \n")
        }
    }

if (is.null(deltaAW) & "angweight" %in% tolower(algorithms)) {
        smoothed.t = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 5, 
            6), silent = TRUE), type = "angweight", iter_scale_ang, 
            delta = deltaAW.start)
        s.dist.t = sum(abs(meshDist(model,smoothed.t, plot = F, 
            method = "m")$dists))
        t.deltaAW = deltaAW.start
        t.deltaAW_junk = NULL
        if (s.dist.t > thr_n) {
            p.v.dist = NULL
            iteration = round(seq(1, deltaAW.iter, length = l.deltaAW.iter))
            for (j in 1:length(iteration)) {
                t.deltaAW_junk = t.deltaAW^(iteration[j]) * deltaAW.f
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "angweight", 
                  iter_scale_ang, delta = t.deltaAW_junk)
                p.m.dist = sum(abs(meshDist(model,p.smoothed, 
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("retrieval deltaAW range estimated the first step at", 
                    iteration[j], "iteration\n")) & break
                }
            }
            N1_step = length(p.v.dist)
            if (N1_step == 1) {
                iteration_2 = iteration[N1_step]:iteration[N1_step + 
                  1]
            }
            else if (N1_step == length(iteration)) {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step]
            }
            else {
                iteration_2 = iteration[N1_step - 1]:iteration[N1_step + 
                  1]
            }
            t.deltaAW = deltaAW.start
            p.m.dist = NULL
            p.v.dist = NULL
            t.deltaAW_junk = NULL
            for (i in 1:length(iteration_2)) {
                t.deltaAW_junk = t.deltaAW^(iteration_2[i]) * 
                  deltaAW.f
                p.smoothed = vcgSmooth(vcgClean(sur_n, sel = c(0:3, 
                  5, 6), silent = TRUE), type = "angweight", 
                  iter_scale_ang, delta = t.deltaAW_junk)
                p.m.dist = sum(abs(meshDist(model,p.smoothed, 
                  plot = F, method = "m")$dists))
                p.v.dist = c(p.v.dist, p.m.dist)
                if (p.m.dist < thr_n) {
                  printf(verbose, paste("upper deltaAW parameter range estimated at", 
                    iteration_2[length(p.v.dist)], "permutation\n")) & 
                    break
                }
            }
        }
        if (min(p.v.dist) < thr_n) {
            deltaAW = seq(0, (deltaAW.start^iteration_2[length(p.v.dist)]) * 
                deltaAW.f, length = deltaAW.levels + 1)[-1]
        }
        if (min(p.v.dist) > thr_n) {
            deltaAW = seq(0, deltaAW.start, length = deltaAW.levels + 
                1)[-1]
            printf(verbose, "upper delta (angweight) parameter range not estimated, will be used the user setting for up range \n")
        }
    }

num_it = rep(1, length(algorithms))
    if ("taubin" %in% tolower(algorithms) == TRUE) {
        num_it[which(tolower(algorithms) == "taubin")] = length(lambda)
    }
    if ("fujilaplace" %in% tolower(algorithms) == TRUE) {
        num_it[which(tolower(algorithms) == "fujilaplace")] = length(deltaFJ)
    }
    if ("angweight" %in% tolower(algorithms) == TRUE) {
        num_it[which(tolower(algorithms) == "angweight")] = length(deltaAW)
    }
    
smooth.lev=NULL    
smooth.sets=NULL
for (j in 1:length(algorithms)) {
        algorithm = tolower(algorithms)[j]
        it_alg = seq(1, num_it[j], 1)
        algorithm = tolower(algorithms)[j]
    for (m in 1:length(it_alg)){
    for (i in 1:iter)  {
                if (algorithm == "angweight") {
                  smoothed = vcgSmooth(sur_n, type = "angweight", 
                    iteration = i, delta = deltaAW[m])
                }
                if (algorithm == "fujilaplace") {
                  smoothed = vcgSmooth(sur_n, type = "fujilaplace", 
                    iteration = i, delta = deltaFJ[m], lambda = deltaFJ[m])
                
                }
                if (algorithm == "taubin") {
                  smoothed = vcgSmooth(sur_n, type = "taubin", iteration = i, 
                    lambda = lambda[m], mu = mu[m])
                }
                if (algorithm %in% c("laplace", "hclaplace")) {
                  smoothed = vcgSmooth(sur_n, type = algorithm, 
                    iteration = i)
                }
                if (algorithm %in% c("hclaplace", "laplace")) {
                  label_output_1 = paste(substr(algorithm,1,3), "00", i, 
                    sep = "")
                }
                if (algorithm %in% c("taubin", "fujilaplace", 
                  "angweight")) {
                  label_output_1 = paste(substr(algorithm,1,3), it_alg[m], 
                    "00", i, sep = "")
                }
                
      thr_iter=sum(abs(meshDist(model,smoothed,plot = F,method = "m")$dists))
      smooth.sets = c(smooth.sets, thr_iter)}
      smooth.lev=c(smooth.lev,label_output_1)
    }}

results=matrix(smooth.sets,nrow=iter)
colnames(results)=smooth.lev
best.pos=which(results==min(results),arr.ind=T)
best.val=results[best.pos]
best.alg=colnames(results)[best.pos[2]]
best.iter=best.pos[1]

if(substr(best.alg,1,3)=="tau"){
print(paste(substr(best.alg,1,3),lambda[as.numeric(substr(best.alg,4,4))],"iter",best.iter))}
if(substr(best.alg,1,3)=="ang"){
print(paste(substr(best.alg,1,3),deltaAW[as.numeric(substr(best.alg,4,4))],"iter",best.iter))}
if(substr(best.alg,1,3)=="fuj"){
print(paste(substr(best.alg,1,3),deltaFJ[as.numeric(substr(best.alg,4,4))],"iter",best.iter))}
if(substr(best.alg,1,3)%in%c("hcl","lap")){
print(paste(substr(best.alg,1,3),"iter",best.iter))}

out = list(matrix=results,noise.dist=thr_n,tau.par=lambda,fuj.par=deltaFJ,ang.par=deltaAW)}  