#' aro.smooth.tool
#'
#' This function find the optimal smoothing algorithm setting giving a mesh a landmark set and a semi-landmark set. 
#' @param model triangular mesh stored as object of class "mesh3d"
#' @param SL.set character: kxm matrix semi-landmark set 
#' @param L.set character: kxm matrix landmark set 
#' @param algorithms character: algorithm types stored in Morpho::vcgSmooth
#' @param iter numeric: number of smoothing iterations (raccomended no more than 10 iteration)
#' @param tarface numeric: target of triangle number
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
#' @param num.cores numeric: number of CPUs cores, if NULL use the number of physical CPUs/cores
#' @param iter_scale_tau numeric: iter used in scale factor estimation for taubin algorithm
#' @param iter_scale_fuj numeric: iter used in scale factor estimation for fujilaplace algorithm
#' @param iter_scale_ang numeric: iter used in scale factor estimation for angweight algorithm
#' @return matrix_result matrix: matrix with all result (loss/retrieval of anatomical information) for smoothing setting iteration
#' @return CS numeric vector: semi-landmark centroid size values in starting, decimated, and best smoothed surfaces
#' @return L_set numeric: matrix for landmark set of best smoothed surface
#' @return SL_set numeric: matrix for landmark set of best smoothed surface
#' @return mesh best smoothed mesh
#' @author Antonio Profico, Alessio Veneziano, Alessandro Lanteri, Paolo Piras
#' @examples
#' \dontrun{
#' #load the example 1: mesh, and L set
#' data(exp.dog.mesh)
#' data(exp.dog.SLset)
#' data(exp.dog.Lset)
#' example=aro.smooth.tool(model=dog.mesh,SL.set=dog.SLset,L.set=dog.Lset,iter=10,tarface=2000) 
#' }
#' 
#' \dontrun{
#' load the example 2: mesh, and L set
#' data(exp.teeth.mesh)
#' data(exp.teeth.SLset)
#' data(exp.teeth.Lset)
#' example=aro.smooth.tool(model=teeth.mesh,SL.set=teeth.SLset,L.set=teeth.Lset,iter=10,tarface=10000,lambda.iter = 350,l.lambda.iter=20,deltaFJ.iter = 350,l.deltaFJ.iter=20,deltaAW.iter = 350,l.deltaAW.iter=20) 
#' }
#' \dontrun{
#' load the example 3: mesh, and L set
#' data(exp.SCP1.mesh)
#' data(exp.SCP1.SLset)
#' data(exp.SCP1.Lset)
#' example=aro.smooth.tool(model=SCP1.mesh,SL.set=SCP1.SLset,L.set=SCP1.Lset,iter=10,tarface=10000,lambda.iter = 150,l.lambda.iter=20,deltaFJ.iter = 150,l.deltaFJ.iter=20,deltaAW.iter = 150,l.deltaAW.iter=20) 
#' }
#' @export

aro.smooth.tool=function(
model,
SL.set,
L.set,
algorithms = c("taubin","angweight","fujilaplace","laplace","hclaplace"),
iter=10,
tarface,
deltaFJ = NULL,
deltaAW =  NULL,
lambda =c(0.01,0.05,0.10,0.15,0.20),
tau =0.01,
lambda.levels = 5,
deltaFJ.levels = 5,
deltaAW.levels = 5,
lambda.start = 0.95,
deltaFJ.start = 0.95,
deltaAW.start = 0.95,
lambda.f=0.99,
deltaFJ.f=0.99,
deltaAW.f=0.99,
lambda.iter = 70,
l.lambda.iter=7,
deltaFJ.iter = 70,
l.deltaFJ.iter=7,
deltaAW.iter = 70,
l.deltaAW.iter=7,
num.cores = NULL,
iter_scale_tau=1,
iter_scale_fuj=1,
iter_scale_ang=1){
  



real = model
dicom = vcgClean(vcgQEdecim(real, tarface = tarface, quality = F,     
                 qthresh = 0.5,silent = TRUE), sel = c(0:3,5,6),silent = TRUE) 
set = as.matrix(L.set) 
land_prj_dicom = t(projRead(set, dicom)$vb[-4, ])
position = aro.clo.points(t(dicom$vb[-4, ]), land_prj_dicom)
  if(length(unique(position))!=dim(land_prj_dicom)[1]){
    stop("landmark too close: increase tarface values or exclude landmark too close")}
    if (TRUE %in% unlist(lapply(list(lambda, deltaFJ, deltaAW),is.null)) == TRUE) {
    m.dist.t = sum(abs(meshDist(dicom, real, plot = F,method="m")$dists))
  }
  if (is.null(lambda) & "taubin" %in% tolower(algorithms)) {
    smoothed.t = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent = TRUE), type = "taubin", 
                           iter_scale_tau, lambda = lambda.start, mu = (-lambda.start-tau)) 
    s.dist.t = sum(abs(meshDist(smoothed.t, real, plot = F,method="m")$dists)) 
    t.lambda = lambda.start 
    t.lambda_junk=NULL
    if (s.dist.t > m.dist.t) {
      p.v.dist = NULL
      iteration=round(seq(1,lambda.iter,length=l.lambda.iter)) 
      for(j in 1:length(iteration)){
        t.lambda_junk = t.lambda^(iteration[j])*lambda.f 
        p.mu = -t.lambda_junk-tau 
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                     type = "taubin", iter_scale_tau, lambda = t.lambda_junk, mu = p.mu) 
        p.m.dist = sum(abs(meshDist(p.smoothed, real, plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
        if (p.m.dist < m.dist.t) {
          printf(verbose,paste("retrieval lambda range at the first step estimated at",iteration[j], "iteration\n")) & break}
      }
      N1_step=length(p.v.dist)
      if(N1_step==1){
        iteration_2=iteration[N1_step]:iteration[N1_step+1]  
      } else if(N1_step==length(iteration)){
        iteration_2=iteration[N1_step-1]:iteration[N1_step]
      } else {iteration_2=iteration[N1_step-1]:iteration[N1_step+1]}
      t.lambda=lambda.start
      p.mu = - t.lambda- tau        
      t.lambda_junk=NULL
      p.m.dist=NULL
      p.v.dist=NULL
      for (i in 1:length(iteration_2)){
        t.lambda_junk = t.lambda^(iteration_2[i])*lambda.f
        p.mu = -t.lambda_junk - tau
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE),type = "taubin", iteration=iter_scale_tau,
                               lambda = t.lambda_junk, mu = p.mu) 
        p.m.dist=sum(abs(meshDist(p.smoothed, real,plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
        if (p.m.dist < m.dist.t){ 
                   printf(verbose,paste("upper lambda parameter range estimated at", 
                               iteration_2[length(p.v.dist)], "iteration_2\n")) & break}}}
    
    if (min(p.v.dist) < m.dist.t){
      lambda = seq(0, (lambda.start^iteration_2[length(p.v.dist)])*lambda.f,length = lambda.levels + 1)[-1]}
    if (min(p.v.dist) > m.dist.t) {
      lambda = seq(0, lambda.start,length = lambda.levels + 1)[-1]
      printf(verbose,"upper lambda parameter range not estimated, will be used the user setting for up range \n")}
    mu = - lambda- tau
  }else{mu=- lambda- tau}

if (is.null(deltaFJ) & "fujilaplace" %in% tolower(algorithms)) {
    smoothed.t = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent = TRUE), type = "fujilaplace", 
                           iteration=iter_scale_fuj, delta = deltaFJ.start)
    s.dist.t = sum(abs(meshDist(smoothed.t, real, plot = F,method="m")$dists))
    t.deltaFJ = deltaFJ.start
    if (s.dist.t > m.dist.t){
      p.v.dist = NULL
      t.deltaFJ_junk=NULL
      iteration=round(seq(1,deltaFJ.iter,length=l.deltaFJ.iter))
      for(j in 1:length(iteration)){
        t.deltaFJ_junk =t.deltaFJ^(iteration[j])*deltaFJ.f
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "fujilaplace", iter_scale_fuj, delta = t.deltaFJ_junk)
        p.m.dist = sum(abs(meshDist(p.smoothed, real, 
                       plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
            if (p.m.dist < m.dist.t){
            printf(verbose,paste("retrieval deltaFJ estimated at the first step at", 
                               iteration[j], "iteration\n")) & break}
      }
      N1_step=length(p.v.dist)
      if(N1_step==1){
        iteration_2=iteration[N1_step]:iteration[N1_step+1]  
      } else if(N1_step==length(iteration)){
        iteration_2=iteration[N1_step-1]:iteration[N1_step]
      } else {iteration_2=iteration[N1_step-1]:iteration[N1_step+1]}
      t.deltaFJ=deltaFJ.start
      t.deltaFJ_junk=NULL
      p.m.dist=NULL
      p.v.dist=NULL
      for (i in 1:length(iteration_2)){
        t.deltaFJ_junk =t.deltaFJ^(iteration_2[i])*deltaFJ.f
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "fujilaplace", iter_scale_fuj, delta = t.deltaFJ_junk) 
        p.m.dist=sum(abs(meshDist(p.smoothed, real, 
                                  plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
          if (p.m.dist < m.dist.t){ 
          printf(verbose,paste("upper deltaFJ parameter range estimated at", 
                               iteration_2[length(p.v.dist)], "permutation\n")) & break
      }}}
    
    if (min(p.v.dist) < m.dist.t){
      deltaFJ = seq(0,  (deltaFJ.start^iteration_2[length(p.v.dist)])*deltaFJ.f, 
                    length = deltaFJ.levels + 1)[-1]}
    if (min(p.v.dist) > m.dist.t) {
      deltaFJ = seq(0,  deltaFJ.start,length = deltaFJ.levels + 1)[-1]
      printf(verbose,"upper delta (fujilaplace) parameter range not estimated, will be used the user setting for up range \n")}
  }
if (is.null(deltaAW)& "angweight" %in% tolower(algorithms)) {
    smoothed.t = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent = TRUE), type = "angweight", 
                           iter_scale_ang, delta = deltaAW.start)
    s.dist.t = sum(abs(meshDist(smoothed.t, real, plot = F,method="m")$dists))
    t.deltaAW = deltaAW.start
    t.deltaAW_junk=NULL
    if (s.dist.t > m.dist.t) {
      p.v.dist = NULL
      iteration=round(seq(1,deltaAW.iter,length=l.deltaAW.iter))
      for(j in 1:length(iteration)){
        t.deltaAW_junk = t.deltaAW^(iteration[j])*deltaAW.f
        
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "angweight", iter_scale_ang, delta = t.deltaAW_junk)
        p.m.dist = sum(abs(meshDist(p.smoothed, real, 
                                    plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
        if (p.m.dist < m.dist.t){
          printf(verbose,paste("retrieval deltaAW range estimated the first step at", 
                               iteration[j], "iteration\n")) & break}
      }
      N1_step=length(p.v.dist)
      if(N1_step==1){
        iteration_2=iteration[N1_step]:iteration[N1_step+1]  
      } else if(N1_step==length(iteration)){
        iteration_2=iteration[N1_step-1]:iteration[N1_step]
      } else {iteration_2=iteration[N1_step-1]:iteration[N1_step+1]}
      
      t.deltaAW=deltaAW.start
      p.m.dist=NULL
      p.v.dist=NULL
      t.deltaAW_junk=NULL
      for (i in 1:length(iteration_2)){
        t.deltaAW_junk = t.deltaAW^(iteration_2[i])*deltaAW.f
        p.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "angweight", iter_scale_ang, delta = t.deltaAW_junk) 
        p.m.dist=sum(abs(meshDist(p.smoothed, real, 
                                  plot = F,method="m")$dists))
        p.v.dist = c(p.v.dist, p.m.dist)
        if (p.m.dist < m.dist.t){
          printf(verbose,paste("upper deltaAW parameter range estimated at", 
                               iteration_2[length(p.v.dist)], "permutation\n")) & break}
      }}
    
    if (min(p.v.dist) < m.dist.t) {
      deltaAW = seq(0,  (deltaAW.start^iteration_2[length(p.v.dist)])*deltaAW.f, 
                    length = deltaAW.levels + 1)[-1]}
    if (min(p.v.dist) > m.dist.t) {
      deltaAW = seq(0,  deltaAW.start, 
                    length = deltaAW.levels + 1)[-1]
      printf(verbose,"upper delta (angweight) parameter range not estimated, will be used the user setting for up range \n")}
  }
num_it = rep(1, length(algorithms))
  if ("taubin" %in% tolower(algorithms) == TRUE) {
    num_it[which(tolower(algorithms) == "taubin")] = length(lambda)}
  if ("fujilaplace" %in% tolower(algorithms) == TRUE) {
    num_it[which(tolower(algorithms) == "fujilaplace")] = length(deltaFJ)}
  if ("angweight" %in% tolower(algorithms) == TRUE) {
    num_it[which(tolower(algorithms) == "angweight")] = length(deltaAW)}
  if (is.null(num.cores)){
    num.cores = detectCores()}
dir.create("sur.tmp") 
nland = dim(set)[1] 
set_dicom = t(dicom$vb[-4, ])[position, ] 
land_prj = projRead(set_dicom, real) 
set_prj_real = t(land_prj$vb[-4, ]) 
thr = array(NA, dim = c(nland, 3, 2))
thr[, , 1] = set_prj_real
thr[, , 2] = set_dicom
dimnames(thr)[[3]] = c(paste("smoo_temp", "real"), paste("smoo_temp","dicom"))
patch=SL.set
nsland = dim(patch)[1]  
first_thr = createAtlas(real, as.matrix(set_prj_real), as.matrix(patch))
mesh2ply(real, paste("sur.tmp/smoo_temp", "real")) 
mesh2ply(dicom, paste("sur.tmp/smoo_temp", "dicom")) 
patched_thr <- placePatch(first_thr, thr, path = "sur.tmp/",fileext = ".ply",silent = TRUE) 
sink("NUL")
sliding_thr <- slider3d(patched_thr, SMvector = c(1:nland), 
                          surp = c(1:dim(patched_thr)[1])[-c(1:nland)], sur.path = "sur.tmp/", 
                          deselect = T, fixRepro = FALSE,mc.cores=num.cores) 
sink()
thr_eu = sum(sqrt(rowSums((patch - sliding_thr$dataslide[(nland + 1):(nland + nsland), , 2])^2))) 
header(verbose, "smoothing & sliding STEP", padding=0)
smooth.sets=array(NA, dim = c(nland + nsland, 3, 1))
registerDoParallel(cores = num.cores)
for(j in 1:length(algorithms)){
    algorithm = tolower(algorithms)[j]
    if (algorithm == "fujilaplace") {
      clean_smoothed_FJ=vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                                  type = algorithm, iteration = 1, delta = deltaFJ[1],lambda = deltaFJ[1]) 
      vertex_FJ = t(clean_smoothed_FJ$vb[-4, ]) 
      land_prj_FJ = t(projRead(set, clean_smoothed_FJ)$vb[-4, ]) 
      position_smoo = aro.clo.points(vertex_FJ, land_prj_FJ)} 
    
    if (algorithm == "angweight") {                
      clean_smoothed_AW=vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                                  type = algorithm, iteration = 1, delta = deltaAW[1],lambda = deltaAW[1])
      vertex_AW = t(clean_smoothed_AW$vb[-4, ])
      land_prj_AW = t(projRead(set, clean_smoothed_AW)$vb[-4, ])
      position_smoo = aro.clo.points(vertex_AW, land_prj_AW)}
    
    if (algorithm == "taubin") {                
      clean_smoothed_TA=vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                                  type = algorithm, iteration = 1, mu = mu[1],lambda = lambda[1])
      vertex_TA = t(clean_smoothed_TA$vb[-4, ])
      land_prj_TA = t(projRead(set, clean_smoothed_TA)$vb[-4, ])
      position_smoo = aro.clo.points(vertex_TA, land_prj_TA)}
    
    if (algorithm == "hclaplace") {                
      clean_smoothed_HC=vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                                  type = algorithm, iteration = 1)
      vertex_HC = t(clean_smoothed_HC$vb[-4, ])
      land_prj_HC = t(projRead(set, clean_smoothed_HC)$vb[-4, ])
      position_smoo = aro.clo.points(vertex_HC, land_prj_HC)}
    
    if (algorithm == "laplace") {                
      clean_smoothed_LA=vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                                  type = algorithm, iteration = 1)
      vertex_LA = t(clean_smoothed_LA$vb[-4, ])
      land_prj_LA = t(projRead(set, clean_smoothed_LA)$vb[-4, ])
      position_smoo = aro.clo.points(vertex_LA, land_prj_LA)}
    it_alg = seq(1, num_it[j], 1)
    algorithm = tolower(algorithms)[j]
    header(verbose, algorithm, padding=0)
    smooth.set=foreach(m = (it_alg),.combine=function(...) abind(..., along=3),.packages=c("Morpho","Rvcg","R.utils"))%:%
      foreach(i =1:iter,.combine=function(...) abind(..., along=3),.packages=c("Morpho","Rvcg","R.utils")) %dopar%{
        header(verbose, paste("algorithm:",algorithm,"setting",it_alg[m],"iter",i,"on",iter), padding=0)
        
        if (algorithm == "angweight") {
          smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "angweight", iteration = i, delta = deltaAW[m])
        }
        if (algorithm == "fujilaplace") {
          
          smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "fujilaplace", iteration = i, delta = deltaFJ[m], 
                               lambda = deltaFJ[m])
          
        }
        if (algorithm == "taubin") {
          smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = "taubin", iteration = i, lambda = lambda[m], 
                               mu = mu[m])
        }
        if (algorithm %in% c("laplace", "hclaplace")) {
          smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6),silent=TRUE), 
                               type = algorithm, iteration = i)
        }
        if(algorithm%in%c("hclaplace","laplace")){
          label_output_1 = paste(algorithm, 
                                 "00", i, sep = "")}
        if(algorithm%in%c("taubin","fujilaplace","angweight")){
          label_output_1 = paste(algorithm, it_alg[m], 
                                 "00", i, sep = "")}
        mesh2ply(smoothed, filename = paste("sur.tmp/",label_output_1, sep = "")) 
        vertex_smoo = t(smoothed$vb[-4, ]) 
        set_smoo =vertex_smoo[position_smoo, ] 
        array_smoo = array(NA, dim = c(nland, 3, 2))
        array_smoo[, , 1] = set_prj_real
        array_smoo[, , 2] = set_smoo
        dimnames(array_smoo)[[3]] = c(paste("smoo_temp","real"), paste(label_output_1))
        patched <- placePatch(first_thr, array_smoo,path = "sur.tmp/")
        sliding <- slider3d(patched, SMvector = c(1:nland), 
                            surp = c(1:dim(patched)[1])[-c(1:nland)], sur.path = "sur.tmp/", 
                            deselect = T, fixRepro = FALSE,mc.cores=num.cores)
        
        smooth.set = array(sliding$dataslide[, , 2], dim = c(nland + nsland, 3, 1))
      }
    smooth.sets=bindArr(smooth.sets,smooth.set,along=3)
  }
  
  set_alg=rep(substr(algorithms,1,3),num_it)
  if ("tau" %in% substr(tolower(algorithms), 1, 3) == TRUE) {
    set_alg[which(set_alg == "tau")] = paste("tau_", lambda, 
                                             sep = "")
  }
  if ("fuj" %in% substr(tolower(algorithms), 1, 3) == TRUE) {
    set_alg[which(set_alg == "fuj")] = paste("fuj_", deltaFJ, 
                                             sep = "")
  }
  if ("ang" %in% substr(tolower(algorithms), 1, 3) == TRUE) {
    set_alg[which(set_alg == "ang")] = paste("ang_", deltaAW, 
                                             sep = "")
  }
  sur_eu_dis=NULL
  for (i in 1:dim(smooth.sets)[[3]]){
    sur_eu_dis[i] = sum(sqrt(rowSums((patch - smooth.sets[(nland + 1):(nland + nsland), , i])^2)))
  }
  result_junk = matrix(sur_eu_dis[-c(1)], nrow = iter)
  result = (1 - (result_junk/thr_eu)) * 100 
  colnames(result) = set_alg
  rownames(result) = c(1:iter)
  best.set.pos = which(result == max(result), arr.ind = TRUE) 
  best.alg.type = substr(colnames(result)[best.set.pos[2]],1, 3) 
  best.alg.iter = as.numeric(substr(rownames(result)[best.set.pos[1]],1, 3)) 
  best.alg.value = result[which(result == max(result), arr.ind = TRUE)] 
  if (best.alg.value < 0) {
    print("neither of algorithm settings entail a recovery of anatomical information")
  }
  if (best.alg.value > 0) {
    print(paste("the type algorithm", best.alg.type, "entails a recovery equals to", 
          round(best.alg.value, 4), "% at the", "iteration number", 
          best.alg.iter))}
  if (best.alg.value > 0) {

  if (best.alg.type %in% c("tau") == TRUE) {
    best.typ.set = as.numeric(substr(colnames(result)[best.set.pos[2]],5,1000))
                  
    best.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6)), type = best.alg.type, 
                              iteration = best.alg.iter, lambda = best.typ.set, 
                              mu = -as.numeric(best.typ.set) - tau)
    mesh2ply(best.smoothed, filename = paste("type_", best.alg.type, 
                                             "_lambda_", best.typ.set, "_mu_", -as.numeric(best.typ.set) - 
                                               tau, "_iter_", best.alg.iter, sep = "")) 
  }
  
  if (best.alg.type %in% c("ang", "fuj") == TRUE) {
    best.typ.set = as.numeric(substr(colnames(result)[best.set.pos[2]],5,1000))
    best.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6)), type = best.alg.type, 
                              iteration = best.alg.iter, delta = best.typ.set)
    mesh2ply(best.smoothed, filename = paste("type_", best.alg.type, 
                                             "_delta_", best.typ.set, "_iter_", best.alg.iter, 
                                             sep = ""))
  }
  
  if (best.alg.type %in% c("lap", "hcl") == TRUE) {
    best.smoothed = vcgSmooth(vcgClean(dicom, sel = c(0:3,5,6)), type = best.alg.type, 
                              iteration = best.alg.iter)
    mesh2ply(best.smoothed, filename = paste("type_", best.alg.type, 
                                             "_iter_", best.alg.iter, sep = ""))
  }}

  
  unlink("sur.tmp", recursive = "TRUE") 
  shifting= sqrt(rowSums((set - set_prj_real)^2)) 
  shifting.mean=mean(shifting) 
  print(paste("mean shifting landmark equals to",shifting.mean,"mm"))
  shifting.range=range(shifting) 
  array.result = smooth.sets[,,-1]
  Lset.best.smoothed = array.result[c(1:nland), , which(max(as.vector(result)) == as.vector(result))] #landmark set of best smoothed deimated mesh
  SLset.best.smoothed = array.result[c((nland + 1):(nland + 
                                                      nsland)), , which(max(as.vector(result)) == as.vector(result))] #semi-landmark set on best smoothed decimeted mesh
  CS.vector = c(cSize(sliding_thr$dataslide[(nland +1):(nland + nsland), , 1]), 
                cSize(sliding_thr$dataslide[(nland +1):(nland + nsland), , 2]), 
                cSize(array.result[(nland +1):(nland + nsland), , which(max(as.vector(result)) == as.vector(result))]))
  names(CS.vector) = c("CS_iso_surface", "CS_decimated", "CS_best_smoothed")
  if (best.alg.value > 0) {
  out = list(LsetMeanShift=shifting.mean,LsetRangeShift=shifting.range,matrix_result = result, CS = CS.vector, L_set = Lset.best.smoothed, 
             SL_set = SLset.best.smoothed, mesh = best.smoothed)
  plot3D(array.result[, , which(max(as.vector(result)) == as.vector(result))], 
         cex = c(rep(1, nland), rep(0.5, nsland)), col = c(rep(2, 
                                                               nland), rep(3, nsland)))
    shade3d(out$mesh, col = 4)}
if (best.alg.value < 0) {
 out = list(LsetMeanShift=shifting.mean,LsetRangeShift=shifting.range,matrix_result = result, CS = CS.vector, L_set = Lset.best.smoothed, 
             SL_set = SLset.best.smoothed)}
  return(out)
}