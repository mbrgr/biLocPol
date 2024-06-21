#### Leave One Plane Out Cross Validation #####

#' lopocv
#' @description Leave One Plane Out Cross Validation.
#' @param Y Raw observations without any transformation. Numerical with dimension n x p.
#' @param h.seq Sequence of bandwidth to compare in cross validation. Numerical vector.
#' @param m Degree of local polynomial estimation. 1 and 2 is possible.
#' @param h.parallel Logical value
#' @param h.parallel.environment Logical value that indicates if the 'makeCluster', 'plan' and 'stopCluster' functions shall be called
#' @param ... Further arguments passend to 'local_polynomial_weights' such as 'parallel' and 'parallel.environment'
#'
#' @return Bandwidth with minimal sup error in cross validation
#' @importFrom future.apply future_apply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom future plan
#' @export
#'
#' @examples 0 # TODO:
lopocv = function(Y, h.seq, m = 2, h.parallel = F, h.parallel.environment = F, na.rm = F, ...){
  p = length(Y[1,])
  n = length(Y[,1])
  Y = apply(Y, 2, function(x){x - mean(x)})

  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, ...)
    max_diff = numeric(n)
    for(l in 1:n){
      Z_l         = tcrossprod(Y[l,], Y[l,])
      Z_minus_l   = observation_transformation(Y[-l,], na.rm = na.rm)
      estimate    = eval_weights(w_h, Z_minus_l)
      max_diff[l] = max( abs(Z_l - estimate)[!as.logical(diag(p))][!is.na(Z_l)] )
    }
    mean(max_diff)
  }

  if (h.parallel) {
    if(h.parallel.environment){
      cl = parallel::makeCluster(parallel::detectCores( ) - 1)
      future::plan(future::cluster)
    }
    mean_sup = future_sapply(h.seq, help, future.seed = T)
    if(h.parallel.environment){
      parallel::stopCluster(cl)
    }
  } else {
    mean_sup = sapply(h.seq, help)
  }
  h.seq[which.min(mean_sup)]

}


#### K-Fold Cross Validation #####

#' k_fold_c
#' @description K-Fold Cross Validation
#'
#' @param Y Raw observations without any transformation. Numerical with dimension n x p.
#' @param h.seq Sequence of bandwidth to compare in cross validation. Numerical vector.
#' @param K Into how many subset shall the data be devided?
#' @param m Degree of local polynomial estimation. 1 and 2 is possible.
#' @param h.parallel Logical value wether the
#' @param h.parallel.environment Logical value that indicates if the 'makeCluster', 'plan' and 'stopCluster' functions shall be called
#' @param ... Further arguments passend to 'local_polynomial_weights'.
#'
#' @return Bandwidth with minimal sup error in cross validation
#' @importFrom future.apply future_apply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom future plan
#' @export
#'
#' @examples 0 # TODO:
k_fold_cv = function(Y, h.seq, K = 2, m = 1, h.parallel = F, h.parallel.environment = F, na.rm = F,...){

  p = length(Y[1,])
  n = length(Y[,1])

  grp = sample(rep(1:K, ceiling(n/K)), n)

  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, ...)
    max_diff = numeric(K)
    for(kk in 1:K){
      test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full", na.rm = T), p, p)
      train_grp    = observation_transformation(Y[grp != kk,], na.rm = na.rm)
      pred         = eval_weights(w_h, train_grp)
      max_diff[kk] = max(abs(test_grp - pred)[!as.logical(diag(p))])
    }
    mean(max_diff)
  }

  if(h.parallel){
    if(h.parallel.environment){
      cl = parallel::makeCluster(parallel::detectCores( ) - 1)
      future::plan(future::cluster)
    }
    mean_sup = future_sapply(h.seq, help, future.seed = T)
    if(h.parallel.environment){
      parallel::stopCluster(cl)
    }
  }else{
    mean_sup = sapply(h.seq, help)
  }
  h.seq[which.min(mean_sup)]
}





