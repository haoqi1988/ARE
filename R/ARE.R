#' Factorize of a Tensor
#'
#' This function factorizes a tensor into the sum of two products. The first term is the product of attributes
#' and attribute relations, and the second term is the product of the indirect relations and
#' meta-relations between indirect relations and direct relations.
#'
#' @param N,K,r,p,entity_names,link_names,attribute_names,indirect_link_names,lamdaA,lamdaR,lamdaW,X,M,error_tolerance 14 inputs
#' @return A list containing A, the attributes, R, attribute relations and W, the meta-relations
#' @export
# start defing the factorization function here ======================
ARE <- function(N, # number of entity
                K, # number of layers/link types
                r, # number of attribute types
                p, # number of indirect link types
                entity_names,
                link_names,
                attribute_names,
                indirect_link_names,
                lamdaA, # weight of A
                lamdaR, # weight of R
                lamdaW, # weight of W
                X,      # here X is a 3-d array, does not need to be a tensor
                M,      # here M is a 3-d array, does not need to be a tensor
                error_tolerance){

  # create random initial R, A, W according to the dims of X and M
  # create matrix W, dim(W) = K*P, fill in random values (here i am using 0's)
  W <- array(rnorm(K*p), c(K,p),dimnames = list(link_names,indirect_link_names))
  # create tensor R, dim(R) = r*r*K, the weights of attribue pairs in predicting direct links,
  R <- array(rnorm(r*r*K),c(r,r,K),dimnames = list(attribute_names,attribute_names,link_names))    # random value is fine
  # create matrix A, dim(A) = N*r, the weights of attribue pairs in predicting direct links
  A <- array(rnorm(N*r), c(N,r), dimnames = list(entity_names,attribute_names))    # random value is fine

  # start the machine learning - estimate R A W based on X M
  # prepare data - this part requires package 'rTensor'
  R <- rTensor::as.tensor(R, drop = F)
  X <- rTensor::as.tensor(X, drop = F)
  M <- rTensor::as.tensor(M, drop = F)
  X3 <- rTensor::unfold(X, 3, c(2,1))
  M3 <- rTensor::unfold(M, 3, c(2,1))
  M3d <- M3@data
  X3d <- X3@data

  min_error_change <- error_tolerance + 1
  min_error_list<- 0
  while (min_error_change >= error_tolerance) {

    ### estimae R, A, W ===================================================
    # estimate W
    R3 <- rTensor::unfold(R, 3, c(2,1))
    R3d <- R3@data
    product <- rTensor::ttm(M, t(A), m=1)
    product <- rTensor::ttm(product, t(A), m=2)
    product3 <- rTensor::unfold(product,3, c(2,1))
    product3 <- product3@data
    W <- (X3d%*%t(M3d) - R3d%*%t(product3)) %*% solve(M3d %*% t(M3d)+ lamdaW*diag(p))

    ### estimate A
    E <- X - rTensor::ttm(M, W, m = 3)
    sumEAR <- 0
    for (k in 1:K) {
      Ek <- E[,,k]@data
      Rk <- R[,,k]@data
      sumEAR <- sumEAR + Ek %*% A %*% t(Rk) + t(Ek) %*% A %*% Rk
    }
    sumRAAR <- 0
    for (k in 1:K){
      Rk <- R[,,k]@data
      sumRAAR <- sumRAAR + Rk %*% t(A) %*% A %*% t(Rk) + t(Rk) %*% t(A) %*% A %*% Rk
    }
    A <- sumEAR %*% solve(sumRAAR+ diag(r))

    ### estimate R
    SVD <- svd(A)
    d <- SVD$d
    U <- SVD$u
    V <- SVD$v
    S <- matrix(d%x%d, nrow = r,ncol = r, byrow = T)
    S <- S/(S^2 + lamdaR)
    for (k in 1:K){
      Ek <- E[,,k]@data
      R[,,k] <- V %*% (S * (t(U) %*% Ek %*% U)) %*% t(V)
    }

    ### calculate errors
    R1A <- rTensor::ttm(R, A, m=1)
    R1A2A <- rTensor::ttm(R1A, A, m=2)
    min_error <- rTensor::fnorm(E - R1A2A)^2 + lamdaA*norm(A, "f")^2 + lamdaR*rTensor::fnorm(R)^2 + lamdaW*norm(W, 'f')^2
    ### iterate to here. Finish==========================================================

    min_error_list <- c(min_error_list,min_error)
    min_error_change <- abs(min_error_list[length(min_error_list)] - min_error_list[length(min_error_list)-1])
  }
  return(list('A' = A,
              'R' = R,
              'W' = W,
              'min_error_list' = min_error_list))
}
# end function definition here ====


# start defining the PlotPRCurve() function ======
PlotPRCurve <- function(A,R,W,X,M){
  M <- rTensor::as.tensor(M, drop = F)
  X_predicted <- rTensor::ttm(rTensor::ttm(R, A, m = 1), A, m = 2) + rTensor::ttm(M, W, m = 3)
  NumOnes <- 0 # compute the number of 1's in X
  for (k in 1:dim(X)[3]){
    NumOnes <- NumOnes + sum(X[,,k])
  }
  fg <- rep(0,NumOnes)   ### vector of fore ground predictions
  bg <- rep(0,N*N*K - NumOnes)  ### first need to find out the number of 0's in the observed data
  ifg <- 1
  ibg <- 1
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[2]){
      for (k in 1:dim(X)[3]){
        if (X[i,j,k] == 1){
          fg[ifg] <- X_predicted[i,j,k]@data
          ifg <- ifg + 1
        } else {
          bg[ibg] <- X_predicted[i,j,k]@data
          ibg <- ibg + 1
        }
      }
    }
  }
  return(PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE))
}
# end defining function =====
