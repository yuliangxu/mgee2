#' mgee2k
#'
#' Corrected GEE2 for ordinal data. This method yields unbiased estimators, but the
#' misclassification parameters are required to known.
#' 
#' \emph{mgee2k} implements the misclassification adjustment method outlined in Chen et al.(2014)
#'  where the misclassification parameters are known. In this case, validation data are not required,
#'  and only the observed data of the outcome and covariates are needed for the implementation.
#'  
#' @export
#' @useDynLib mgee2
#' @inheritParams ordGEE2
#' @param formula a formula object which specifies the relationship between 
#'   the response and covariates for the observed data.
#' @param data a dataframe or matrix object for the observed data set.
#' @param gamMat a matrix object which records the misclassification parameter gamma for response Y.
#' @param varphiMat a matrix object which records the misclassification parameter phi for covariate X.
#' @param misvariable a character object which names the error-prone covariate W.
#' @return  A list with component
#'     \item{beta}{the coefficients in the order as those specified in the formula for the response and covariates.}
#'     \item{alpha}{the oefficients for paired responses global odds ratios. The number of 
#'     alpha coefficients corresponds to the paired responses odds ratio structure selected
#'      in corstr. When corstr="exchangeable", only one baseline alpha 
#'      is fitted. When corstr="log-linear", baseline, first order, 
#'      second order (interaction) terms are fitted.}
#'     \item{variance}{variance-covariance matrix of the estimator of all parameters.}
#'     \item{convergence}{a logical variable; TRUE if the model converges.}
#'     \item{iteration}{the number of iterations for the estimates of the model parameters to converge.}
#'     \item{differ}{a list of difference of estimation for convergence}   
#'     \item{call}{Function called}
#' @examples
#'   if(0){
#'   data(obs1)
#'   obs1$visit <- as.factor(obs1$visit)
#'   obs1$treatment <- as.factor(obs1$treatment)
#'   obs1$S <- as.factor(obs1$S)
#'   obs1$W <- as.factor(obs1$W)
#'   ## set misclassification parameters to be known.
#'   varphiMat <- gamMat <- log( cbind(0.04/0.95, 0.01/0.95,
#'                                     0.95/0.03, 0.02/0.03,
#'                                     0.04/0.01, 0.95/0.01) )
#'   mgee2k.fit = mgee2k(formula = S~W+treatment+visit, id = "ID", data = obs1,
#'                     corstr = "exchangeable", misvariable = "W", gamMat = gamMat, 
#'                     varphiMat = varphiMat)
#'   }
#' @references
#' Z. Chen, G. Y. Yi, and C. Wu.  Marginal analysis of longitudinal ordinal data with misclassification inboth response and covariates. \emph{Biometrical Journal}, 56(1):69-85, Oct. 2014
#' 
#' Xu, Yuliang, Shuo Shuo Liu, and Y. Yi Grace. 2021. “mgee2: An R Package for Marginal Analysis of Longitudinal Ordinal Data with Misclassified Responses and Covariates.” \emph{The R Journal} 13 (2): 419.

mgee2k <- function(formula, id, data, corstr="exchangeable", misvariable,
                  gamMat, varphiMat, maxit=50, tol=1e-3)  {

  ID <- data[, as.character(id)]
  S.nam <- as.character(formula)[2]
  #fac.terms <- as.character(formula)[3]
  #misvariable <- strsplit(fac.terms, split="\\+")[[1]][1]

  # ID <- data[, as.character(id)]
  # S.nam <- strsplit(formula, "~")[[1]][1]
  # fac.terms <- strsplit(formula, "~")[[1]][2]
  # misvariable <- strsplit(fac.terms, split="\\+")[[1]][1]

  N <- length(ID)
  n <- length(unique(ID))
  clsize.vec <- as.vector(table(ID))
  m <- clsize.vec[1]
  K <- length(unique(data[,S.nam])) - 1
  K_x <- length(unique(data[,misvariable])) - 1
  f <- kronecker(ID, rep(1,K))

  DM <- getDM(formula=formula, data=data)
  DM <- cbind(matrix(rep(diag(K), sum(clsize.vec)),ncol=K, byrow=T),
              DM[kronecker(1:sum(clsize.vec),rep(1,K)),-1])
  p_b <- ncol(DM)

  Resp <- as.factor(getResp(formula=formula, data=data))
  Resp <- as.numeric(levels(Resp))[Resp]
  S.Mat <- matrix(0, nrow=N, ncol=K)
  for (k in 1:K) {
    S.Mat[, k] <- Resp==k
  }
  S <- as.vector(t(S.Mat))
  W.c <- as.factor(data[,misvariable])
  W.c <- as.numeric(levels(W.c))[W.c]
  W.Mat <- matrix(0, nrow=N, ncol=K_x)
  for (k in 1:K_x) {
    W.Mat[, k] <- W.c==k
  }
  W <- as.vector(t(W.Mat))

  ## Get unbiased surrogates
  #L <- matrix(rep(1,N),ncol=1)
  L_ij <- rbind(1)
  Ytilde <- as.vector( apply(t(S.Mat), 2,
                             FUN=get.UnbSurr_ij, L_ij=L_ij, gamMat=gamMat, K=K) )
  Ztilde <- get_Z(Ytilde, n, m, K)

  L.x_ij <- rbind(1)
  X.ast.Mat <- apply(t(W.Mat), 2, L_ij=L.x_ij,
                     FUN=get.UnbSurr_ij, gamMat=varphiMat, K=K_x)
  X.ast <- as.vector(X.ast.Mat)
  X.ast.ext.Mat <- rbind(1-apply(X.ast.Mat, 2, FUN="sum"), X.ast.Mat)
  X.ast.ext.Mat.spl <- lapply(split(data.frame(t(X.ast.ext.Mat)),
                                    f=ID), FUN="as.matrix")

  naigee2.fit <- ordGEE2(formula=formula, id=id,
                         corstr=corstr, data=data)
  beta_est.init <- naigee2.fit$beta
  alpha_est.init <- naigee2.fit$alpha
  p_a <- length(alpha_est.init)
  theta_est.init <- c(beta_est.init, alpha_est.init)
  p_t <- p_b + p_a

  X_i.ENUM <- NULL
  for (j in 1:m) {
    X_i.ENUM <- rbind( X_i.ENUM, rep(kronecker(0:K_x,
                                               rep(1,(1+K_x)^(m-j))), (1+K_x)^(j-1)) )
  }
  j1.indx <- NULL
  k1.indx <- NULL
  j2.indx <- NULL
  k2.indx <- NULL
  for(j1 in 1:(m-1)) {
    for(j2 in (j1+1):m) {
      j1.indx <- c(j1.indx, rep(j1, K^2))
      k1.indx <- c(k1.indx, kronecker(1:K, rep(1, K)))
      j2.indx <- c(j2.indx, rep(j2,K^2))
      k2.indx <- c(k2.indx, kronecker(rep(1,K), 1:K))
    }
  }
  if (corstr=="exchangeable") {
    #alpha_est.init <- 0.5
    #p_a <- length(alpha_est.init)
    assocDM <- matrix(rep(1, length(Ztilde)), ncol=p_a)
  } else {
    if (corstr=="log-linear") {
      #alpha_est.init <- rep(0.35, 1+(K-1)+K*(K-1)/2)
      #p_a <- length(alpha_est.init)
      assocDM_i <- cbind(rep(1, length(j1.indx)),
                         (as.numeric(k1.indx==2) +
                            as.numeric(k2.indx==2)),
                         as.numeric((k1.indx==2)&(k2.indx==2)))
      assocDM <- matrix(rep(t(assocDM_i), n), ncol=p_a, byrow=T)
    }
  }

  wgt <- sapply(X.ast.ext.Mat.spl,
                FUN=get.prod.X_i.ast.ENUM, X_i.ENUM=X_i.ENUM, m_i=m)

  nENUM <- (1+K_x)^m
  XDM_ENUM <- matrix(0, nrow=m*K, ncol=nENUM*K_x)
  X_l <- matrix(0, nrow=m, ncol=K_x)
  for (enum in 1:nENUM)
  {
    for (r in 1:K_x)
      X_l[, r] <- X_i.ENUM[, enum] == r
    XDM_ENUM[, 1:K_x+(enum-1)*K_x] <- matrix(
      kronecker(as.vector(X_l), rep(1, K)), ncol=K_x)
  }

  wgt <- t(wgt)

  beta <- beta_est.init
  alpha <- alpha_est.init
  theta <- theta_est.init

  dif <- 1
  differ <- c()  ## added to record the difference
  iter <- 0

  res <- list()
  res$beta <- rep(NA, p_b)
  res$alpha <- rep(NA, p_a)
  res$variance <- matrix(0, nrow=p_t, ncol=p_t)
  res$convergence <- FALSE
  res$iteration = 0

  while(iter<maxit & dif >tol) {
    beta_est.o <- beta
    alpha_est.o <- alpha
    theta_est.o <- theta

    nrow_DM_i <- m*K
    nrow_assocDM_i <- K^2 * (m*(m-1))/2

    U = rep(0, p_t)
    M = matrix(0, nrow=p_t, ncol=p_t)
    Sigma = matrix(0, nrow=p_t, ncol=p_t)
    for (i in 1:n) {
      ##  ---- Pass to C  ---- ##
      U_i <-  rep(0, p_t)
      M_i <- matrix(0, nrow=p_t, ncol=p_t)
      Sigma_i <- matrix(0, nrow=p_t, ncol=p_t)

      mgee2_i <- .C("Cgetmgee2_i",
                    as.double(DM[(i-1)*nrow_DM_i + 1:nrow_DM_i, ]),
                    as.double(Ytilde[(i-1)*nrow_DM_i + 1:nrow_DM_i]),

                    as.double(assocDM[(i-1)*nrow_assocDM_i + 1:nrow_assocDM_i, ]),
                    as.double(Ztilde[(i-1)*nrow_assocDM_i + 1:nrow_assocDM_i]),
                    as.double(wgt[i, ]),
                    as.double(XDM_ENUM),
                    as.integer(m),
                    as.integer(K),
                    as.integer(K_x),
                    as.integer(p_b),
                    as.integer(p_a),
                    beta = as.double(beta),
                    alpha = as.double(alpha),
                    U_i = as.double( U_i ),
                    M_i = as.double(M_i),
                    Sigma_i = as.double(Sigma_i)
      )
      U <- U + mgee2_i$U_i
      M <- M + matrix(mgee2_i$M_i, nrow=p_t, byrow=FALSE)
      Sigma <- Sigma + matrix(mgee2_i$Sigma_i,
                              nrow=p_t, byrow=FALSE)
    }
    U_beta <- 1/n * U[1:p_b]
    U_alpha <- 1/n * U[(p_b+1):p_t]
    M1 <- 1/n * M[1:p_b, 1:p_b]
    M2 <- as.matrix( 1/n * M[(p_b+1):p_t, (p_b+1):p_t] )
    if (any(is.na(M1)) | any(is.na(M2))) {
      return(res)
    } else {
      if ((abs(det(M1)) < 1e-15) | (abs(det(M2)) < 1e-10)) {
        #return (res)
        inv.M1 <- ginv(M1)
        inv.M2 <- ginv(M2)
      } else {
        if (any(eigen(M1)$values<=1e-5) | any(eigen(M2)$values<=1e-5)) {
          inv.M1 <- solve(M1)
          inv.M2 <- solve(M2)
        } else {
          inv.M1 <- chol2inv(chol(M1))
          inv.M2 <- chol2inv(chol(M2))
        }
      }
    }
    beta <- beta_est.o + 0.8 * inv.M1%*%U_beta
    alpha <- alpha_est.o + 0.5 * inv.M2%*%U_alpha
    theta <- c(beta, alpha)
    dif <- max(abs(theta/theta_est.o-1))
    differ <- c(differ, dif)  ##  keep the difference
    iter <- iter + 1
  }

  ## ---------------------
  ##   return the result
  ## ---------------------

  beta_nam <- c(unlist(lapply("Y>=", 1:K, FUN=paste, sep="")),
                dimnames(DM)[[2]][-(1:K)])
  alpha_nam <- NULL
  if (corstr=="exchangeable") {
    alpha_nam <- "Delta"
  } else {
    if (corstr=="log-linear") {
      alpha_nam <- c("Delta", "Delta_2", "Delta_22")
    }
  }

  Gamma <- 1/n * M
  Sigma <- 1/n * Sigma

  res$beta <- as.vector(beta)
  res$alpha <- as.vector(alpha)
  if (abs(det(Gamma)) < 1e-15) {
    res$variance <- ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))/n
  } else {
    inv.Gamma <- matrix(0, nrow=p_t, ncol=p_t)
    inv.Gamma[1:p_b, 1:p_b] <- inv.M1
    inv.Gamma[(p_b+1):p_t, (p_b+1):p_t] <- inv.M2
    inv.Gamma[(p_b+1):p_t, 1:p_b] <- -( inv.M2%*%
                                          Gamma[(p_b+1):p_t, 1:p_b]%*%inv.M1 )
    res$variance <- inv.Gamma%*% Sigma %*% t(inv.Gamma)/n
  }
  res$convergence <- (iter<maxit & dif < tol)
  res$iteration = iter
  res$differ <- differ ## return the difference
  names(res$beta) <- beta_nam
  names(res$alpha) <- alpha_nam
  dimnames(res$variance) <- list(c(beta_nam, alpha_nam),
                                 c(beta_nam, alpha_nam))
  # Yuliang: include call in the output
  res$call <- match.call()
  class(res) <- "mgee2"
  return(res)
}### end of mgee2()
