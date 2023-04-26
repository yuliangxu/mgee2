#' ordGEE2
#'
#' This function provides a naive approach to estimate the data without any correction
#' or misclassification parameters. This may lead to biased estimation for response
#' parameters.
#' 
#' In addition to developing the package \emph{mgee2} to implement the methods of Chen et al.(2014) which accommodate 
#' misclassification effects in inferential procedures, we also implement the naive method of ignoring the feature of misclassification, 
#' and call the resulting function \emph{ordGEE2}. 
#' This function can be used together with the precedingly described \emph{mgee2k} or \emph{mgee2v} to evaluate the impact of not 
#' addressing misclassification effects
#' @export
#' @useDynLib mgee2
#' @param formula a formula object: a symbolic description of the model with error-prone
#'   response, error-prone covariates and other covariates.
#' @param id a character object which records individual id in the data.
#' @param data a dataframe or matrix of the observed data, including id, error-prone ordinal response
#'   error-prone ordinal covaritaes, other covariates.
#' @param corstr a character object. The default value is "exchangeable", 
#'   corresponding to the structure where the association between two paired 
#'   responses is considered to be a constant. The other option is "log-linear" 
#'   which  indicates the log-linear association between two paired responses.
#' @param maxit an integer which specifies the maximum number of iterations. The default is 50.
#' @param tol a numeric object which indicates the tolerance threshold. The default is 1e-3.
#' 
#'
#' @examples
#'   data(obs1)
#'   obs1$Y <- as.factor(obs1$Y)
#'   obs1$X <- as.factor(obs1$X)
#'   obs1$visit <- as.factor(obs1$visit)
#'   obs1$treatment <- as.factor(obs1$treatment)
#'   obs1$S <- as.factor(obs1$S)
#'   obs1$W <- as.factor(obs1$W)
#'   naigee.fit = ordGEE2(formula = S~W+treatment+visit, id = "ID",
#'                        data = obs1, corstr = "exchangeable")
#'
#' @return  A list with component
#'     \item{beta}{the coefficients in the order of 1) all non-baseline levels for response,
#'     2) covariates - same order as specified in the formula}
#'     \item{alpha}{the coefficients for paired responses global odds ratios. Number of alpha
#'     coefficients corresponds to the paired responses odds ratio structure selected in "corstr";
#'     when corstr="exchangeable", only one baseline alpha is fitted.}
#'     \item{variance}{variance-covariance matrix of all fitted parameters}
#'     \item{convergence}{a logical variable, TRUE if the model converges}
#'     \item{iteration}{number of iterations for the model to converge}
#'     \item{differ}{a list of difference of estimation for convergence}   ##
#'     \item{call}{Function called}
#' @references
#' Z. Chen, G. Y. Yi, and C. Wu. Marginal analysis of longitudinal ordinal data with misclassification inboth response and covariates. \emph{Biometrical Journal}, 56(1):69-85, Oct. 2014
#' 
#' Xu, Yuliang, Shuo Shuo Liu, and Y. Yi Grace. 2021. “mgee2: An R Package for Marginal Analysis of Longitudinal Ordinal Data with Misclassified Responses and Covariates.” \emph{The R Journal} 13 (2): 419.


ordGEE2 <- function(formula, id, data, corstr="exchangeable",
                    maxit=50, tol=1e-3)  {


  ID <- data[, as.character(id)]
  N <- length(ID)
  n <- length(unique(ID))
  clsize.vec <- as.numeric(table(ID))
  m <- clsize.vec[1]

  Resp <- as.factor(getResp(formula=formula, data=data))
  Resp <- as.numeric(levels(Resp))[Resp]
  K <- length(unique(Resp))-1
  Y.Mat <- matrix(0, nrow=N, ncol=K)
  for (k in 1:K) {
    Y.Mat[, k] <- Resp==k
  }

  Y <- as.vector(t(Y.Mat))
  Z <- get_Z(Y, n=n, m=m, K=K)
  f <- kronecker(ID, rep(1,K))

  # Yuliang: getDM - set contrasts=NULL for model.matrix in getDM
  DM <- getDM(formula=formula, data=data)
  DM <- cbind(matrix(rep(diag(K), sum(clsize.vec)),ncol=K, byrow=T),
              DM[kronecker(1:sum(clsize.vec),rep(1,K)),-1])
  p_b <- ncol(DM)

  # Yuliang: change formula type to class:formula
  #fac.terms <- strsplit(formula, "~")[[1]][2] # covariates formula
  #form.k <- paste("as.numeric(Resp>=K)", fac.terms, sep="~")
  form.k <- paste("as.numeric(Resp>=K)", as.character(formula)[3], sep="~")
  init.glm <- glm(formula=eval(parse(text=form.k)),
                  family=binomial,data=data)
  beta_est.init <- coefficients(init.glm)
  if(K>1){
    for(k in (K-1):1){
      form.k <- paste("as.numeric(Resp>=k)", as.character(formula)[3], sep="~")
      init.glm <- glm(formula=eval(parse(text=form.k)),
                      family=binomial,data=data)
      beta_est.init <- c(coefficients(init.glm)[1], beta_est.init)
    }
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
  indx.tab <- cbind(j1.indx, k1.indx, j2.indx, k2.indx)
  indx.tab.t <- t(indx.tab)
  if (corstr=="exchangeable") {
    alpha_est.init <- 0.5
    p_a <- 1
    assocDM <- matrix(rep(1, length(Z)), ncol=p_a)
  } else {
    if (corstr=="log-linear") {
      p_a <- K + K * (K - 1)/2
      alpha_est.init <- c(log(3), rep(0, p_a-1))

      assocDM_i <- cbind(rep(1, length(j1.indx)),
                         (as.numeric(k1.indx==2) +
                            as.numeric(k2.indx==2)),
                         as.numeric((k1.indx==2)&(k2.indx==2)))
      assocDM <- matrix(rep(t(assocDM_i), n), ncol=p_a, byrow=T)
    }
  }

  p_t <- p_b + p_a
  theta_est.init <- c(beta_est.init, alpha_est.init)

  beta <- beta_est.init
  alpha <- alpha_est.init
  theta <- theta_est.init

  dif <- 1
  differ <- c()  ## added to record the difference
  iter <- 0

  res <- list()
  res$beta <- numeric(p_b)
  res$alpha <- numeric(p_a)
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

      ordgee2_i <- .C("Cgetordgee2_i",
                      as.double(DM[(i-1)*nrow_DM_i + 1:nrow_DM_i, ]),
                      as.double(Y[(i-1)*nrow_DM_i + 1:nrow_DM_i]),
                      as.double(assocDM[(i-1)*nrow_assocDM_i + 1:nrow_assocDM_i, ]),
                      as.double(Z[(i-1)*nrow_assocDM_i + 1:nrow_assocDM_i]),
                      as.integer(m),
                      as.integer(K),
                      as.integer(p_b),
                      as.integer(p_a),
                      beta = as.double(beta),
                      alpha = as.double(alpha),
                      U_i = as.double( U_i ),
                      M_i = as.double(M_i),
                      Sigma_i = as.double(Sigma_i)
      )
      U <- U + ordgee2_i$U_i
      M <- M + matrix(ordgee2_i$M_i, nrow=p_t, byrow=FALSE)
      Sigma <- Sigma + matrix(ordgee2_i$Sigma_i,
                              nrow=p_t, byrow=FALSE)
    }
    U_beta <- 1/n * U[1:p_b]
    U_alpha <- 1/n * U[(p_b+1):p_t]
    M1 <- 1/n * M[1:p_b, 1:p_b]
    M2 <- as.matrix( 1/n * M[(p_b+1):p_t, (p_b+1):p_t] )
    if (any(is.na(M1)) | any(is.na(M2))) {
      return(res)
    } else {
      if ((abs(det(M1)) < 1e-15) | (abs(det(M2)) < 1e-15)) {
        return (res)
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
    alpha <- alpha_est.o + 0.8 * inv.M2%*%U_alpha
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
}
