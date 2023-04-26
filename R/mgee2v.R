#' mgee2v
#'
#' Corrected GEE2 for ordinal data, with validation subsample
#'
#'The function \emph{mgee2v} does not require the misclassification parameters to be known, 
#'but require the availability  of validation data.
#' Similar to \emph{mgee2k}, the function \emph{mgee2v} needs the data set to be structured by individual id, i=1,...,n, and visit time, j_i=1,...,m_i. 
#' The data set should contain the observed response and covariates S and W. 
#' To indicate whether or not a subject is in the validation set, an indicator variable 
#' \emph{delta} should be added in the data set, and we use a column named \emph{valid.sample.ind} for this purpose. 
#' The column name of the error-prone covariate W should also be specified in \emph{misvariable}. 
#' @export
#' @useDynLib mgee2
#' @inheritParams mgee2k
#' @param valid.sample.ind  a string object which names the indicator variable delta. 
#' When a data point belongs to the validation set, delta = 1; otherwise 0.
#' @param y.mcformula a string object which indicates the misclassification formula
#'  between true response Y and surrogate(observed) response S.
#' @param x.mcformula a string object which indicates the  misclassification formula
#'   between true error-prone covariate X and surrogate W.
#' @return  A list with component
#'     \item{beta}{the coefficients in the order of 1) all non-baseline levels for response,
#'     2) covariates - same order as specified in the formula}
#'     \item{alpha}{the coefficients for paired responses global odds ratios. Number of alpha
#'     coefficients corresponds to the paired responses odds ratio structure selected in "corstr";
#'     when corstr="exchangeable", only one baseline alpha is fitted.}
#'     \item{variance}{variance-covariance matrix of all fitted parameters}
#'     \item{convergence}{a logical variable, TRUE if the model converges}
#'     \item{iteration}{number of iterations for the model to converge}
#'     \item{call}{Function called}
#' @examples
#'   if(0){
#'   data(obs1)
#'   obs1$Y <- as.factor(obs1$Y)
#'   obs1$X <- as.factor(obs1$X)
#'   obs1$visit <- as.factor(obs1$visit)
#'   obs1$treatment <- as.factor(obs1$treatment)
#'   obs1$S <- as.factor(obs1$S)
#'   obs1$W <- as.factor(obs1$W)
#'   mgee2v.fit = mgee2v(formula = S~W+treatment+visit, id = "ID", data = obs1,
#'                       y.mcformula = "S~1", x.mcformula = "W~1", misvariable = "W",
#'                       valid.sample.ind = "delta",
#'                       corstr = "exchangeable")
#'   }
#' @references
#' Z. Chen, G. Y. Yi, and C. Wu. Marginal analysis of longitudinal ordinal data with misclassification inboth response and covariates. \emph{Biometrical Journal}, 56(1):69-85, Oct. 2014
#' 
#' Xu, Yuliang, Shuo Shuo Liu, and Y. Yi Grace. 2021. “mgee2: An R Package for Marginal Analysis of Longitudinal Ordinal Data with Misclassified Responses and Covariates.” \emph{The R Journal} 13 (2): 419.

mgee2v <- function(formula, id, data, corstr="exchangeable", misvariable = "W", valid.sample.ind = "delta",
                   y.mcformula, x.mcformula, maxit=50, tol=1e-3)  {

  #fixed.gamma <- FALSE

  ## Initial estimates
  naigee2.fit <- ordGEE2(formula=formula, id=id,
                         corstr=corstr, data=data)
  beta_est.init <- naigee2.fit$beta
  alpha_est.init <- naigee2.fit$alpha
  p_b <- length(beta_est.init)
  p_a <- length(alpha_est.init)

  ID <- data[, as.character(id)]
  S.nam <- as.character(formula)[2]
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

  S.c <- as.factor(data[, S.nam])
  Y.levs <- levels(S.c)
  W.c <- as.factor(data[,misvariable])
  X.levs <- levels(W.c)
  delta <- data[,valid.sample.ind]

  ## Estimate misclassification parameters (need modification)
  vss <- data[delta==1, ]
  gamMat.est <- matrix(0, nrow=1, ncol=K*(K+1))
  gamMat.var <- matrix(0, nrow=1, ncol=K*(K+1))
  for (k in 0:K) {
    for (l in 1:K) {
      kl.indx <- ( (vss$Y == Y.levs[k+1]) &
                     (vss[, S.nam] == Y.levs[l+1]) )
      k0.indx <- ( (vss$Y == Y.levs[k+1]) &
                     (vss[, S.nam] == Y.levs[1]) )
      vss.k <- vss[kl.indx | k0.indx, ]
      if ( (sum(kl.indx) * sum(k0.indx)) != 0 ) {
        glm.k <- glm(y.mcformula, family=binomial, data=vss.k)
        gamMat.est[, k*K+l] <- coef(glm.k)
        gamMat.var[, k*K+l] <- summary(glm.k)$cov.un
      } else {
        if ( abs(k - l) < abs(k - 0) ) {
          gamMat.est[, k*K+l] <- rep(
            log(1e10)/nrow(gamMat.est), nrow(gamMat.est))
        } else {
          gamMat.est[, k*K+l] <- rep(
            log(1e-10)/nrow(gamMat.est), nrow(gamMat.est))
        }
      }
    }
  }

  varphiMat.est <- matrix(0, nrow=1, ncol=K_x*(K_x+1))
  varphiMat.var <- matrix(0, nrow=1, ncol=K_x*(K_x+1))
  for (k in 0:K_x) {
    for (l in 1:K_x) {
      kl.indx <- ( (vss$X == X.levs[k+1]) &
                     (vss[, misvariable] == X.levs[l+1]) )
      k0.indx <- ( (vss$X == X.levs[k+1]) &
                     (vss[, misvariable] == X.levs[1]) )
      vss.k <- vss[kl.indx | k0.indx, ]
      if ( (sum(kl.indx) * sum(k0.indx)) != 0 ) {
        glm.k <- glm(x.mcformula, family=binomial, data=vss.k)
        varphiMat.est[, k*K_x+l] <- coef(glm.k)
        varphiMat.var[, k*K_x+l] <- summary(glm.k)$cov.un
      } else {
        if ( abs(k - l) < abs(k - 0) ) {
          varphiMat.est[, k*K_x+l] <- rep(
            log(.05/.01)/nrow(varphiMat.est), nrow(varphiMat.est))
        } else {
          varphiMat.est[, k*K_x+l] <- rep(
            log(.01/.05)/nrow(varphiMat.est), nrow(varphiMat.est))
        }
      }
    }
  }

  S.Mat <- matrix(0, nrow=N, ncol=K)
  Y.Mat <- matrix(0, nrow=N, ncol=K)
  for (k in 1:K) {
    S.Mat[, k] <- as.numeric(S.c == Y.levs[k+1])
    Y.Mat[, k] <- as.numeric(data$Y == Y.levs[k+1])
  }
  S <- as.vector(t(S.Mat))
  W.Mat <- matrix(0, nrow=N, ncol=K_x)
  X.Mat <- matrix(0, nrow=N, ncol=K_x)
  for (k in 1:K_x) {
    W.Mat[, k] <- as.numeric(W.c == X.levs[k+1])
    X.Mat[, k] <- as.numeric(data$X == X.levs[k+1])
  }

  ## -----------------------------------------
  ## Get J_i.all and Q_i.all
  ## -----------------------------------------
  gamMat <- gamMat.est
  varphiMat <- varphiMat.est
  nrow_DM_i <- m*K
  p_g <- prod(dim(gamMat))
  p_vp <- prod(dim(varphiMat))

  Q_g_i.all <- matrix(0, nrow=p_g, ncol=n)
  J_g <- matrix(0, nrow=p_g, ncol=p_g)
  Q_vp_i.all <- matrix(0, nrow=p_vp, ncol=n)
  J_vp <- matrix(0, nrow=p_vp, ncol=p_vp)

  ## -----------------------------------------
  ## Get corrected versions of Y, X
  ## -----------------------------------------
  L_ij <- rbind(1)
  Ytilde.Mat <- t( apply(t(S.Mat), 2,
                         FUN=get.UnbSurr_ij, L_ij=L_ij, gamMat=gamMat, K=K) )
  Ytilde.Mat[delta==1, ] <- Y.Mat[delta==1, ]
  Ytilde <- as.vector(t(Ytilde.Mat))
  Ztilde <- get_Z(Ytilde, n, m, K)
  P.spl <- list()
  invPast.spl <- list()
  L_x_ij <- rbind(1)
  Ytilde.Mat.ext.spl <- list()
  G.spl <- list()
  invGast.spl <- list()
  Xtilde.ext.Mat.spl <- list()
  DM.spl <- list()
  for (i in 1:n)
  {
    P_i <- matrix(0, K+1, m*(K+1))
    invPast_i <- matrix(0, K, m*K)
    Ytilde_i.Mat <- matrix(0, m, K)
    G_i <- matrix(0, K_x+1, m*(K_x+1))
    invGast_i <- matrix(0, K_x, m*K_x)
    Xtilde_i.Mat <- matrix(0, m, K_x)
    for (j in 1:m)
    {
      ## correct Y
      if (delta[(i-1)*m+j]!=1)
      {
        P_ij <- get.mcPrMat_ij(S.Mat[(i-1)*m+j, ],
                               L_ij=L_ij, gamMat=gamMat,
                               K=K, adjacent=FALSE)
        P_i[, (j-1)*(K+1)+(1:(K + 1))] <- P_ij
        invPast_ij <- solve(t(P_ij[2:(K+1), 2:(K+1)])
                            - P_ij[1, 2:(K+1)])
        invPast_i[, (j-1)*K+(1:K)] <- invPast_ij
        Ytilde_i.Mat[j, ] <- (invPast_ij%*%
                                (S.Mat[(i-1)*m+j, ] - P_ij[1, 2:(K+1)]))
      } else {
        Ytilde_i.Mat[j, ]<- Y.Mat[(i-1)*m+j, ]
        ## Get Q_g_i and J_g
        Y_ij_ext <- c(1 - sum(Y.Mat[(i-1)*m+j, ]),
                      Y.Mat[(i-1)*m+j, ])
        S_ij <- S.Mat[(i-1)*m+j, ]
        k <- which(Y_ij_ext==1) - 1
        tau_ij <- as.vector(exp(rbind(L_ij)%*%
                                  gamMat[, k * K + 1:K]))
        tau_ij <- tau_ij/(1+sum(tau_ij))
        dtau_ij_dgamT <- matrix(0, nrow=K, ncol=p_g)
        dtau_ij_dgamT[, k * K + 1:K] <- kronecker(
          diag(tau_ij*(1-tau_ij)), rbind(L_ij))
        D_g_ij <- t(dtau_ij_dgamT)
        V_g_ij <- diag(tau_ij) - outer(tau_ij, tau_ij)
        Q_g_i.all[, i] <- Q_g_i.all[, i] +
          D_g_ij%*%solve(V_g_ij)%*%(S_ij - tau_ij)
        J_g <- J_g - D_g_ij%*%solve(V_g_ij)%*%dtau_ij_dgamT
      }

      ## Correct X
      if (delta[(i-1)*m+j]!=1)
      {
        G_ij <- get.mcPrMat_ij(W.Mat[(i-1)*m+j, ],
                               L_ij=L_x_ij, gamMat=varphiMat,
                               K=K_x, adjacent=FALSE)
        G_i[, (j-1)*(K_x+1)+(1:(K_x + 1))] <- G_ij
        invGast_ij <- solve(t(G_ij[2:(K_x+1), 2:(K_x+1)])
                            - G_ij[1, 2:(K_x+1)])
        invGast_i[, (j-1)*K_x+(1:K_x)] <- invGast_ij
        Xtilde_i.Mat[j, ] <- (invGast_ij%*%
                                (W.Mat[(i-1)*m+j, ] - G_ij[1, 2:(K_x+1)]))
      } else {
        Xtilde_i.Mat[j, ]<- X.Mat[(i-1)*m+j, ]

        ## Get Q_vp_i and J_vp
        X_ij_ext <- c(1 - sum(X.Mat[(i-1)*m+j, ]),
                      X.Mat[(i-1)*m+j, ])
        W_ij <- W.Mat[(i-1)*m+j, ]
        r <- which(X_ij_ext==1) - 1
        pi_ij <- as.vector(exp(rbind(L_x_ij)%*%
                                 varphiMat[, r * K_x + 1:K_x]))
        pi_ij <- pi_ij/(1+sum(pi_ij))
        dpi_ij_dvpT <- matrix(0, nrow=K_x, ncol=p_vp)
        dpi_ij_dvpT[, r * K_x + 1:K_x] <- kronecker(
          diag(pi_ij*(1-pi_ij)), rbind(L_x_ij))
        D_vp_ij <- t(dpi_ij_dvpT)
        V_vp_ij <- diag(pi_ij) - outer(pi_ij, pi_ij)
        Q_vp_i.all[, i] <- Q_vp_i.all[, i] +
          D_vp_ij%*%solve(V_vp_ij)%*%(W_ij - pi_ij)
        J_vp <- J_vp - D_vp_ij%*%solve(V_vp_ij)%*%dpi_ij_dvpT
      }
    }
    P.spl[[i]] <- P_i
    invPast.spl[[i]] <- invPast_i
    Ytilde.Mat.ext.spl[[i]] <- cbind(1- apply(Ytilde_i.Mat, 1, sum),
                                     Ytilde_i.Mat)
    G.spl[[i]] <- G_i
    invGast.spl[[i]] <- invGast_i
    Xtilde.ext.Mat.spl[[i]] <- cbind(1 - apply(Xtilde_i.Mat, 1, sum),
                                     Xtilde_i.Mat)
    DM.spl[[i]] <- DM[(i-1)*nrow_DM_i + 1:nrow_DM_i, ]
  }
  J_g <- J_g/n
  J_vp <- J_vp/n

  X_i.ENUM <- NULL
  for (j in 1:m) {
    X_i.ENUM <- rbind( X_i.ENUM, rep(kronecker(0:K_x,
                                               rep(1,(1+K_x)^(m-j))), (1+K_x)^(j-1)) )
  }

  ## Get the weight for each possible X_i
  wgt <- sapply(Xtilde.ext.Mat.spl,
                FUN=get.prod.X_i.ast.ENUM, X_i.ENUM=X_i.ENUM, m_i=m)
  wgt <- t(wgt)

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

  ## Association
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
    assocDM <- matrix(rep(1, length(Ztilde)), ncol=p_a)
  } else {
    if (corstr=="log-linear") {
      assocDM_i <- cbind(rep(1, length(j1.indx)),
                         (as.numeric(k1.indx==2) +
                            as.numeric(k2.indx==2)),
                         as.numeric((k1.indx==2)&(k2.indx==2)))
      assocDM <- matrix(rep(t(assocDM_i), n), ncol=p_a, byrow=T)
    }
  }
  nrow_assocDM_i <- K^2 * (m*(m-1))/2

  theta_est.init <- c(beta_est.init, alpha_est.init)
  p_t <- p_b + p_a

  ## -------------------------------------
  ## The core of the estimation procedure
  ## -------------------------------------
  beta <- beta_est.init
  alpha <- alpha_est.init
  theta <- theta_est.init
  if (abs(det(J_g))<1e-10) {
    invJ_g <- ginv(J_g)
  } else {
    invJ_g <- solve(J_g)
  }
  if (abs(det(J_vp))<1e-10) {
    invJ_vp <- ginv(J_vp)
  } else {
    invJ_vp <- solve(J_vp)
  }
  dif <- 1
  differ <- c()  ## added to record the difference
  iter <- 0
  res <- list()
  res$beta <- rep(NA, p_b)
  res$alpha <- rep(NA, p_a)
  res$variance <- matrix(0, nrow=p_t, ncol=p_t)
  res$gamMat <- gamMat
  res$varphiMat <- varphiMat
  res$gam.variance <- gamMat.var
  res$varphi.variance <- varphiMat.var
  res$convergence <- FALSE
  res$iteration = 0

  while(iter<maxit & dif >tol) {
    beta_est.o <- beta
    alpha_est.o <- alpha
    theta_est.o <- theta
    U = rep(0, p_t)
    M = matrix(0, nrow=p_t, ncol=p_t)
    Lambda_g <- matrix(0, nrow=p_t, ncol=p_g)
    Lambda_vp <- matrix(0, nrow=p_t, ncol=p_vp)
    U_i <- rep(0, p_t)
    M_i = matrix(0, nrow=p_t, ncol=p_t)
    Lambda_g_i <- matrix(0, nrow=p_t, ncol=p_g)
    Lambda_vp_i <- matrix(0, nrow=p_t, ncol=p_vp)
    U_i.all <- matrix(0, nrow=p_t, ncol=n)
    for (i in 1:n) {
      #Ytilde_i <- t(Ytilde.Mat.spl[[i]])
      DM_i <- DM.spl[[i]]
      Ytilde_i <- Ytilde[(i-1)*(m*K) + 1:(m*K)]
      assocDM_i <- assocDM[(i-1)*nrow_assocDM_i +
                             1:nrow_assocDM_i, ]
      Ztilde_i <- Ztilde[(i-1)*nrow_assocDM_i +
                           1:nrow_assocDM_i]
      Ytilde_i_Mat_ext <- Ytilde.Mat.ext.spl[[i]]
      invPast_ij_all <- invPast.spl[[i]]
      P_ij_all <- P.spl[[i]]
      wgt_i <- wgt[i, ]
      Xtilde_i_Mat_ext <- Xtilde.ext.Mat.spl[[i]]
      invGast_ij_all <- invGast.spl[[i]]
      G_ij_all <- G.spl[[i]]
      delta_i <- delta[(i-1)*m + 1:m]
      ##  ---- Pass to C  ---- ##
      ## Be careful of the dimensions
      mgee2v_i <- .C("Cgetmgee2v_i",
                     as.double(DM_i),
                     as.double(Ytilde_i),
                     as.double(assocDM_i),
                     as.double(Ztilde_i),
                     as.double(Ytilde_i_Mat_ext),
                     as.double(invPast_ij_all),
                     as.double(P_ij_all),
                     as.double(L_ij),
                     as.double(wgt_i),
                     as.double(X_i.ENUM),
                     as.double(XDM_ENUM),
                     as.double(Xtilde_i_Mat_ext),
                     as.double(invGast_ij_all),
                     as.double(G_ij_all),
                     as.double(L_x_ij),
                     as.double(delta_i),
                     as.integer(m),
                     as.integer(K),
                     as.integer(K_x),
                     as.integer(p_b),
                     as.integer(p_a),
                     as.integer(p_g),
                     as.integer(p_vp),
                     as.double(beta),
                     as.double(alpha),
                     U_i = as.double( U_i ),
                     M_i = as.double(M_i),
                     Lambda_g_i = as.double(Lambda_g_i),
                     Lambda_vp_i = as.double(Lambda_vp_i)
      )

      U <- U + (mgee2v_i$U_i)/n
      M <- M + matrix(mgee2v_i$M_i, nrow=p_t, byrow=FALSE)/n
      Lambda_g <- Lambda_g + (mgee2v_i$Lambda_g_i)/n
      Lambda_vp <- Lambda_vp + (mgee2v_i$Lambda_vp_i)/n
      U_i.all[, i] <- mgee2v_i$U_i
    }
    U_beta <- U[1:p_b]
    U_alpha <- U[(p_b+1):p_t]
    M1 <- M[1:p_b, 1:p_b]
    M2 <- as.matrix( M[(p_b+1):p_t, (p_b+1):p_t] )
    if (any(is.na(M1)) | any(is.na(M2))) {
      res$iteration <- iter
      return(res)
    } else {
      if ((abs(det(M1)) < 1e-15) | (abs(det(M2)) < 1e-15)) {
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
  Omega_i.all <- U_i.all - (Lambda_g%*%invJ_g%*%Q_g_i.all +
                              Lambda_vp%*%invJ_vp%*%Q_vp_i.all)
  Lambda <- cbind(Lambda_g, Lambda_vp)

  Sigma <- (Omega_i.all%*%t(Omega_i.all))/n
  Gamma <- M

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
  res$iteration <- iter
  res$differ <- differ ## return the difference
  names(res$beta) <- beta_nam
  names(res$alpha) <- alpha_nam
  dimnames(res$variance) <- list(c(beta_nam, alpha_nam),
                                 c(beta_nam, alpha_nam))
  # Yuliang: include call in the output
  res$call <- match.call()
  class(res) <- "mgee2"
  return(res)
}### end of mgee2v()
