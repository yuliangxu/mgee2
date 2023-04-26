## Get block matrix
blockmat <- function(mat1, mat2)
{
    nrows <- nrow(mat1) + nrow(mat2)
    ncols <- ncol(mat1) + ncol(mat2)
    outmat <- matrix(0, nrows, ncols)
    outmat[1:nrow(mat1), 1:ncol(mat1)] <- mat1
    outmat[(1+nrow(mat1)):nrows, (1+ncol(mat1)):ncols] <- mat2
    return(outmat)
}

## Get unbiased surrogate vector ${\bf X}_ij^{\ast}$
## from observed surrogate vector
get.mcPrMat_ij <- function(S_ij, L_ij, gamMat, K, adjacent=FALSE) {
    if (!adjacent) {
        EXP.Mat <- cbind(rep(1,1+K), matrix(
           as.vector(exp(rbind(L_ij)%*%gamMat)), ncol=K, byrow=TRUE))
        trPrMat <- EXP.Mat/apply(EXP.Mat, 1, FUN=sum)
    }  else {
        print("Need modification")
    }
    return(trPrMat)
}
get.UnbSurr_ij <- function(S_ij, L_ij, gamMat, K, adjacent=FALSE) {
    if (!adjacent) {
        EXP.Mat <- cbind(rep(1,1+K), matrix(
           as.vector(exp(rbind(L_ij)%*%gamMat)), ncol=K, byrow=TRUE))
        trPrMat <- EXP.Mat/apply(EXP.Mat, 1, FUN=sum)
    }  else {
        print("Need modification")
    }
    pi_ij0 <- trPrMat[1, 2:(K+1)]
    Q_ij <- t(trPrMat[2:(K+1), 2:(K+1)]) - pi_ij0
    Y_ij.ast <- solve(Q_ij)%*%(S_ij-pi_ij0)
    return(Y_ij.ast)
}





## U_i
get.U_i <- function(m_i, Y_i, DM_i, beta, j1.indx, k1.indx,
       j2.indx, k2.indx, u_i, alpha, K) {
    p_b <- length(beta)
    p_a <- length(alpha)
    p_t <- p_b+p_a
    indx.tab <- cbind(j1.indx, k1.indx, j2.indx, k2.indx)
    indx.tab.t <- t(indx.tab)
    ## cumulative probs, marginal probs, and their derivatives
    lambda_i <- as.vector(expit(DM_i%*%beta))
    lambda_i.Mat <- matrix(lambda_i, ncol=K, byrow=T)
    mu_i.Mat <- lambda_i.Mat - cbind(lambda_i.Mat[, 2:K], 0)
    mu_i <- as.vector(t(mu_i.Mat))
    dlambda_i_dbetaT <- lambda_i*(1-lambda_i)*DM_i
    dlambda_i.2_dbetaT <- matrix(0, nrow=length(lambda_i),ncol=p_b)
    indx.2 <- kronecker( K*(0:(m_i-1)),2:K, "+")
    dlambda_i.2_dbetaT[indx.2-1, ] <- dlambda_i_dbetaT[indx.2, ]
    dmu_i_dbetaT <- dlambda_i_dbetaT - dlambda_i.2_dbetaT
    D_1i <- t(dmu_i_dbetaT)
    var_i <- mu_i*(1-mu_i)

    ## bivariate cumulatives, and their derivatives
    rows.Mat <- matrix(1:(K*K), ncol=K, byrow=T)
    psi_i <- as.vector(exp(u_i%*%alpha))
    lambda_i_j1k1 <- lambda_i.Mat[cbind(j1.indx, k1.indx)]
    lambda_i_j2k2 <- lambda_i.Mat[cbind(j2.indx, k2.indx)]
    a_i <- 1-(1-psi_i)*(lambda_i_j1k1+lambda_i_j2k2)
    d_i <- a_i^2-4*(psi_i-1)*psi_i*lambda_i_j1k1*lambda_i_j2k2
    f_i <- a_i - sqrt(d_i)
    g_i <- 2*(psi_i-1)
    zeta_i <- f_i/g_i
    zeta_i[psi_i==1] <- (lambda_i_j1k1*lambda_i_j2k2)[psi_i==1]

    dpsi_i_dalphaT <- psi_i*u_i
    da_i_dpsi_i <- lambda_i_j1k1 + lambda_i_j2k2
    dd_i_dpsi_i <- ( 2*a_i*da_i_dpsi_i -4*(2*psi_i-1)*
        (lambda_i_j1k1*lambda_i_j2k2) )
    df_i_dpsi_i <- da_i_dpsi_i - 1/(2*sqrt(d_i))*dd_i_dpsi_i
    dg_i_dpsi_i <- 2
    dzeta_i_dpsi_i <- df_i_dpsi_i/g_i - dg_i_dpsi_i*f_i/g_i^2
    dzeta_i_dalphaT <- dzeta_i_dpsi_i*dpsi_i_dalphaT

    dlambda_iprod_dbetaT <- ( dlambda_i_dbetaT[(j1.indx-1)*K+k1.indx, ] *
        lambda_i_j2k2 + dlambda_i_dbetaT[(j2.indx-1)*K+k2.indx, ] *
        lambda_i_j1k1 )
    da_i_dbetaT <- ( -(1-psi_i)*(dlambda_i_dbetaT[(j1.indx-1)*K+k1.indx, ] +
        dlambda_i_dbetaT[(j2.indx-1)*K+k2.indx, ]) )
    dd_i_dbetaT <- ( 2*a_i*da_i_dbetaT -
        4*(psi_i-1)*psi_i*dlambda_iprod_dbetaT )
    df_i_dbetaT <- da_i_dbetaT  - 1/(2*sqrt(d_i))*dd_i_dbetaT
    dzeta_i_dbetaT <- df_i_dbetaT/g_i


    ## Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, V_2i
    xi_i <- rep(0, length(zeta_i))
    Z_i <- rep(0, length(zeta_i))
    V_1i <- diag(length(Y_i))
    V_2i <- diag(length(zeta_i))
    inv.V_2i <- diag(length(xi_i))
    dxi_i_dalphaT <- matrix(0, nrow=length(zeta_i), ncol=p_a)
    dxi_i_dbetaT <- matrix(0, nrow=length(zeta_i), ncol=p_b)
    Y_i.Mat <- matrix(Y_i, ncol=K, byrow=TRUE)
    mu_i.Mat <- matrix(mu_i, ncol=K, byrow=TRUE)
    for (j1 in 1:m_i) {
        Y_ij1 <- Y_i.Mat[j1,]
        mu_ij1 <- mu_i.Mat[j1,]
        V_1i[K*(j1-1)+1:K, K*(j1-1)+1:K] <- ( diag(mu_ij1, nrow=K) -
            outer(mu_ij1, mu_ij1) )
        if (j1!=m_i) {
            for (j2 in (j1+1):m_i) {
                Y_ij2 <- Y_i.Mat[j2,]
                mu_ij2 <- mu_i.Mat[j2,]
                j1j2.indx <- j1.indx==j1 & j2.indx==j2

                zeta_i.j1j2 <- zeta_i[j1j2.indx]
                dzeta_i.j1j2_dalphaT <- as.matrix(
                    dzeta_i_dalphaT[j1j2.indx, ], ncol=p_a)
                dzeta_i.j1j2_dbetaT <- dzeta_i_dbetaT[j1j2.indx, ]
                xi_i.j1j2 <- zeta_i.j1j2
                dxi_i.j1j2_dalphaT <- dzeta_i.j1j2_dalphaT
                dxi_i.j1j2_dbetaT <- dzeta_i.j1j2_dbetaT
                ## case i) k_1<K & k_2<K
                rows.indx <- as.vector(t(rows.Mat[cbind(1:(K-1),1:(K-1))]))
                rows.k1p1.indx <- as.vector(t(rows.Mat[cbind(2:K,1:(K-1))]))
                rows.k2p1.indx <- as.vector(t(rows.Mat[cbind(1:(K-1),2:K)]))
                rows.k1p1.k2p1.indx <- as.vector(t(rows.Mat[cbind(2:K,2:K)]))
                xi_i.j1j2[rows.indx] <- ( zeta_i.j1j2[rows.indx] -
                    zeta_i.j1j2[rows.k1p1.indx] -
                    zeta_i.j1j2[rows.k2p1.indx] +
                    zeta_i.j1j2[rows.k1p1.k2p1.indx]  )
                dxi_i.j1j2_dalphaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dalphaT[rows.indx,] -
                    dzeta_i.j1j2_dalphaT[rows.k1p1.indx,] -
                    dzeta_i.j1j2_dalphaT[rows.k2p1.indx,] +
                    dzeta_i.j1j2_dalphaT[rows.k1p1.k2p1.indx,] )
                dxi_i.j1j2_dbetaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dbetaT[rows.indx,] -
                    dzeta_i.j1j2_dbetaT[rows.k1p1.indx,] -
                    dzeta_i.j1j2_dbetaT[rows.k2p1.indx,] +
                    dzeta_i.j1j2_dbetaT[rows.k1p1.k2p1.indx,] )
                ## case ii) k_1==K & k_2<K
                rows.indx <- as.vector(t(rows.Mat[K,1:(K-1)]))
                rows.k2p1.indx <- as.vector(t(rows.Mat[K,2:K]))
                xi_i.j1j2[rows.indx] <- ( zeta_i.j1j2[rows.indx] -
                    zeta_i.j1j2[rows.k2p1.indx] )
                dxi_i.j1j2_dalphaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dalphaT[rows.indx,] -
                    dzeta_i.j1j2_dalphaT[rows.k2p1.indx,] )
                dxi_i.j1j2_dbetaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dbetaT[rows.indx,] -
                    dzeta_i.j1j2_dbetaT[rows.k2p1.indx,] )
                ## case iii) k_1<K & k_2==K
                rows.indx <- as.vector(t(rows.Mat[1:(K-1), K]))
                rows.k1p1.indx <- as.vector(t(rows.Mat[2:K,K]))
                xi_i.j1j2[rows.indx] <- ( zeta_i.j1j2[rows.indx] -
                    zeta_i.j1j2[rows.k1p1.indx] )
                dxi_i.j1j2_dalphaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dalphaT[rows.indx,] -
                    dzeta_i.j1j2_dalphaT[rows.k1p1.indx,] )
                dxi_i.j1j2_dbetaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dbetaT[rows.indx,] -
                    dzeta_i.j1j2_dbetaT[rows.k1p1.indx,] )
                ## case iv) k_1==K & k_2==K
                rows.indx <- as.vector(t(rows.Mat[K, K]))
                xi_i.j1j2[rows.indx] <- zeta_i.j1j2[rows.indx]
                dxi_i.j1j2_dalphaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dalphaT[rows.indx,] )
                dxi_i.j1j2_dbetaT[rows.indx, ] <- (
                    dzeta_i.j1j2_dbetaT[rows.indx,] )

                ## Combine results
                Z_i[j1j2.indx] <- kronecker(Y_ij1, Y_ij2)
                xi_i[j1j2.indx] <- xi_i.j1j2
                dxi_i_dalphaT[j1j2.indx, ] <- dxi_i.j1j2_dalphaT
                dxi_i_dbetaT[j1j2.indx, ] <- dxi_i.j1j2_dbetaT

                xi_i.j1j2.Mat <- matrix(xi_i.j1j2, ncol=K, byrow=TRUE)
                cov_1ij1j2 <- xi_i.j1j2.Mat - outer(mu_ij1, mu_ij2)
                V_1i[K*(j1-1)+1:K, K*(j2-1)+1:K] <- cov_1ij1j2
                V_1i[K*(j2-1)+1:K, K*(j1-1)+1:K] <-  t(cov_1ij1j2)

                V_2ij1j2 <- diag(xi_i.j1j2) - outer(xi_i.j1j2, xi_i.j1j2)
                V_2i[j1j2.indx, j1j2.indx] <- V_2ij1j2

                if (abs(det(as.matrix(V_2ij1j2)))<1e-20) {
                    message("V_2ij1j2 is singular")
                    inv.V_2i[j1j2.indx, j1j2.indx] <- ginv(V_2ij1j2)
                }  else {
                    inv.V_2i[j1j2.indx, j1j2.indx] <- solve(V_2ij1j2)
                }
            }
        }
    }
    if (abs(det(as.matrix(V_1i)))<1e-20) {
        message("V_1i is singular")
        inv.V_1i <- ginv(V_1i)
    } else    inv.V_1i <- solve(V_1i)
    D_2i <- t(dxi_i_dalphaT)
    U_1i <- D_1i%*%inv.V_1i%*%(Y_i-mu_i)
    U_2i <- D_2i%*%inv.V_2i%*%(Z_i-xi_i)

    ## --------------------------
    ## Get M
    ## --------------------------
    M_1i <- -D_1i%*%inv.V_1i%*%t(D_1i)
    M_12i <- matrix(0, nrow=p_b, ncol=p_a)
    M_2i <- -D_2i%*%inv.V_2i%*%t(D_2i)
    M_21i <- -D_2i%*%inv.V_2i%*%dxi_i_dbetaT
    M_i <- matrix(0, nrow=p_t, ncol=p_t)
    M_i[1:p_b, 1:p_b] <- M_1i
    M_i[(p_b+1):p_t, 1:p_b] <- M_21i
    M_i[(p_b+1):p_t, (p_b+1):p_t] <- M_2i

    U_i_unls <- c(U_1i, U_2i, M_i)
    return(U_i_unls)
}


get.xi_i <- function(j1.indx, k1.indx, j2.indx, k2.indx, zeta_i, K) {
    indx.tab <- cbind(j1.indx, k1.indx, j2.indx, k2.indx)
    indx.tab.t <- t(indx.tab)
    xi_i <- NULL
    for (j in 1:length(zeta_i)) {
        if((k1.indx[j]<K) && (k2.indx[j]<K))  {
            vec <- indx.tab.t[, j]
            vec.k1plus1 <- vec + c(0,1,0,0)
            jp1 <- apply(indx.tab.t==vec.k1plus1, 2, prod)==1
            vec.k2plus1 <- vec + c(0,0,0,1)
            jp2 <- apply(indx.tab.t==vec.k2plus1, 2, prod)==1
            vec.k1k2plus1 <- vec + c(0,1,0,1)
            jp3 <- apply(indx.tab.t==vec.k1k2plus1, 2, prod)==1
            xi_i <- c(xi_i, zeta_i[j]-zeta_i[jp1]
                -zeta_i[jp2]+zeta_i[jp3])
        } else {
            if ((k1.indx[j]==K) && (k2.indx[j]<K))  {
                vec <- indx.tab.t[, j]
                vec.k2plus1 <- vec + c(0,0,0,1)
                jp2 <- apply(indx.tab.t==vec.k2plus1, 2, prod)==1
                xi_i <- c(xi_i, zeta_i[j]-zeta_i[jp2])
            } else {
                if ((k1.indx[j]<K) && (k2.indx[j]==K))  {
                    vec <- indx.tab.t[, j]
                    vec.k1plus1 <- vec + c(0,1,0,0)
                    jp1 <- apply(indx.tab.t==vec.k1plus1, 2, prod)==1
                    xi_i <- c(xi_i, zeta_i[j]-zeta_i[jp1])
                } else {
                    xi_i <- c(xi_i, zeta_i[j])
                }
            }
        }
    }
    return(xi_i)
}


get.xi_i.ext <- function(j1.indx.ext, k1.indx.ext, j2.indx.ext,
          k2.indx.ext, mu_i.Mat, indx.tab.t, xi_i, K) {
    indx.tab.ext <- cbind(j1.indx.ext, k1.indx.ext,
        j2.indx.ext, k2.indx.ext)
    indx.tab.ext.t <- t(indx.tab.ext)

    xi_i.ext <- rep(0, ncol(indx.tab.ext.t))

    for (j in 1:ncol(indx.tab.ext.t)) {
        j1 <- indx.tab.ext.t[1,j]
        k1 <- indx.tab.ext.t[2,j]
        j2 <- indx.tab.ext.t[3,j]
        k2 <- indx.tab.ext.t[4,j]
        if ((k1==0) & (k2!=0)) {
            xi_i.ext[j] <- ( mu_i.Mat[j2,k2] -
               sum(xi_i[indx.tab.t[1,]==j1 &
                  indx.tab.t[3,]==j2 & indx.tab.t[4,]==k2]) )
        } else {
            if ((k1!=0) & (k2==0)) {
                xi_i.ext[j] <- ( mu_i.Mat[j1,k1] -
                    sum(xi_i[indx.tab.t[1,]==j1 &
                    indx.tab.t[2,]==k1 & indx.tab.t[3,]==j2 ]) )
            } else {
                if ((k1==0) & (k2==0)) {
                    xi_i.ext[j] <- ( 1-sum(mu_i.Mat[c(j1,j2),])+
                        sum(xi_i[indx.tab.t[1,]==j1 &
                        indx.tab.t[3,]==j2]) )
                } else {
                    xi_i.ext[j] <- xi_i[indx.tab.t[1,]==j1 &
                        indx.tab.t[2,]==k1 & indx.tab.t[3,]==j2 &
                        indx.tab.t[4,]==k2]
                }
            }
        }
    }
    return(xi_i.ext)
}


## Get polytomous response vector from a categorical response
getPolyResp <- function(Y.c, K, ref="0") {
    Resp.Mat <- matrix(0, nrow=length(Y.c), ncol=1+K)
    Resp.Mat[cbind(1:length(Y.c),
        as.numeric(levels(Y.c))[Y.c]+1)] <- rep(1,length(Y.c))
    Vec <- as.vector(t(Resp.Mat[,2:(1+K)]))
    return(Vec)
}

## -----------------------------------------------------------------
##   Function to simulate true correlated binary data
## -----------------------------------------------------------------
get.truedat <- function(Yform, n, m, K, K_x, OR.str="exchangeable",
        regcoef1, regcoef2, alpha)
{
    ##------------------------
    ## Generate the covariates
    ##------------------------
    ID <- as.vector(outer(rep(1,m), 1:n))
    ID.treat1 <- sample(1:n, size=round(n/2))
    treat.indx <- numeric(n)
    treat.indx[ID.treat1] <- rep(1, length(ID.treat1))
    treatment <- as.vector(outer(rep(1,m), treat.indx))
    treatment <- as.factor(treatment)
    visit <- rep(1:m, n)
    visit <- as.factor(visit)


    ## X~multinomial  (time-variant)
    X <- as.vector( sapply(rep(1, n*m),
       FUN=sample, x=0:2, prob=c(.5, .3, .2)) )
    X <- as.factor(X)

    ## -------------------------
    ## Generate Y
    ## -------------------------
    Y <- rep(NA, length(ID))
    Pr.Y <- rep(NA, length(ID))
    truedat <- data.frame(ID=ID, Y=Y, X=X, treatment=treatment, visit=visit)

    origDM <- getDM(formula=Yform, data=truedat)

    lambda.1 <- as.vector(expit(origDM%*%regcoef1))
    lambda.2 <- as.vector(expit(origDM%*%regcoef2))
    mu.0 <- 1-lambda.1
    mu.1 <- lambda.1 - lambda.2
    mu.2 <- lambda.2


    j1.indx <- NULL
    k1.indx <- NULL
    j2.indx <- NULL
    k2.indx <- NULL
    j1.indx.ext <- NULL
    k1.indx.ext <- NULL
    j2.indx.ext <- NULL
    k2.indx.ext <- NULL
    for(j1 in 1:(m-1)) {
        j1.indx <- c(j1.indx, rep(j1, (m-j1)*K^2))
        k1.indx <- c(k1.indx, kronecker(1:K, rep(1, (m-j1)*K)))
        j2.indx <- c(j2.indx, kronecker(rep(1,K),
            kronecker((j1+1):m, rep(1, K))))
        k2.indx <- c(k2.indx, kronecker(rep(1, (m-j1)*K), 1:K))

        j1.indx.ext <- c(j1.indx.ext, rep(j1, (m-j1)*(1+K)^2))
        k1.indx.ext <- c(k1.indx.ext,
            kronecker(0:K, rep(1, (m-j1)*(1+K))))
        j2.indx.ext <- c(j2.indx.ext, kronecker(rep(1,(1+K)),
            kronecker((j1+1):m, rep(1, 1+K))))
        k2.indx.ext <- c(k2.indx.ext,
            kronecker(rep(1, (m-j1)*(1+K)), 0:K))
    }
    indx.tab <- cbind(j1.indx, k1.indx, j2.indx, k2.indx)
    indx.tab.t <- t(indx.tab)

    for (i in 1:n) {
        indx <- ID==i

        mu_i.0 <- mu.0[indx]
        mu_i.1 <- mu.1[indx]
        mu_i.2 <- mu.2[indx]
        mu_i.Mat <- cbind(mu_i.1, mu_i.2)
        lambda_i.Mat <- cbind(lambda.1[indx], lambda.2[indx])
        lambda_i <- as.vector(t(lambda_i.Mat))

        if (OR.str=="exchangeable") {
            u_i <- matrix(1, nrow=length(j1.indx), ncol=1)
        } else {
            if (OR.str=="log-linear") {
                u_i <- cbind(rep(1, length(j1.indx)),
                    (as.numeric(k1.indx==2) + as.numeric(k2.indx==2)),
                    as.numeric((k1.indx==2)&(k2.indx==2)))
            }
        }
        psi_i <- as.vector(exp(u_i%*%alpha))

        lambda_i_j1k1 <- lambda_i.Mat[cbind(j1.indx, k1.indx)]
        lambda_i_j2k2 <- lambda_i.Mat[cbind(j2.indx, k2.indx)]
        a_i <- 1-(1-psi_i)*(lambda_i_j1k1+lambda_i_j2k2)
        d_i <- a_i^2-4*(psi_i-1)*psi_i*lambda_i_j1k1*lambda_i_j2k2
        zeta_i <- (a_i - sqrt(d_i))/(2*(psi_i-1))
        zeta_i[psi_i==1] <- (lambda_i_j1k1*lambda_i_j2k2)[psi_i==1]
        xi_i <- get.xi_i(j1.indx=j1.indx, k1.indx=k1.indx,
            j2.indx=j2.indx, k2.indx=k2.indx, zeta_i=zeta_i, K=K)

        ## following quantities include level 0 (.ext mean extended)
        mu_i.ext.Mat <- cbind(mu_i.0, mu_i.1, mu_i.2)

        xi_i.ext <- get.xi_i.ext(j1.indx.ext=j1.indx.ext,
            k1.indx.ext=k1.indx.ext, j2.indx.ext=j2.indx.ext,
            k2.indx.ext=k2.indx.ext, mu_i.Mat=mu_i.Mat,
            indx.tab.t=indx.tab.t, xi_i=xi_i, K=K)
        mu_i.ext.j1k1 <- mu_i.ext.Mat[cbind(j1.indx.ext,k1.indx.ext+1)]
        var_i.ext.j1k1 <- (1-mu_i.ext.j1k1)*mu_i.ext.j1k1
        mu_i.ext.j2k2 <- mu_i.ext.Mat[cbind(j2.indx.ext,k2.indx.ext+1)]
        var_i.ext.j2k2 <- (1-mu_i.ext.j2k2)*mu_i.ext.j2k2
        rho_i.ext <- (xi_i.ext - mu_i.ext.j1k1*mu_i.ext.j2k2)/
            sqrt(var_i.ext.j1k1*var_i.ext.j2k2)

        Y_i.ENUM <- NULL
        for (j in 1:m) {
           Y_i.ENUM <- rbind(Y_i.ENUM,
               rep( kronecker(0:K, rep(1, (1+K)^(m-j))), (1+K)^(j-1) ))
        }
        Pr_i.ENUM <- rep(0, ncol(Y_i.ENUM))
        y.j1.indx <- NULL
        y.j2.indx <- NULL
        for(j in 1:(m-1)) {
            y.j1.indx <- c(y.j1.indx, rep(j, length((j+1):m)))
            y.j2.indx <- c(y.j2.indx, (j+1):m)
        }
        for (enum in 1:ncol(Y_i.ENUM)) {
            mu_i.ext.enum <- mu_i.ext.Mat[cbind(1:m,Y_i.ENUM[, enum]+1)]
            sum_2 <- 0
            for (j1 in 1:(m-1)) {
                k1 <- Y_i.ENUM[j1, enum]
                mu_i.ext.j1k1.enum <- mu_i.ext.Mat[j1 , k1+1]
                var_i.ext.j1k1.enum <- mu_i.ext.j1k1.enum*(1
                    -mu_i.ext.j1k1.enum)
                for (j2 in (j1+1):m) {
                    k2 <- Y_i.ENUM[j2,enum]
                    mu_i.ext.j2k2.enum <- mu_i.ext.Mat[j2 , k2+1]
                    var_i.ext.j2k2.enum <- mu_i.ext.j2k2.enum*(1
                        -mu_i.ext.j2k2.enum)
                    rho_i.ext.j1k1.j2k2.enum <- rho_i.ext[j1.indx.ext==j1 &
                        k1.indx.ext==k1 & j2.indx.ext==j2 & k2.indx.ext==k2]
                    sum_2 <- ( sum_2+rho_i.ext.j1k1.j2k2.enum*(1-
                        mu_i.ext.j1k1.enum)*(1-mu_i.ext.j2k2.enum)/sqrt(
                        var_i.ext.j1k1.enum*var_i.ext.j2k2.enum) )
                }
            }
            Pr_i.ENUM[enum] <- prod(mu_i.ext.enum)*(1+sum_2)
        }
        Y_i <- Y_i.ENUM[, sample(size=1, 1:ncol(Y_i.ENUM),
            prob=Pr_i.ENUM)]
        Y[indx] <- Y_i
    }
    truedat$Y <- as.factor(Y)

    return(truedat)
}


get.Y.mc.probs_ij <- function(Y_ij, L_ij, gamMat, K) {
    mc.probs_ij <- c(1, as.vector(exp(rbind(L_ij)%*%
        gamMat))[K*Y_ij+1:K])
    mc.probs_ij <- mc.probs_ij/sum(mc.probs_ij)
    return(mc.probs_ij)
}
get.X.mc.probs_ij <- function(X_ij, L.x_ij, varphiMat, K_x) {
    mc.probs_ij <- c(1, as.vector(exp(rbind(L.x_ij)
        %*%varphiMat))[K_x*X_ij+1:K_x])
    mc.probs_ij <- mc.probs_ij/sum(mc.probs_ij)
    return(mc.probs_ij)
}

get.obsdat <- function(n, Y.mcform, X.mcform, gamMat, varphiMat,
        truedat, validation="random30")
{
    obsdat <- data.frame(truedat, S=NA, W=NA)
    ID <- truedat$ID
    if (is.numeric(truedat$Y)) {
        Y <- truedat$Y
    } else {
        Y <- as.numeric(levels(truedat$Y))[truedat$Y]
    }
    if (is.numeric(truedat$X)) {
        X <- truedat$X
    } else {
        X <- as.numeric(levels(truedat$X))[truedat$X]
    }

    K <- max(Y)
    K_x <- max(X)

    ## Get surrogate response
    Y.spl <- split(Y, f=1:length(Y))
    Y.mc.DV.spl <- split(getDM(formula=Y.mcform, data=obsdat),
        f=1:length(Y))
    S.probs.all.Mat <- matrix(mapply(Y.spl, Y.mc.DV.spl,
        FUN=get.Y.mc.probs_ij,  SIMPLIFY = TRUE,
        MoreArgs=list(gamMat=gamMat, K=K)),
        nrow=1+K, byrow=FALSE)
    S <- apply(S.probs.all.Mat, 2, FUN=sample,
        size=1, x=0:K, replace=FALSE)
    obsdat$S <- as.factor(S)

    ## Get surrogate covariate
    X.spl <- split(X, f=1:length(X))
    X.mc.DV.spl <- split(getDM(formula=X.mcform, data=obsdat),
        f=1:length(X))
    W.probs.all.Mat <- matrix(mapply(X.spl, X.mc.DV.spl,
        FUN=get.X.mc.probs_ij,  SIMPLIFY = TRUE,
        MoreArgs=list(varphiMat=varphiMat, K_x=K_x)),
        nrow=1+K_x, byrow=FALSE)
    W <- apply(W.probs.all.Mat, 2, FUN=sample,
        size=1, x=0:K_x, replace=FALSE)
    obsdat$W <- as.factor(W)

    if(is.null(validation)) {
        obsdat$Y <- rep(NA, nrow(obsdat))
        obsdat$X <- rep(NA, nrow(obsdat))
    } else{
        if (validation=="random30") {
            obsdat$Y <- truedat$Y
            obsdat$X <- truedat$X
            obsdat$delta <- rbinom(n=nrow(obsdat), size=1, prob=0.3)
            obsdat$Y[obsdat$delta==0] <- NA
            obsdat$X[obsdat$delta==0] <- NA
        } else {
            if(validation=="nonrandom40") {
                obsdat$Y <- truedat$Y
                obsdat$X <- truedat$X
                probs <- c(0.45, 0.27, 0.18, 0.75, 0.45, 0.3)[as.numeric(obsdat$treatment=="1") * 3 + as.numeric(obsdat$visit)]
                obsdat$delta <- rbinom(n=nrow(obsdat), size=1, prob=probs)
                obsdat$Y[obsdat$delta==0] <- NA
                obsdat$X[obsdat$delta==0] <- NA
            } else {
                if (validation=="random20") {
                   obsdat$Y <- truedat$Y
                   obsdat$X <- truedat$X
                   obsdat$delta <- rbinom(n=nrow(obsdat), size=1, prob=0.2)
                   obsdat$Y[obsdat$delta==0] <- NA
                   obsdat$X[obsdat$delta==0] <- NA
                }
            }
        }
    }
    return(obsdat)
}




## Get the corrected estimating functions
get.U_i.ast <- function(m_i, Y.ast_i, X.ast.ext.Mat_i,
       DM_i, beta, j1.indx, k1.indx, j2.indx, k2.indx,
       u_i, alpha, K, K_x, X_i.ENUM, e.Mat) {
    U_i.ast <- 0
    for (enum in 1:ncol(X_i.ENUM)) {
        DM_i.enum <- DM_i
        DM_i.enum[,(K+1):(K+K_x)] <- e.Mat[
            kronecker(X_i.ENUM[,enum]+1, rep(1,K_x)), ]
        U_i.enum <- get.U_i(m_i=m_i,Y_i=Y.ast_i, DM_i=DM_i.enum,
            beta=beta, j1.indx=j1.indx, k1.indx=k1.indx,
            j2.indx=j2.indx, k2.indx=k2.indx,
            u_i=u_i, alpha=alpha, K=K)
        X.ast.ext_i.enum <- X.ast.ext.Mat_i[
            cbind(1:m_i,X_i.ENUM[,enum]+1)]
        U_i.ast <- U_i.ast + U_i.enum * prod(X.ast.ext_i.enum)
    }
    return(U_i.ast)
}

get.DM.vec.ENUM <- function(i, X_i.ENUM, DM.spl, K, K_x, e.Mat) {
            DM_i <- DM.spl[[i]]
            l.DM <- prod(dim(DM_i))
            n.ENUM <- ncol(X_i.ENUM)
            DM_i.vec.ENUM <- numeric(l.DM*n.ENUM)
            for (enum in 1:n.ENUM) {
                DM_i[, (K+1):(K+K_x)] <- e.Mat[kronecker(
                    X_i.ENUM[,enum]+1, rep(1,K_x)), ]
                DM_i.vec.ENUM[l.DM*(enum-1)+1:l.DM] <- as.vector(DM_i)
            }
            return(DM_i.vec.ENUM)
}

get.prod.X_i.ast.enum <- function(X.ast.ext.Mat_i, X_i.enum, m_i) {
        return(prod(X.ast.ext.Mat_i[cbind(1:m_i, X_i.enum+1)]))
}

get.prod.X_i.ast.ENUM <- function(X.ast.ext.Mat_i, X_i.ENUM, m_i) {
        prod.X_i.ast.ENUM <- apply(X_i.ENUM, 2,
            FUN=get.prod.X_i.ast.enum,
            X.ast.ext.Mat_i=X.ast.ext.Mat_i, m_i=m_i)
        return(prod.X_i.ast.ENUM)
}

get.U.ast_i <- function(m_i, Y.ast_i, pX.DM_i.vec.ENUM,
       beta, j1.indx, k1.indx, j2.indx, k2.indx,
       u_i, alpha, K, K_x, n.ENUM) {
    prodX.ast_i.ENUM <- pX.DM_i.vec.ENUM[1:n.ENUM]
    DM_i.vec.ENUM <- matrix(pX.DM_i.vec.ENUM[-(1:n.ENUM)],
        ncol=n.ENUM, byrow=FALSE)
    U_i.ENUM <- apply(DM_i.vec.ENUM, 2,
        FUN=get.U_i.vec, m_i=m_i, Y_i=Y.ast_i, beta=beta,
        j1.indx=j1.indx, k1.indx=k1.indx,
        j2.indx=j2.indx, k2.indx=k2.indx,
        u_i=u_i, alpha=alpha, K=K)
    U_i.ast <- apply(t(U_i.ENUM)*prodX.ast_i.ENUM, 2, sum)

    return(U_i.ast)
}

get.U_i.vec <- function(m_i, Y_i, DM_i.vec, beta, j1.indx, k1.indx,
       j2.indx, k2.indx, u_i, alpha, K)  {
    DM_i <- matrix(DM_i.vec, ncol=length(beta), byrow=FALSE)
    return(get.U_i(m_i=m_i, Y_i=Y_i, DM_i=DM_i, beta=beta,
       j1.indx=j1.indx, k1.indx=k1.indx,
       j2.indx=j2.indx, k2.indx=k2.indx, u_i=u_i,
       alpha=alpha, K=K))
}




get_Z_i <- function(Y_i, m, K)
{
    nel <- (K^2) * (m * (m - 1))/2
    Y_i_MAT <- matrix(Y_i, nrow=K)
    Z_i <- NULL
    for (j1 in 1:(m-1))
    {
        for (j2 in (j1+1):m)
        {
            Z_i <- c(Z_i, kronecker(
                Y_i_MAT[, j1], Y_i_MAT[, j2]))
        }
    }
    return(Z_i)
}

get_Z <- function(Y, n, m, K)
{
    len_Zi <- (K^2) * (m * (m - 1))/2
    Z <- numeric(len_Zi * n)
    for (i in 1:n)
    {
        Y_i <- Y[(i-1)*(m*K) + 1:(m*K)]
        Z[(i-1)*len_Zi + 1:len_Zi] <- get_Z_i(Y_i, m, K)
    }
    return(Z)
}




## Inference about the error processes
get.errpar <- function(obsform, id, obsdat, OR.str="exchangeable",
     Y.mcform, X.mcform, p_b, p_a, maxit=50, tol=1e-3)  {

    ID <- obsdat[, as.character(id)]
    S.nam <- strsplit(obsform, "~")[[1]][1]
    fac.terms <- strsplit(obsform, "~")[[1]][2]
    W.nam <- strsplit(fac.terms, split="\\+")[[1]][1]

    N <- length(ID)
    n <- length(unique(ID))
    clsize.vec <- as.vector(table(ID))
    m <- clsize.vec[1]
    K <- length(unique(obsdat[,S.nam])) - 1
    K_x <- length(unique(obsdat[,W.nam])) - 1

    S.c <- as.factor(obsdat[, S.nam])
    Y.levs <- levels(S.c)
    W.c <- as.factor(obsdat[,W.nam])
    X.levs <- levels(W.c)
    delta <- obsdat$delta

    ## Estimate misclassification parameters (need modification)
    vss <- obsdat[delta==1, ]
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
                glm.k <- glm(Y.mcform, family=binomial, data=vss.k)
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
                 (vss[, W.nam] == X.levs[l+1]) )
            k0.indx <- ( (vss$X == X.levs[k+1]) &
                 (vss[, W.nam] == X.levs[1]) )
            vss.k <- vss[kl.indx | k0.indx, ]
            if ( (sum(kl.indx) * sum(k0.indx)) != 0 ) {
                glm.k <- glm(X.mcform, family=binomial, data=vss.k)
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

    res <- list()
    res$gamMat <- gamMat.est
    res$varphiMat <- varphiMat.est
    res$gam.variance <- gamMat.var
    res$varphi.variance <- varphiMat.var

    return (res)
}

