

## ===================== Sunfunctions =========================== ##
expit <- function(x) {
	return(exp(x)/(1+exp(x)))
}

logit <- function(x) {
	return(log(x/(1-x)))
}


## Get a design matrix, given the formula and data set
getDM <- function(formula, data)
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    #m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$cOR.str <- m$Mv <- m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
    #if (is.null(m$id))
    #    m$id <- as.name("id")
    #if (!is.null(m$na.action) && m$na.action != "na.omit") {
    #    warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
    #    m$na.action <- as.name("na.omit")
    #}
    m$na.action <- as.name("na.pass")
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    y <- as.matrix(model.extract(m, "response"))
    # Yuliang: contrasts.arg = NULL?
    x <- model.matrix(Terms, m, contrasts.arg = NULL)

    return(DM=x)
}


## Get the response vector, given the formula and data set
getResp <- function(formula, data)
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    #m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$cOR.str <- m$Mv <- m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
    #if (is.null(m$id))
    #    m$id <- as.name("id")
    if (!is.null(m$na.action) && m$na.action != "na.omit") {
        warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
        m$na.action <- as.name("na.omit")
    }
    m$na.action <- as.name("na.pass")
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    y <- as.matrix(model.extract(m, "response"))
    y <- as.vector(y)
    return(Resp=y)
}


### ===================================================================================
### 		Correlated binary responses
### ===================================================================================

## ---------------------------------------------------
##    outer() with a vector argument
## ---------------------------------------------------
outerX <- function(x, FUN="*") {
	xx <- outer(x,x, FUN)
 	return(xx)
}


## ----------------------------------------------------------
##  Function to get the lower triangle elements
##    (not including diagonal) in a matrix.
## ----------------------------------------------------------
LowerTri <- function(x, vec=FALSE, FUN="*", ...) { ##lower.tri=TRUE, diag=FALSE,
	if (length(dim(x))==2) {
		n <- nrow(x)
		if (n==ncol(x)) {
			lowerTri <- NULL
			for(j in 1:(n-1)) lowerTri <- c(lowerTri, x[(j+1):n, j])
			return(lowerTri)
		} else stop("The matrix is not a square matrix")
	} else if(vec) {
		FUN <- match.fun(FUN)
		xx <- outer(x,x, FUN, ...)
		n <- nrow(xx)
		lowerTri <- NULL
		for(j in 1:(n-1)) lowerTri <- c(lowerTri, xx[(j+1):n, j])
		return(lowerTri)
	}
}

UpperTri <- function(x, vec=FALSE, FUN="*", ...) {
	if (length(dim(x))==2) {
		n <- nrow(x)
		if (n==ncol(x)) {
			upperTri <- NULL
			for(j in 1:(n-1)) upperTri <- c(upperTri, x[j, (j+1):n])
			return(upperTri)
		} else stop("The matrix is not a square matrix")
	} else if(vec) {
		FUN <- match.fun(FUN)
		xx <- outer(x,x, FUN, ...)
		n <- nrow(xx)
		upperTri <- NULL
		for(j in 1:(n-1)) upperTri <- c(upperTri, xx[j, (j+1):n])
		return(upperTri)
	}
}

## ----------------------------------------------------------
##	Function to get correlation matrix from $\rho$
## ----------------------------------------------------------
rho2CorrM <- function(x, lower.tri=TRUE, diag=FALSE)
{
	m <- round(1/2+sqrt(2*length(x)+1/4))
	if (length(x)==(m*(m-1)/2)) {
		CorrM <- diag(nrow=m)
		for(j in 1:(m-1)) {
		CorrM[(j+1):m, j] <- x[(m*(j-1)-j*(j-1)/2+1): (m*j-j*(j+1)/2)]
		CorrM[j, (j+1):m] <- x[(m*(j-1)-j*(j-1)/2+1): (m*j-j*(j+1)/2)]
		}
		return(CorrM)
	} else stop("Length of input vector is not correct")
}


## ----------------------------------------------------------
##	Enumerate all $2^m$ possible binary strings
## ----------------------------------------------------------
binstr.ENUM <- function(m)
{
	x <- matrix(nrow=m, ncol=2^m)
	for (j in 1:m) {
		x[j,] <- as.vector(outer(c(rep(0,2^(j-1)),rep(1,2^(j-1))),
			rep(1,2^(m-j))))
	}
	return(x)
}


## --------------------------------------------
##    Function to summarize results
## --------------------------------------------
getpropmea <- function(pars, est.mat, var.mat=NULL)
{
   bias <- apply(est.mat, 2, mean) - pars
   empvar <- apply(est.mat, 2, var)
   rb <- bias/pars*100
   if(!is.null(var.mat)) {
      amv <- apply(var.mat, 2, mean)
      bv <-  amv - empvar
      rbv <- bv/empvar*100
      nrep <- nrow(est.mat)
      cp <- apply(abs((est.mat-outer(rep(1,nrep),pars)))<
         sqrt(var.mat)*1.96, 2 , mean)
   } else {
      amv <- NULL
      bv <-  NULL
      rbv <- NULL
      nrep <- nrow(est.mat)
      cp <- NULL
   }
   prop <- data.frame(
      RB=rb,
      EV=empvar,
      AMV=amv,
      RBV=rbv,
      CP=cp
   )
   return(prop)
}


sum.table <- function(pars, estMAT, varMAT=NULL)
{
   flag.var <- !is.null(varMAT)
   bad.est.indx <- apply(is.na(estMAT)|(abs(estMAT)>5), 1, FUN="any")
   bad.indx <- bad.est.indx
   if(flag.var) {
      bad.var.indx <- apply(is.na(varMAT)|(varMAT<=0)|(varMAT>10),
         1, FUN="any")
      bad.indx <- bad.est.indx|bad.var.indx
      var.mat <- varMAT[!bad.indx,]
   }
   est.mat <- estMAT[!bad.indx,]

   bias <- apply(est.mat, 2, mean) - pars
   empvar <- apply(est.mat, 2, var)
   rb <- bias/pars*100
   if(flag.var) {
      amv <- apply(var.mat, 2, mean)
      bv <-  amv - empvar
      rbv <- bv/empvar*100
      nrep <- nrow(est.mat)
      cp <- apply(abs((est.mat-outer(rep(1,nrep),pars)))<
         sqrt(var.mat)*1.96, 2 , mean)
   } else {
      amv <- NULL
      bv <-  NULL
      rbv <- NULL
      nrep <- nrow(est.mat)
      cp <- NULL
   }
   prop <- data.frame(
      RB=rb,
      EV=empvar,
      AMV=amv,
      #RBV=rbv,
      CP=cp
   )
   cat(paste(nrep, "good simulation replicates"), "\n")
   return(prop)
}


get.sumtab <- function(pars, res.out, digit=4)
{
   p_t <- length(pars)
   conv.indx <- ( res.out[, ncol(file.out)-1]==1 &
       #apply(abs(res.out[,1:p_t]) < 3, 1, prod) &
       apply(res.out[,(p_t+1):(2*p_t)] < 5, 1, prod) &
       apply(res.out[,(p_t+1):(2*p_t)] > 0, 1, prod) )
   nconv <- sum(conv.indx)
   est.out <- res.out[conv.indx, 1:p_t]
   var.out <- res.out[conv.indx, (p_t+1):(p_t+p_t)]

   bias <- apply(est.out, 2, mean) - pars
   RB <- bias/pars*100
   EV <- apply(est.out, 2, var)
   AMV <- apply(var.out, 2, mean)
   biasV <-  AMV - EV
   RBV <- biasV/EV*100
   CP <- apply(abs((est.out-outer(rep(1,nconv),pars)))<
         sqrt(var.out)*1.96, 2 , mean)
   sumtab <- data.frame(
       #bias, biasV,
       RB, EV, AMV,
       #RBV,
       CP)
   return(round(sumtab, digit))
}



