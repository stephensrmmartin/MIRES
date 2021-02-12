##' Generates unidimensional data for testing the HMRE/MIRES approach.
##'
##' \code{mipattern} is a list specifying a pattern of MI.
##' The first entry should be a string specifying one of: constant, random, none, items, params, or custom.
##' The other entries depend on the specification desired, as described below.
##' \describe{
##'   \item{constant}{ (2) Numeric: All RE SDs are set to this value.}
##'   \item{random}{ (2) Numeric: All RE SDs are generated between 0 and this value.}
##'   \item{none}{ All RE SDs set to zero (full invariance).}
##'   \item{items}{ (2) Numeric: Value of RE SDs for (3) items. E.g., "items", .4, 4 would set all parameters for items 1 to 4 to have an RE SD of .4.}
##'   \item{params}{ (2) Numeric: Value of RE SD for parameter type (3), where (3) is an integer (0: Loadings, 1: Residual SDs, 2: Intercepts). E.g., "params", .4, 2 would set all RE-SDs of intercepts to be .4.}
##'   \item{custom}{ (2) Numeric: Specify all 3J RE-SDs manually in order of loadings, residual log SDs, and intercepts.}
##' }
##' Note that this is \emph{not} the generative model specified by MIRES, but a convenience function for meeting the bare assumptions while generating MI or non-MI data.
##' 
##' @title Unidimensional data generation.
##' @param J Integer. Number of indicators.
##' @param K Integer. Number of groups.
##' @param n Integer. Number of observations within group.
##' @param fixed Named List. lambda, resid_log, nu, in that order.
##' @param mipattern List. See details.
##' @param etadist (Default: NULL). NULL, "std", or a list of two. If NULL (Default), all groups have latent scores distributed standard normal. If "std", means are standard normal, and (log) latent SDs are also standard normal (i.e., standard log normal; product to 1). If a list, slot one provides the K means, slot two provides the K log SDs. These should have a mean of zero (sum to zero) and a mean of 1 (product to 1), respectively.
##' @return List of meta(data), params, data, and a data frame.
##' @author Stephen Martin
##' @importFrom mvtnorm rmvnorm
##' @keywords internal
datagen_uni <- function(J, K, n, fixed, mipattern, etadist = NULL) {
    N <- K * n # Total rows. Assumes equal number of observations per group.

    group <- rep(1:K, each = n)

    etadist <- etadist %IfNull% "zeroes"
    eta_mean <- rep(0, K)
    eta_sd <- rep(1, K)
    if(etadist == "std") {
        eta_mean <- rnorm(K)
        eta_sd <- exp(rnorm(K))
    }
    if(is.list(etadist)) {
        eta_mean <- etadist[[1]]
        eta_sd <- etadist[[2]]
    }

    paramSD <- switch(mipattern[[1]],
                      constant = rep(mipattern[[2]], J * 3),
                      random = runif(J * 3, 0, mipattern[[2]]),
                      none = rep(0, J * 3),
                      items = rep(c(
                          rep(mipattern[[2]], mipattern[[3]]),
                          rep(0, J - mipattern[[3]])
                      ), 3),
                      params = {
                          p <- rep(0, 3*J)
                          p[(1 + mipattern[[3]] * J) : (mipattern[[3]] * J + J)] <- rep(mipattern[[2]], J)
                          p
                      },
                      custom = mipattern[[2]]
                      )

    # Generate REs [Assumption: No correlation]
    # DONE : Change REs so that the minimum RE lambdas for each j + fixed[j] >= 0; as per HMRE assumption.
    # Maybe just take samples from abs(rmvnorm(K, fixefs, sigma)) to enforce the positiveness.
    # Then subtract off the REs to get the randoms.
    ## random <- rmvnorm(K, sigma = diag(paramSD^2))
    random_coefs <- rmvnorm(K, unlist(fixed), diag(paramSD^2))
    random_coefs[, 1:J] <- abs(random_coefs[, 1:J])
    random <- t(t(random_coefs) - unlist(fixed))

    lambda <- t(fixed$lambda)
    resid_log <- fixed$resid
    nu <- fixed$nu


    # Generate observations [Assumption: All means are zero. TODO: Fix this]
    eta <- rnorm(N, 0, 1)*eta_sd[group] + eta_mean[group]

    # TODO : Define nu, lambda, resid_log; figure out interface for that.
    y <- rep(1, N) %*% t(nu) + eta %*% lambda + # Fixed effects
        random[group, (2*J + 1) : (3 * J)] + # Random intercepts
        eta * random[group, 1 : J] # Random loadings

    err <- matrix(rnorm(N*J, 0, exp(rep(1, N) %*% t(resid_log) + # Fixed resid
                               random[group, (J + 1) : (2 * J)] # Random resid
                               )),
                  N, J) # Matrixify
    y <- y + err
    colnames(y) <- paste0("x_", 1:J)

    data <- list(N = N, J = J, n = n, K = K, x = y, group = group)
    params <- list(lambda = fixed$lambda,
                   resid = fixed$resid,
                   nu = fixed$nu,
                   random = random,
                   eta = eta,
                   eta_mean = eta_mean,
                   eta_sd = eta_sd,
                   lambda_random = random[, 1 : J],
                   resid_random = random[, (J + 1) : (2 * J)],
                   nu_random = random[, (2 * J + 1) : (3 * J)],
                   random_sigma = paramSD,
                   RE_cor = diag(1, 3*J, 3*J) # No RE cor currently
                   ) 
    meta <- list(N = N, J = J, n = n, K = K, F = 1, mipattern = mipattern)
    df <- cbind(as.data.frame(y), group = group)

    list(data = data,
         params = params,
         meta = meta,
         df = df)
}

##' @title Paper simulation function (For historical purposes)
##' @param J Integer. Number of indicators.
##' @param K Integer. Number of groups.
##' @param n Integer. Observations per group.
##' @param paramSDPattern List.
##' @return List.
##' @author Stephen Martin
##' @keywords internal
generateData <- function(J,K,n,paramSDPattern){
  # Set hyperdata
  J <- J # Number of items
  K <- K # Number of groups
  N <- K*n # Total N
  group <- rep(1:K,each=n)
  
  # Generate fixed nu, (positive) loading, (log) residual variance
  nu <- rnorm(J,0,.2)
  lambda <- runif(J,.4,.98)
  resid <- sqrt(1 - lambda^2)
  resid_log <- log(resid)
  
  # Generate random nu, loading, (log) residual variance
  paramSD <- switch(paramSDPattern[[1]],
     constant = rep(paramSDPattern[[2]],J*3),
     random = runif(J*3,0,paramSDPattern[[2]]),
     none = rep(0,J*3),
     items = rep(c(rep(paramSDPattern[[2]],paramSDPattern[[3]]),rep(0,J-paramSDPattern[[3]])),3),
     params = {p <- rep(0,3*J); p[(1 + paramSDPattern[[3]]*J):(paramSDPattern[[3]]*J + J)] <- rep(paramSDPattern[[2]],J); p}
  )
  # Generate REs
  paramRE <- rmvnorm(K,sigma=diag(paramSD^2))
  
  # Generate observations
  theta <- rnorm(N,0,1) # Standardized
  ds <- data.frame(subject = 1:N,group=group,theta=theta)
  genResponse <- function(group,theta,nu,lambda,resid_log,paramRE){
    J <- length(lambda)
    lambda_k <- lambda + paramRE[group,1:J]
    resid_log_k <- resid_log + paramRE[group,(J+1):(2*J)]
    nu_k <- nu + paramRE[group,(2*J + 1):(3*J)]
    y <- nu_k + theta*lambda_k + rnorm(J,0,exp(resid_log_k))
    t(y)
  }
  ds[,paste0('x',1:J)] <- t(mapply(genResponse,ds$group,ds$theta, MoreArgs=list(nu=nu,lambda=lambda,resid_log=resid_log,paramRE=paramRE)))
  list(K=K,J=J,n=n,N=N,paramRE = paramRE,paramSD = paramSD,lambda=lambda,resid_log=resid_log,nu=nu,theta=theta,ds=ds)
  
}

generateDataMeans <- function(J,K,n,paramSDPattern,thetaMeans,thetaSDs){
  # Set hyperdata
  J <- J # Number of items
  K <- K # Number of groups
  N <- K*n # Total N
  group <- rep(1:K,each=n)
  
  # Generate fixed nu, (positive) loading, (log) residual variance
  nu <- rnorm(J,0,.2)
  lambda <- runif(J,.4,.98)
  resid <- sqrt(1 - lambda^2)
  resid_log <- log(resid)
  
  # Generate random nu, loading, (log) residual variance
  paramSD <- switch(paramSDPattern[[1]],
     constant = rep(paramSDPattern[[2]],J*3),
     random = runif(J*3,0,paramSDPattern[[2]]),
     none = rep(0,J*3),
     items = rep(c(rep(paramSDPattern[[2]],paramSDPattern[[3]]),rep(0,J-paramSDPattern[[3]])),3),
     params = {p <- rep(0,3*J); p[(1 + paramSDPattern[[3]]*J):(paramSDPattern[[3]]*J + J)] <- rep(paramSDPattern[[2]],J); p}
  )
  # Generate REs
  paramRE <- rmvnorm(K,sigma=diag(paramSD^2))
  
  # Generate observations
  theta <- rnorm(N,thetaMeans[group],thetaSDs[group]) # Standardized
  ds <- data.frame(subject = 1:N,group=group,theta=theta)
  genResponse <- function(group,theta,nu,lambda,resid_log,paramRE){
    J <- length(lambda)
    lambda_k <- lambda + paramRE[group,1:J]
    resid_log_k <- resid_log + paramRE[group,(J+1):(2*J)]
    nu_k <- nu + paramRE[group,(2*J + 1):(3*J)]
    y <- nu_k + theta*lambda_k + rnorm(J,0,exp(resid_log_k))
    t(y)
  }
  ds[,paste0('x',1:J)] <- t(mapply(genResponse,ds$group,ds$theta, MoreArgs=list(nu=nu,lambda=lambda,resid_log=resid_log,paramRE=paramRE)))
  list(K=K,J=J,n=n,N=N,paramRE = paramRE,paramSD = paramSD,lambda=lambda,resid_log=resid_log,nu=nu,theta=theta,ds=ds)
  
}
