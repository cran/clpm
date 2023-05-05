#' @importFrom stats model.response model.weights model.matrix terms model.frame
#' @importFrom stats .checkMFClasses .getXlevels coef delete.response
#' @importFrom stats is.empty.model lm.wfit na.pass napredict offset pt sd var vcov

#' @export
clpm <- function (formula, data, subset, na.action, weights,  
    contrasts = NULL, lambda = NULL, control = clpm.control(), 
    ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", 
        "weights", "na.action", "offset"), 
        names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf)

    if (any((w <- model.weights(mf)) < 0)) {
        stop("negative 'weights'")
    }
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    if (is.null(w)) {
        w <- rep.int(1, nrow(mf))
        alarm <- FALSE
    }
    else {
        alarm <- (w == 0)
        sel <- which(!alarm)
        mf <- mf[sel, ]
        y <- y[sel]
        w <- w[sel]
        w <- w/mean(w)
    }
    if (any(alarm)) {
        warning("observations with null weight will be dropped")
    }
    if ((n <- nrow(mf)) == 0) {
        stop("zero non-NA cases")
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(NA_real_, 
            0, ncol(y)) else numeric(), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
            0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- rep(0, n)
            z$residuals <- y
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- lpm.fit(y, x, w, lambda, control)
    }
    z$na.action <- attr(mf, "na.action")
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$model <- mf
    class(z) <- c("clpm", "lm")
    z
}

#' @export
vcov.clpm <- function(object, ...) object$covar


L_lambda <- function(beta, X, y, w, lambda, hessian = TRUE){
  
  
  
  p <- drop(X%*%cbind(beta))
  p_1 <- p - 1
  eps <- y - p
  omega0 <- as.numeric(p < 0)
  omega1 <- as.numeric(p > 1)
  
  
  
  
  Loss.i <- eps^2 + lambda*(omega1*(p_1)^2 + omega0*p^2)
  Loss <- sum(w*Loss.i)
  grad.i <- (-2*X)*(eps - lambda*(omega1*(p_1) + omega0*p))
  grad <- .colSums(w*grad.i, nrow(X), ncol(X))
  hess <- crossprod(2*w*X, X*(1 + lambda*(omega0 + omega1)))
  
  list(Loss = Loss, grad = grad, hess = hess, grad.i = grad.i, pred = p)
}

lpm.newton <- function(beta, X, y, w, lambda, tol, maxit, safeit){
  
  
  L <- L_lambda(beta, X, y, w, lambda, hessian=FALSE)
  g <- L$grad
  conv <- FALSE
  eps <- 0.1
  
  # Preliminary safe iterations, only g is used
  
  for(i in 1:safeit){
#print(paste(i, max(abs(g))))
    if(conv | max(abs(g)) < tol){break}
    cond <- FALSE # Ho migliorato la funzione? Faccio un "while" che si ferma quando il passo ? stato migliorativo
    
    while(!cond){
      new.beta <- beta - g*eps   ##### ITERAZIONE BASATA SUL GRADIENTE!!!!
      if(max(abs(g*eps)) < tol){conv <- TRUE; break} # Non ho migliorato, ma non mi muovo pi? di "tol", fammi uscire!
      L1 <- L_lambda(new.beta,X, y, w, lambda, hessian=FALSE)
      g1 <- L1$grad
      cond <- (L1$Loss < L$Loss) # Ho migliorato! Bravo. Ricalcola il gradiente e fai un'altra iterazione.
      eps <- eps*0.5
    }
    
    if(conv){break}
    g <- g1
    L <- L1
    beta <- new.beta
    eps <- eps*2
  }
  
  # Newton-Raphson
  
  alg <- "nr"
  conv <- FALSE
  eps <- 1
  L <- L_lambda(beta,X, y, w, lambda, hessian = TRUE)
  h <- L$hess
  g <- L$grad
  L<-L
  
  for(i in 1:maxit){
#print(paste(i, max(abs(g)), alg, eps,L$Loss))
    if(conv | max(abs(g)) < tol){break}
    
    H1 <- try(chol(h), silent = FALSE)
    err <- inherits(H1, "try-error")
    
    if(!err){
      if(alg == "gs"){alg <- "nr"; eps <- 1}
      delta <- chol2inv(H1)%*%g # Newton-Raphson
    }
    else{
      if(alg == "nr"){alg <- "gs"; eps <- 0.1}
      delta <- g # Gradient search
    }
    
    cond <- FALSE
    while(!cond){
      new.beta <- beta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      L1 <- L_lambda(new.beta,X, y, w, lambda, hessian = TRUE)
      h1 <- L1$hess
      g1 <- L1$grad
      
      cond <- (L1$Loss < L$Loss)
      eps <- eps*0.9
      
    }
    
    if(conv){break}
    L <- L1
    g <- g1
    h <- h1
    beta <- new.beta
    eps <- min(eps*10,1)
  }
  
  
  list(beta = beta, lambda = lambda, converged = (i < maxit), n.it = i, L = L, ok = (alg == "nr"))
}

#' @export
clpm.control <- function(tol = 1e-10, maxit, trace = FALSE){
    if(missing(maxit)){maxit <- NULL}
    if(missing(tol)) {tol <- 1e-10}
    list(tol = tol, maxit = maxit, trace = trace)
}


lpm.fit <- function(y, x, w, lambda, control){

  delta <- 1e-3
  maxit <- (if(!is.null(control$maxit)) control$maxit else ncol(x)*10)
  tol <- control$tol
  trace <- control$trace

  sX <- apply(x, 2, sd) #standard deviation per ogni variabile (intercetta,x1,x2...,xn) della model matrix
  mX <- colMeans(x)     #media di ogni variabile della model matrix
  int <- (sX == 0)      #qui indico che l'intercetta avr? sd=0, utile per dopo per vedere se c'?..
  
  
  if(any(int)){                 #any...se c'? un'intercetta
    which.int <- which(int)           #l'elemento che ? l'intercetta...(in genere posizione 1)
    mX[which.int] <- 0; sX[which.int] <- 1  #do io la media e la sd all'intercetta
  }
  else{mX <- mX*0}              #se non c'? intercetta..non le posso centrare ma posso tuttavia dividere per la deviazione standard
  XX <- scale(x, center = mX, scale = sX) #standardizzo le X, con centro media e scala la sd
  
  
  m0 <- lm.wfit(XX, y, w) #i beta sono quindi gi? standardizzati
  beta0 <- m0$coef  

  #SINGOLARITA'...beta0 potrebbe avere degli NA i quali non vogliamo passare all'algortimo sotto...
  q <- length(beta0)
  na <- (is.na(beta0))
  beta0 <- beta0[!na]
  X <- x[,!na, drop = FALSE]
  XX <- XX[,!na, drop = FALSE]
  
  #QUESTO DUE CODICI SUCCESSIVI ANDRANNO POI ANDRANNO IN OUTPUT!!!!!!!!!!!!!!!!!
  #beta <- rep(NA, q) vettore di soli NA di lunghezza beta0
  #beta[!na] <- fit$beta quelli non NA me li rimette, gli altri non stimati per connilearit? mi metter? NA in output
  #MA FORSE POTREBBERO NON SERIVIRE..MA VEDIAMO COME FAREMO L'OUTPUT
 
  # Rimuovere singolarit? (ma tenere traccia) ---> usare solo X[, !is.na(beta0), drop = FALSE]
  
  
  
  # Funzione che standardizza le X e tiene traccia di medie e sd delle x. Eseguirai i calcoli
  # sulle x standardizzate, e poi dovrai ri-convertire i coefficienti. Occhio a modelli
  # senza intercetta!! y = beta0 + beta1*(x - mu)/sigma ----> come sono i coefficienti 
  # nel modello y = gamma0 + gamma1*x?
  
  # matrice hessiana semidefinita positiva! Altrimenti error (se non trovo
  # una soluzione ben definita per nessuno dei lambda proposti dall'utente).
  if(!is.null(lambda)){
    for(j in lambda){
      fit <- lpm.newton2(beta0, XX, y, w,  tol = tol, maxit = maxit, safeit = 25, lambda = j, gamma = 1e-5)
      ok <- (fit$ok & all(fit$L$pred > -delta & fit$L$pred < 1 + delta))
      if(ok){break}
    }
  }
  else{ # restituire la soluzione "migliore", che potebbe comunque avere crossing.
    lambda <- 1
    ok <- FALSE
    while(!ok){
      fit <- lpm.newton2(beta0, XX, y, w,  tol = tol, maxit = maxit, safeit = 25, lambda = lambda, gamma = 1e-5)
      ok <- (fit$ok & all(fit$L$pred > -delta & fit$L$pred < 1 + delta))
      lambda <- lambda*1.25 # in realt? forse devo RIDURRE il lambda, se fullrank = FALSE!
      # La soluzione deve essere un compromesso fra no-crossing e full rank.
      # Serviranno delle misure di sicurezza, in lpm.newton2: a qualche valore
      # di lambda, potremmo non farcela.
    }
  }
  
  # list(X = X, y = y, w = weights, beta0 = beta0)
  
  # De-standardizza i coefficienti, calcola standard errors, corrompe un oggetto di classe "lm"
  # e restituisce il tutto.
  
  
  beta <- betaU <- rep(0, ncol(x)) # betaU = unconstrained
  beta[!na] <- fit$beta; betaU[!na] <- beta0
  if(any(int)){
    beta[!int] <- beta[!int]/sX[!int]
    beta[int] <- beta[int] - sum((beta*mX)[!int])

    betaU[!int] <- betaU[!int]/sX[!int]
    betaU[int] <- betaU[int] - sum((betaU*mX)[!int])
  }
  else{beta <- beta/sX; betaU <- betaU/sX}
  beta[na] <- betaU[na] <- NA
  names(beta) <- names(betaU) <- colnames(x)

  # Matrice di covarianza asintotica
u <- list(beta = beta[!na], X = X, y = y, w = w)
  L <- L_lambda(betaU[!na], X, y, w, lambda = 0, hessian = TRUE)
a <- L
A <- fit
W <- w
BetaU <- betaU
  Omega <- crossprod(L$grad.i, L$grad.i*w)
o <- Omega
  V <- chol2inv(chol(L$hess%*%chol2inv(chol(Omega))%*%L$hess))
#V <- chol2inv(chol(L$hess))
  covar <- matrix(NA, q, q)
  covar[!na, !na] <- V
  rownames(covar) <- colnames(covar) <- colnames(x)
  
  z <- list(coefficients = beta, covar = covar, 
        residuals = y - fit$L$pred, rank = m0$rank, 
        fitted.values = fit$L$pred, assign = attr(x, "assign"), qr = m0$qr, 
        df.residual = length(y) - m0$rank, obj.function = L$Loss, gradient = fit$L$grad, 
        convergence = fit$converged, 
        n.it = fit$n.it, control = control, lambda = fit$lambda)
}

lpm.newton2 <- function(beta, XX, y, w, lambda, tol, maxit, safeit, gamma = 1e-5){
  n <- length(y)
  safeit <- 20
  fit <- lpm.newton(beta, XX, y, w, lambda = lambda, tol = tol, maxit = maxit, safeit = safeit)
  g <- max(abs(fit$L$grad))/n # (g=gradiente medio)....diviso per n perch? il gradiente ? una somma di vari contributi al gradiente quindi sar? sempre pi? grande se aumento le osservazioni
  while(g > gamma){ # SISTEMA DI SICUREZZA!!
    safeit <- safeit + 10
    fit <- lpm.newton(fit$beta, XX, y, w, lambda=lambda, tol = tol, maxit = maxit, safeit = safeit)
    g <- max(abs(fit$L$grad))/n
    gamma<- gamma*2
  }
  fit
}

#' @export
summary.clpm <- function (object, correlation = FALSE, symbolic.cor = FALSE, 
                         ...) 
{
  #Q <- getFromNamespace("qr.lm", "stats")
  Q <- function (x, ...) {x$qr}
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {  
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, 
                                                                 c("Estimate", "Std. Error", "t value", 
                                                                   "Pr(>|t|)")))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
    if (correlation) 
      ans$correlation <- ans$cov.unscaled
    return(ans)
  }
  
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- Q(object)
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 
      1e-30) 
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  #    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  #    se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  se <- sqrt(diag(object$covar))[Qr$pivot[p1]]
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
                            `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), 
                                                                  rdf, lower.tail = FALSE))
  ans$aliased <- is.na(z$coefficients)
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
                        numdf = p - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  #    ans$cov.unscaled <- R
  ans$covar <- object$covar[Qr$pivot[p1], Qr$pivot[p1]]
  #    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
  #        1)]
  dimnames(ans$covar) <- dimnames(ans$coefficients)[c(1,1)]
  if (correlation) {
    #        ans$correlation <- (R * resvar)/outer(se, se)
    ans$correlation <- ans$covar/outer(se, se)
    #        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    dimnames(ans$correlation) <- dimnames(ans$covar)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  class(ans) <- c("summary.lm", "summary.lpm")
  ans
}

#' @export
predict.clpm <- function (object, newdata, se.fit = FALSE, ...) {
  #Q <- getFromNamespace("qr.lm", "stats")
  Q <- function (x, ...) {x$qr}
  tt <- terms(object)
  if (!inherits(object, "lm")) 
    warning("calling predict.lm(<fake-lm-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  n <- length(object$residuals)
  p <- object$rank
  p1 <- seq_len(p)
  piv <- if (p) 
    Q(object)$pivot[p1]
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
    warning("prediction from a rank-deficient fit may be misleading")
  beta <- object$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  if (!is.null(offset)) 
    predictor <- predictor + offset
  
  if (se.fit) {
    
    V <- vcov(object)
    V <- V[piv, piv, drop = FALSE]
    x <- X[,piv, drop = FALSE]
    se <- sqrt(diag(x%*%V%*%t(x)))
  }
  
  if (missing(newdata) && !is.null(na.act <- object$na.action)) {
    predictor <- napredict(na.act, predictor)
    if (se.fit) 
      se <- napredict(na.act, se)
  }
  if (se.fit) 
    list(fit = predictor, se.fit = se)
  else predictor
}



