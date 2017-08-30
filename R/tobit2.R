
inf.omit = function(t) t[is.finite(t)]
#logit = function(t) log(t/(1-t))
#logitinv = function(t) exp(t)/(exp(t)+1)
#logitmueta = function(t) exp(t)/(1+exp(t))^2
tobit2 <- function(link=list("identity","identity","log","logit2" ), censoring = FALSE, rho=NULL, eps = 1e-3) {
## Extended family for Gaussian location scale model...
## so mu is mu1 and tau=1/sig is mu2
## tau = 1/(b + exp(eta)) eta = log(1/tau - b)
## 1. get derivatives wrt mu, tau
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
## the first derivatives of the log likelihood w.r.t
## the first and second parameters...
  env <- new.env(parent = .GlobalEnv)
     assign(".censoring", censoring, envir = env)
     assign("fl1", fl1, envir = env)
     assign("fl2", fl2, envir = env)
     assign("fl3", fl3, envir = env)
     assign("fl4", fl4, envir = env)
     assign("inf.omit", inf.omit, envir = env)

  ## first deal with links and their derivatives...
  if (length(link)<3+is.null(rho)) stop("tobitlss requires 4 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity","sqrt"), c("inverse", "log", "identity","sqrt"), "log", "logit2")
  stats <- list()

for (i in 1:3){
  if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
  stop(link[[i]]," link not available for mu parameter of gaulss")
  fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
  fam <- fix.family.link(fam)
  stats[[i]]$d2link <- fam$d2link
  stats[[i]]$d3link <- fam$d3link
  stats[[i]]$d4link <- fam$d4link
}
logit = function(t) log(t/(1-t))
logitinv = function(t) exp(t)/(exp(t)+1)
logitmueta = function(t) exp(t)/(1+exp(t))^2
if (is.null(rho)){
  if (link[[4]] %in% okLinks[[4]]) { ## creating the logb link
    stats[[4]] <- list()
    stats[[4]]$valideta <- function(eta) TRUE 
    stats[[4]]$link = link[[4]]
    stats[[4]]$linkfun <- function(mu){ t = ( (mu+1)/2);  log(t/(1-t)) } #.Call("C_logit_link", (mu+1)/2) #eval(parse(text=paste("function(mu) log(1/mu -",b,")")))
    stats[[4]]$linkinv <- function(eta) 2*exp(eta)/(exp(eta)+1) -1 #.Call("C_logit_linkinv", eta) -1 #eval(parse(text=paste("function(eta) 1/(exp(eta) +",b,")")))
    stats[[4]]$mu.eta <- function(eta) 2*exp(eta)/(1+exp(eta))^2  #.Call("C_logit_mu_eta", eta)

    environment(stats[[4]]$linkfun) <- environment(stats[[4]]$linkinv) <- environment(stats[[4]]$mu.eta) <- asNamespace("stats")

    stats[[4]]$d2link <-function(mu) 0.25*(1/(1 - (mu/2+0.5))^2 - 1/(mu/2+0.5)^2 )
    stats[[4]]$d3link <- function (mu) 0.125*( 2/(1 - (mu/2+0.5))^3 + 2/(mu/2+0.5)^3)
    stats[[4]]$d4link <- function (mu) 0.0625*(6/(1 - (mu/2+0.5))^4 - 6/(mu/2+0.5)^4)
  } else stop(link[[4]]," link not available for  parameter of gaulss")
}

  residuals <- function(object,type=c("deviance","pearson","response")) {
      type <- match.arg(type)
      rsd <- object$y-object$fitted[,1]
      if (type=="response") return(rsd) else
      return( (rsd/object$fitted[,3])) ## (y-mu)/sigma 
    }
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## in principle the following seems reasonable, but because no
    ## price is paid for the high null variance, it leads to silly
    ## % deviance explained...
    #er <- fitNull(G$y,G$family,G$w,G$offset,nlp=length(attr(G$X,"lpi")),tol=1e-7)
    #object$null.deviance <- sum(((object$y-er$mu[,1])*er$mu[,2])^2*G$w)
inf.omit = function(t) t[is.finite(t)]

    object$null.deviance <- sum(inf.omit(((object$y-mean(object$y))/object$fitted[,2])^2))
  })

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
  ## function defining the gamlss Gaussian model log lik. 
  ## N(mu,sigma^2) parameterized in terms of mu and log(sigma)
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.

    if (!is.null(offset)) offset[[5]] <- 0

   jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    if (!is.null(offset[[1]])) eta <- eta + offset[[1]]
    mu1 <- family$linfo[[1]]$linkinv(eta)
#mu1 is uncensored

    eta2 <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]]
    if (!is.null(offset[[2]])) eta2 <- eta2 + offset[[2]]
    mu2 <-  family$linfo[[2]]$linkinv(eta2)
#mu2 is for censoring

    eta3 <- X[,jj[[3]],drop=FALSE]%*%coef[jj[[3]]]
    if (!is.null(offset[[3]])) eta3 <- eta3 + offset[[3]]
    sigma <-  family$linfo[[3]]$linkinv(eta3)

if (is.null(family$rho)){  
  eta4 <- X[,jj[[4]],drop=FALSE]%*%coef[jj[[4]]]
    if (!is.null(offset[[4]])) eta4 <- eta4 + offset[[4]]
    rho <-  family$linfo[[4]]$linkinv(eta4)
estRho = TRUE
} else {
rho = family$rho
estRho = FALSE
}

#lower = .lower
#upper = .upper
censoring = family$censoring
eps = family$eps
l =  (2*log(sigma)+2*((y-mu1)^2)/(2*sigma^2)   -2*   log(1- pnorm((-mu2 - rho * (y-mu1)/sigma )/sqrt(1-rho^2 + eps)) ))
#l =  2*log(sigma)  -2*log(dnorm((y - mu1)/sigma)) -2*   log(1-pnorm((mu2 + rho * (y-mu1)/sigma - upper)/sqrt(1-rho^2)) - pnorm((-mu2 - rho * (y-mu1)/sigma + lower)/sqrt(1-rho^2)) )
l[censoring] =  -2*log(pnorm(-mu2))[censoring]
l = -sum(l)
    n <- length(y)
    l1 <- matrix(0,n,3+estRho)
#mu1 mu2 sigma rho
    if (deriv>0) {

tmp = fl1(y, mu1, mu2, sigma, rho, censoring, eps)
tmp = -tmp
if (!estRho) tmp = tmp[,-4]
l1 = tmp
tmp = matrix(fl2(y,mu1,mu2,sigma,rho, censoring, eps), n)[, c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 11L, 12L, 16L)]
tmp = -tmp
#remove derivatives that are to do with rho
if (!estRho) tmp = tmp[,-c(4,7,9,10)]
l2 = tmp

      ## need some link derivatives for derivative transform
if (!estRho){
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta2), family$linfo[[3]]$mu.eta(eta3))
      g2 <- cbind(family$linfo[[1]]$d2link(mu1),family$linfo[[2]]$d2link(mu2),family$linfo[[3]]$d2link(sigma) )
} else {

      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta2), family$linfo[[3]]$mu.eta(eta3), family$linfo[[4]]$mu.eta(eta4)  )
      g2 <- cbind(family$linfo[[1]]$d2link(mu1),family$linfo[[2]]$d2link(mu2),family$linfo[[3]]$d2link(sigma),family$linfo[[4]]$d2link(rho) )
}
    }
 
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
tmp = matrix(fl3(y,mu1,mu2,sigma,rho, censoring, eps), n)[,c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 11L, 12L, 16L, 22L, 23L, 24L, 27L, 
28L, 32L, 43L, 44L, 48L, 64L) ]
tmp = -tmp
if (!estRho) tmp = tmp[,-c(4,7,9,10,13,15,16,18,19,20)]
l3 = tmp

if (estRho){
      g3 <- cbind(family$linfo[[1]]$d3link(mu1),family$linfo[[2]]$d3link(mu2),family$linfo[[3]]$d3link(sigma),family$linfo[[4]]$d3link(rho) )
} else {
      g3 <- cbind(family$linfo[[1]]$d3link(mu1),family$linfo[[2]]$d3link(mu2),family$linfo[[3]]$d3link(sigma))
}
    }

    if (deriv>3) {

tmp = matrix(fl4(y, mu1,mu2,sigma,rho, censoring, eps), n)[, c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 11L, 12L, 16L, 22L, 23L, 24L, 27L, 
28L, 32L, 43L, 44L, 48L, 64L, 86L, 87L, 88L, 91L, 92L, 96L, 107L, 
108L, 112L, 128L, 171L, 172L, 176L, 192L, 256L)]
tmp = -tmp
if (!estRho) tmp = tmp[,-c(4,7,9,10,13,15,16,18,19,20,23,25,26,28,29,30,32,33,34,35)]
l4 = tmp

if (estRho){
      g4 <- cbind(family$linfo[[1]]$d4link(mu1),family$linfo[[2]]$d4link(mu2),family$linfo[[3]]$d4link(sigma),family$linfo[[4]]$d4link(rho) )
} else {
      g4 <- cbind(family$linfo[[1]]$d4link(mu1),family$linfo[[2]]$d4link(mu2),family$linfo[[3]]$d4link(sigma))
}
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- l; ret



  } ## end ll gaulss

  initialize <- expression({
  ## idea is to regress g(y) on model matrix for mean, and then 
  ## to regress the corresponding log absolute residuals on 
  ## the model matrix for log(sigma) - may be called in both
  ## gam.fit5 and initial.spg... note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.

        n <- rep(1, nobs)
  mustart <- drop(y)
 y = drop(y)
   use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
	if (!is.null(offset)) offset[[5]] <- 0
        jj <- attr(x,"lpi")

        start <- rep(0,ncol(x))
#	start[jj[[3]][1]] = max(0, log(sd(y[!censoring])))
	#experimentally fill in sigma?

#browser()
      }

#      n <- rep(1, nobs)
#      ## should E be used unscaled or not?..
#      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
#      if (is.null(start)) {
#        jj <- attr(x,"lpi")
#	#offs <- attr(x,"offset")
#	if (!is.null(offset)) offset[[3]] <- 0
#        start <- rep(0,ncol(x))
#        yt1 <- if (family$link[[1]]=="identity") y else 
#               family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
#	if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
#        x1 <- x[,jj[[1]],drop=FALSE]
#        e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
#        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
#        if (use.unscaled) {
#          qrx <- qr(rbind(x1,e1))
#          x1 <- rbind(x1,e1)
#          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
#          startji[!is.finite(startji)] <- 0       
#        } else startji <- pen.reg(x1,e1,yt1)
#        start[jj[[1]]] <- startji
#        lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
#	if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
#        x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
#        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
#        if (use.unscaled) {
#          x1 <- rbind(x1,e1)
#          startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
#          startji[!is.finite(startji)] <- 0
#        } else startji <- pen.reg(x1,e1,lres1)
#        start[jj[[2]]] <- startji
#      }
#
#


  }) ## initialize gaulss
     environment(ll) <-environment(residuals)<- environment(postproc) <- environment(initialize) <- env

  structure(list(family="tobitlss",ll=ll,link=paste(link),nlp=3+is.null(rho),
    tri = trind.generator(3+is.null(rho)), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2, ## can use full Newton here
    censoring=censoring,  rho = rho, eps = eps),class = c("general.family","extended.family","family"))
} ## end gaulss


#source("~/tobit4derivatives.R")
