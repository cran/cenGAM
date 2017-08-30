tobit1 <- function(link="identity", left.threshold= -Inf, right.threshold=Inf, theta = NULL,  initial.theta = 0){

  #left and right are censorship functions, possibly one per observation

    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    okLinks <- c("identity")
    if (linktemp %in% okLinks)
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for gaussian family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }

 Theta <-  NULL; n.theta <- 1
    if (!is.null(theta)&&theta!=0) {
       if (theta>0) {
           iniTheta <- Theta <- log(theta) ## fixed theta supplied
           n.theta <- 0 ## signal that there are no theta parameters to estimate
       } else iniTheta <- log(-theta) ## initial theta supplied
    } else iniTheta <- initial.theta ##  inital log theta value

  env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    assign(".left.threshold", left.threshold, envir = env)
    assign(".right.threshold", right.threshold, envir = env)
    getright.threshold <- function() get(".right.threshold") 
    getleft.threshold <- function() get(".left.threshold") 

    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") 
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) {
        th <- get(".Theta")
        rep( exp(2*th) , length(mu))
    } # probably not very important
  validmu <- function(mu) all(is.finite(mu))


    dev.resids <- function(y, mu, wt,theta=NULL) {
    ## '-'2*loglik instead of deviance in REML/ML expression
#	d1 = (wt ==1)#y[,1] <= y[,2]
#	d2 = (wt== 3)#y[,1] >= y[,3]
right.threshold = get(".right.threshold")
left.threshold = get(".left.threshold")
k = (y>left.threshold) + (y>= right.threshold)+1
      if (is.null(theta)) theta <- get(".Theta")
      lth = theta
      theta <- exp(theta) ## note log theta supplied
res=rep(NA, length(y))
res[k==2] = 2*(lth- log(dnorm((y - mu)[k==2] /theta))) 
res[k==1] = - 2* log(pnorm( (y - mu)[k==1] /theta))
res[k==3] = - 2* log(pnorm( ( mu - y )[k==3] /theta))
res *wt
#res=	2*( (!d1) * (!d2)*(lth- log(dnorm((y - mu)/theta))) - (d1 * log(pnorm( (y - mu)/theta))) - (d2 * log(pnorm( ( mu - y )/theta))))  #???

    }

Dd = function(y, mu, theta, wt, level=0) {
    ## derivatives of the -2*loglik...
right.threshold = get(".right.threshold")
left.threshold = get(".left.threshold")
k = (y>left.threshold) + (y>= right.threshold)+1

	d1 = (k==1)
	d2 = (k==3)
y1 = y[d1]
y2 = y[d2]
mu1 = mu[d1]
mu2 = mu[d2]
       ltheta <- theta
      theta <- exp(theta)
	del = function(x) dnorm(x)/pnorm(x)
	erfc = function(x) 2*pnorm(-x*sqrt(2))
r = list()
r$Dmu <- 2*( (!d1) * (!d2)* (mu-y )/theta^2 )
r$Dmu[d1] <- 2*((del((y1 - mu1)/theta)/theta))
r$Dmu[d2] <- 2*( - (  del ( ( mu2 - y2 )/theta)/theta))
r$Dmu=wt*r$Dmu
r$Dmu2 <-2*( (!d1) * (!d2)/theta^2   )
r$Dmu2[d1] <-2*(- ((sqrt(2/pi)*(mu1-y1)*exp(-(mu1-y1)^2/(2*theta^2)))/(theta^3*erfc((mu1-y1)/(sqrt(2)*theta)))-(2*exp(-(mu1-y1)^2/theta^2))/(pi*theta^2*erfc((mu1-y1)/(sqrt(2)*theta))^2)))
r$Dmu2[d2] <-2*(- ((sqrt(2/pi)*(y2-mu2)*exp(-(y2-mu2)^2/(2*theta^2)))/(theta^3*erfc((y2-mu2)/(sqrt(2)*theta)))-(2*exp(-(y2-mu2)^2/theta^2))/(pi*theta^2*erfc((y2-mu2)/(sqrt(2)*theta))^2) ))
r$Dmu2=wt*r$Dmu2
r$EDmu2 <- r$Dmu2
      if (level>0) { ## quantities needed for first derivatives

r$Dth =  2*( -(!d1) * (!d2)* ((y-mu)^2 - theta^2)/theta^2 )
r$Dth[d1] =  2*(- (sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))))
r$Dth[d2] =  2*(- (sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))))
r$Dth=wt*r$Dth
r$Dmuth = 2*( (!d1) * (!d2)*  (2*(y-mu) )/theta^2 )
r$Dmuth[d1] = 2*(((sqrt(2/pi)*(mu1-y1)^2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))-(2*(mu1-y1)*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))))
r$Dmuth[d2] = 2*(-((sqrt(2/pi)*(y2-mu2)^2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))-(2*(y2-mu2)*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))))
r$Dmuth=wt*r$Dmuth
r$Dmu3 =  rep(0, length(mu))

r$Dmu3[d1] = 2*(((sqrt(2/pi)*(mu1-y1)^2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-5*ltheta))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))-(6*(mu1-y1)*exp(-exp(-2*ltheta)*(mu1-y1)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))+(4*sqrt(2)*exp(-3/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)) )
r$Dmu3[d2] =- 2*((sqrt(2/pi)*(y2-mu2)^2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-5*ltheta))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))-(6*(y2-mu2)*exp(-exp(-2*ltheta)*(y2-mu2)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))+(4*sqrt(2)*exp(-3/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3))
r$Dmu3=wt*r$Dmu3
r$Dmu2th = rep(2*( ( -2 )/theta^2), length(y))
r$Dmu2th[d2]=-2*(-(2*(y2-mu2)^2*exp(-(y2-mu2)^2/theta^2)/theta^4)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2)+(sqrt(2/pi)*(y2-mu2)*exp(-1/2*(y2-mu2)^2/theta^2)/theta^3*((y2-mu2)^2/theta^2-3))/(erfc(((y2-mu2)/theta)/sqrt(2)))+(4*sqrt(2)*(y2-mu2)*exp(-3/2*(y2-mu2)^2/theta^2)/theta^3)/(pi^(3/2)*erfc(((y2-mu2)/theta)/sqrt(2))^3)-(2*exp(-(y2-mu2)^2/theta^2)/theta^2*(2*(y2-mu2)^2/theta^2-2))/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2))
r$Dmu2th[d1]=-2*(-(2*(mu1-y1)^2*exp(-(mu1-y1)^2/theta^2)/theta^4)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2)+(sqrt(2/pi)*(mu1-y1)*exp(-1/2*(mu1-y1)^2/theta^2)/theta^3*((mu1-y1)^2/theta^2-3))/(erfc(((mu1-y1)/theta)/sqrt(2)))+(4*sqrt(2)*(mu1-y1)*exp(-3/2*(mu1-y1)^2/theta^2)/theta^3)/(pi^(3/2)*erfc(((mu1-y1)/theta)/sqrt(2))^3)-(2*exp(-(mu1-y1)^2/theta^2)/theta^2*(2*(mu1-y1)^2/theta^2-2))/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2))
r$Dmu2th=wt*r$Dmu2th
#r$Dmu2th = 2*( (!d1) * (!d2)*  ( -2 )/theta^2 -d1*((2*(mu-y)^2*exp(-exp(-2*ltheta)*(mu-y)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^2)+(sqrt(2/pi)*(mu-y)*exp(-1/2*exp(-2*ltheta)*(mu-y)^2-3*ltheta)*(exp(-2*ltheta)*(mu-y)^2-3))/(erfc((exp(-ltheta)*(mu-y))/sqrt(2)))+(4*sqrt(2)*(mu-y)*exp(-3/2*exp(-2*ltheta)*(mu-y)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^3)-(2*exp(-exp(-2*ltheta)*(mu-y)^2-2*ltheta)*(2*exp(-2*ltheta)*(mu-y)^2-2))/(pi*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^2))-d2*((2*(y-mu)^2*exp(-exp(-2*ltheta)*(y-mu)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^2)+(sqrt(2/pi)*(y-mu)*exp(-1/2*exp(-2*ltheta)*(y-mu)^2-3*ltheta)*(exp(-2*ltheta)*(y-mu)^2-3))/(erfc((exp(-ltheta)*(y-mu))/sqrt(2)))+(4*sqrt(2)*(y-mu)*exp(-3/2*exp(-2*ltheta)*(y-mu)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^3)-(2*exp(-exp(-2*ltheta)*(y-mu)^2-2*ltheta)*(2*exp(-2*ltheta)*(y-mu)^2-2))/(pi*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^2) ))
}
      if (level>1) { ## whole damn lot
r$Dmu4 = rep(0, length(mu))
#r$Dmu4[d1] = -2*((sqrt(2/pi)*(mu1-y1)^3*exp(-1/2*(mu1-y1)^2/theta^2)/theta^7)/(erfc(((mu1-y1)/theta)/sqrt(2)))-(14*(mu1-y1)^2*exp(-(mu1-y1)^2/theta^2)/theta^6)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2)-(3*sqrt(2/pi)*(mu1-y1)*exp(-1/2*(mu1-y1)^2/theta^2)/theta^5)/(erfc(((mu1-y1)/theta)/sqrt(2))) + (24*sqrt(2)*(mu1-y1)*exp(-3/2*(mu1-y1)^2/theta^2)/theta^5)/(pi^(3/2)*erfc(((mu1-y1)/theta)/sqrt(2))^3)+(8*exp(-(mu1-y1)^2/theta^2)/theta^4)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2)-(24*exp(-2*(mu1-y1)^2/theta^2)/theta^4)/(pi^2*erfc(((mu1-y1)/theta)/sqrt(2))^4))
r$Dmu4[d1] = -2*((sqrt(2/pi)*(mu1-y1)^3*exp(-1/2*(mu1-y1)^2/theta^2)/theta^7)/(erfc(((mu1-y1)/theta)/sqrt(2)))-(14*(mu1-y1)^2*exp(-(mu1-y1)^2/theta^2)/theta^6)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2)-(3*sqrt(2/pi)*(mu1-y1)*exp(-1/2*(mu1-y1)^2/theta^2)/theta^5)/(erfc(((mu1-y1)/theta)/sqrt(2))) + (24*sqrt(2)*(mu1-y1)*exp(-3/2*(mu1-y1)^2/theta^2)/theta^5)/(pi^(3/2)*erfc(((mu1-y1)/theta)/sqrt(2))^3)+(8*exp(-(mu1-y1)^2/theta^2)/theta^4)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2)-((sqrt(24)*exp(-1*(mu1-y1)^2/theta^2)/theta^2)/(pi*erfc(((mu1-y1)/theta)/sqrt(2))^2))^2)
#r$Dmu4[d2]= -2*((sqrt(2/pi)*(y2-mu2)^3*exp(-1/2*(y2-mu2)^2/theta^2)/theta^7)/(erfc(((y2-mu2)/theta)/sqrt(2)))-(14*(y2-mu2)^2*exp(-(y2-mu2)^2/theta^2)/theta^6)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2)-(3*sqrt(2/pi)*(y2-mu2)*exp(-1/2*(y2-mu2)^2/theta^2)/theta^5)/(erfc(((y2-mu2)/theta)/sqrt(2)))+(24*sqrt(2)*(y2-mu2)*exp(-3/2*(y2-mu2)^2/theta^2)/theta^5)/(pi^(3/2)*erfc(((y2-mu2)/theta)/sqrt(2))^3)+(8*exp(-(y2-mu2)^2/theta^2)/theta^4)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2)-(24*exp(-2 * (y2-mu2)^2/theta^2)/theta^4)/(pi^2*erfc(((y2-mu2)/theta)/sqrt(2))^4))
r$Dmu4[d2]= -2*((sqrt(2/pi)*(y2-mu2)^3*exp(-1/2*(y2-mu2)^2/theta^2)/theta^7)/(erfc(((y2-mu2)/theta)/sqrt(2)))-(14*(y2-mu2)^2*exp(-(y2-mu2)^2/theta^2)/theta^6)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2)-(3*sqrt(2/pi)*(y2-mu2)*exp(-1/2*(y2-mu2)^2/theta^2)/theta^5)/(erfc(((y2-mu2)/theta)/sqrt(2)))+(24*sqrt(2)*(y2-mu2)*exp(-3/2*(y2-mu2)^2/theta^2)/theta^5)/(pi^(3/2)*erfc(((y2-mu2)/theta)/sqrt(2))^3)+(8*exp(-(y2-mu2)^2/theta^2)/theta^4)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2)-((sqrt(24)*exp(-1 * (y2-mu2)^2/theta^2)/theta^2)/(pi*erfc(((y2-mu2)/theta)/sqrt(2))^2))^2)
r$Dmu4=wt*r$Dmu4

# 2*(d1*((sqrt(2/pi)*(mu-y)^3*exp(-1/2*exp(-2*ltheta)*(mu-y)^2-7*ltheta))/(erfc((exp(-ltheta)*(mu-y))/sqrt(2)))+(14*(mu-y)^2*exp(-exp(-2*ltheta)*(mu-y)^2-6*ltheta))/(pi*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^2)+(3*sqrt(2/pi)*(mu-y)*exp(-1/2*exp(-2*ltheta)*(mu-y)^2-5*ltheta))/(erfc((exp(-ltheta)*(mu-y))/sqrt(2)))-(24*sqrt(2)*(mu-y)*exp(-3/2*exp(-2*ltheta)*(mu-y)^2-5*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^3)-(8*exp(-exp(-2*ltheta)*(mu-y)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^2)+(24*exp(-2*exp(-2*ltheta)*(mu-y)^2-4*ltheta))/(pi^2*erfc((exp(-ltheta)*(mu-y))/sqrt(2))^4))+d2*((sqrt(2/pi)*(y-mu)^3*exp(-1/2*exp(-2*ltheta)*(y-mu)^2-7*ltheta))/(erfc((exp(-ltheta)*(y-mu))/sqrt(2)))+(14*(y-mu)^2*exp(-exp(-2*ltheta)*(y-mu)^2-6*ltheta))/(pi*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^2)+(3*sqrt(2/pi)*(y-mu)*exp(-1/2*exp(-2*ltheta)*(y-mu)^2-5*ltheta))/(erfc((exp(-ltheta)*(y-mu))/sqrt(2)))-(24*sqrt(2)*(y-mu)*exp(-3/2*exp(-2*ltheta)*(y-mu)^2-5*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^3)-(8*exp(-exp(-2*ltheta)*(y-mu)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^2)+(24*exp(-2*exp(-2*ltheta)*(y-mu)^2-4*ltheta))/(pi^2*erfc((exp(-ltheta)*(y-mu))/sqrt(2))^4)))
r$Dth2 = 2*( 2*(!d1) * (!d2)* (y-mu)^2 /theta^2)
r$Dth2[d1] =2*((2*(mu1-y1)^2*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)-(sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*(exp(-2*ltheta)*(mu1-y1)^2-1))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))))
r$Dth2[d2] = 2*((2*(y2-mu2)^2*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)-(sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*(exp(-2*ltheta)*(y2-mu2)^2-1))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))))
r$Dth2=wt*r$Dth2
r$Dmuth2 = 2*( -4*(!d1) * (!d2)* (y-mu) /theta^2)
r$Dmuth2[d1] = -2*((4*(mu1-y1)*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta)*(exp(-2*ltheta)*(mu1-y1)^2-1))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)-sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*((4*(mu1-y1)^2*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)-(sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*(exp(-2*ltheta)*(mu1-y1)^2-1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)+(sqrt(2/pi)*(2*(mu1-y1)^2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)-exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*(exp(-2*ltheta)*(mu1-y1)^2-1)^2))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))))
r$Dmuth2[d2] =2*((4*(y2-mu2)*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta)*(exp(-2*ltheta)*(y2-mu2)^2-1))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)-sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*((4*(y2-mu2)^2*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)-(sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*(exp(-2*ltheta)*(y2-mu2)^2-1))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)+(sqrt(2/pi)*(2*(y2-mu2)^2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)-exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*(exp(-2*ltheta)*(y2-mu2)^2-1)^2))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))))
r$Dmuth2=wt*r$Dmuth2
r$Dmu2th2= 2*( 4*(!d1) * (!d2) /theta^2)

r$Dmu2th2[d1]=2*((2*(4*exp(-exp(-2*ltheta)*(mu1-y1)^2-6*ltheta)*(mu1-y1)^2-2*exp(-exp(-2*ltheta)*(mu1-y1)^2-4*ltheta))*(mu1-y1)^2)/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)-(8*exp(-exp(-2*ltheta)*(mu1-y1)^2-4*ltheta)*((2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*sqrt(2/pi)*(mu1-y1)^2)/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3+(2*(mu1-y1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)*(mu1-y1))/pi+2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)*sqrt(2/pi)*((2*exp(-2*ltheta)*(mu1-y1)^2)/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))+(exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*sqrt(2/pi)*(exp(-2*ltheta)*(mu1-y1)^2-1)*(mu1-y1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2+(exp(-2*ltheta)*(mu1-y1)^2-1)/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))))*(mu1-y1)+(sqrt(2/pi)*(exp(-2*ltheta)*(mu1-y1)^2-1)*(exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)-exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-5*ltheta)*(mu1-y1)^2)*(mu1-y1))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))+(2*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta)*(((12*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^4)-(2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)*sqrt(2/pi)*(mu1-y1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)*(mu1-y1)^2+(8*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*sqrt(2/pi)*(mu1-y1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3+2/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2))/pi-exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*sqrt(2/pi)*((exp(-2*ltheta)*(mu1-y1)^2-1)*((4*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)-(exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)*sqrt(2/pi)*(mu1-y1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)*(mu1-y1)+(6*exp(-2*ltheta)*(mu1-y1))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))+(2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*sqrt(2/pi)*(3*exp(-2*ltheta)*(mu1-y1)^2-1))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2))
r$Dmu2th2[d2]=2*((2*(4*exp(-exp(-2*ltheta)*(y2-mu2)^2-6*ltheta)*(y2-mu2)^2-2*exp(-exp(-2*ltheta)*(y2-mu2)^2-4*ltheta))*(y2-mu2)^2)/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)-(8*exp(-exp(-2*ltheta)*(y2-mu2)^2-4*ltheta)*((2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*sqrt(2/pi)*(y2-mu2)^2)/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3+(2*(y2-mu2))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)*(y2-mu2))/pi+2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)*sqrt(2/pi)*((2*exp(-2*ltheta)*(y2-mu2)^2)/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))+(exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*sqrt(2/pi)*(exp(-2*ltheta)*(y2-mu2)^2-1)*(y2-mu2))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2+(exp(-2*ltheta)*(y2-mu2)^2-1)/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))))*(y2-mu2)+(sqrt(2/pi)*(exp(-2*ltheta)*(y2-mu2)^2-1)*(exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)-exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-5*ltheta)*(y2-mu2)^2)*(y2-mu2))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))+(2*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta)*(((12*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^4)-(2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)*sqrt(2/pi)*(y2-mu2))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)*(y2-mu2)^2+(8*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*sqrt(2/pi)*(y2-mu2))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3+2/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2))/pi-exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*sqrt(2/pi)*((exp(-2*ltheta)*(y2-mu2)^2-1)*((4*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)-(exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)*sqrt(2/pi)*(y2-mu2))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)*(y2-mu2)+(6*exp(-2*ltheta)*(y2-mu2))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))+(2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*sqrt(2/pi)*(3*exp(-2*ltheta)*(y2-mu2)^2-1))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2))
r$Dmu2th2=wt*r$Dmu2th2

r$Dmu3th = rep(0, length(mu))

r$Dmu3th[d1] = 2*(((sqrt(2/pi)*(mu1-y1)*((mu1-y1)^3*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-7*ltheta)-3*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-5*ltheta)))/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2)))-sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta)*(3*((4*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)-(sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)+(mu1-y1)*((sqrt(2/pi)*(mu1-y1)^2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-5*ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2-(12*(mu1-y1)*exp(-exp(-2*ltheta)*(mu1-y1)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2+(12*sqrt(2)*exp(-3/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^4)))+3*sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)*((mu1-y1)*((4*exp(-exp(-2*ltheta)*(mu1-y1)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^3)-(sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)+(2*sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2)+3*sqrt(2/pi)*(exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-3*ltheta)-(mu1-y1)^2*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-5*ltheta))*((sqrt(2/pi)*(mu1-y1)*exp(-1/2*exp(-2*ltheta)*(mu1-y1)^2-ltheta))/erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))^2+1/(erfc((exp(-ltheta)*(mu1-y1))/sqrt(2))))))


r$Dmu3th[d2] = -2*((sqrt(2/pi)*(y2-mu2)*((y2-mu2)^3*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-7*ltheta)-3*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-5*ltheta)))/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))-sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta)*(3*((4*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)-(sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)+(y2-mu2)*((sqrt(2/pi)*(y2-mu2)^2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-5*ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2-(12*(y2-mu2)*exp(-exp(-2*ltheta)*(y2-mu2)^2-4*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)-(sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2+(12*sqrt(2)*exp(-3/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/(pi^(3/2)*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^4)))+3*sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)*((y2-mu2)*((4*exp(-exp(-2*ltheta)*(y2-mu2)^2-2*ltheta))/(pi*erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^3)-(sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)+(2*sqrt(2/pi)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2)+3*sqrt(2/pi)*(exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-3*ltheta)-(y2-mu2)^2*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-5*ltheta))*((sqrt(2/pi)*(y2-mu2)*exp(-1/2*exp(-2*ltheta)*(y2-mu2)^2-ltheta))/erfc((exp(-ltheta)*(y2-mu2))/sqrt(2))^2+1/(erfc((exp(-ltheta)*(y2-mu2))/sqrt(2)))))
}
#for (i in 1:length(r)){
#r[[i]] = pmin(r[[i]], 1e200)
#r[[i]] = pmax(r[[i]], -1e200)
#}
r
} #end dd, thank god


    aic <- function(y, mu, theta=NULL, wt, dev) {
	right.threshold = get(".right.threshold")	
	left.threshold = get(".left.threshold")
	k = (y>left.threshold) + (y>= right.threshold)+1
      if (is.null(theta)) theta <- get(".Theta")
      lth = theta
      theta <- exp(theta) ## note log theta supplied

	2*sum(	(lth- log(dnorm((y - mu)[k==2]/theta)))*wt[k==2] ,  -( log(pnorm( (y - mu)[k == 1]/theta))) *wt[k==1],  - (log(pnorm( ( mu - y )[k==3]/theta))) * wt[k==3]) 
    }


    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik
       list(ls=0,## saturated log likelihood
            lsth1=0,  ## first deriv vector w.r.t theta - last element relates to scale
            lsth2=0) ##Hessian w.r.t. theta
     }


    initialize <- expression({
        n <- rep(1, nobs)
#print(y)
	y = drop(y)
if (((length(family$right.threshold) >  1) && (length(family$right.threshold) != length(y) )) || ((length(family$left.threshold) >  1) && (length(family$left.threshold) != length(y) ))
){
 warning("Mismatch between length of threshold argument and supplied number of observations. Ensure threshold arguments omit values corresponding to NA in dataset or preprocess data?")
}
#	censid = 
	y=pmax(y, family$left.threshold)
	y=pmin(y, family$right.threshold)
  	mustart <- rep(mean(drop(y)), length(y))
#print(family$getTheta(T))
#	if (family$n.theta > 0 ) # family$putTheta(log(sd(y[(y>family$left.threshold)& (y<family$right.threshold)])))
#	    assign(".Theta", log(sd(y[(y>family$left.threshold)& (y<family$right.threshold)]))  , envir = environment(family$putTheta))
    })

    rd <- function(mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      rnorm(mu,sd=Theta,mean=mu)
    }
    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      qnorm(p,sd=Theta,mean=mu)
    }


#   postproc <- expression({
#      object$family$family <-
#      paste("Censored normal (sigma=",round(object$family$getTheta(TRUE),3),")",sep="")
#      object$family$sigma.estimate =  object$family$getTheta(TRUE)
#      object$family$theta.estimate =  object$family$getTheta(FALSE)
#    })


#new postproc


    postproc <- function(family,...){#expression({

	posr <- list()
      posr$family <-
      paste("Censored normal (sigma=",round(family$getTheta(TRUE),3),")",sep="")
      attr(posr$family, "sigma.estimate") =  family$getTheta(TRUE)
      attr(posr$family, "theta.estimate") =  family$getTheta(FALSE)
	posr
    }#)

#     environment(dev.resids) <- environment(aic) <- environment(getTheta) <-
     environment(dev.resids) <-environment(Dd)<- environment(aic) <- environment(getTheta) <-
     environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) <- env
    structure(list(family = "censored normal", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf, left.threshold = left.threshold, right.threshold = right.threshold),
        class = c("extended.family","family"))
} ##

