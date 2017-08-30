\name{cenGAM-package}
\alias{cenGAM-package}
\alias{cenGAM}
\docType{package}
\title{Censored response additive modelling for mgcv}
\description{
\code{cenGAM} provides {\code{\link{tobit1}}} and {\code{\link{tobit2}}} families for generalized additive modelling with the \code{gam} function in the \code{mgcv} package.

Under the Tobit I model, the user supplies a left and/or right threshold (optionally differing between observations, and potentially infinite to denote no censorship) and the response is assumed to be censored if it falls over this threshold. Under the Tobit II model, a more generalised model is assumed where we fit a second additive model that gives whether each observation is censored, with a possible correlation between the error in the response and the censorship.

Seen help for the individual functions for more details.
}
\author{
Zhou Fang <zhou.fang@bioss.ac.uk>

with contributions and help from Javier Palarea.
}
\keyword{ package }
\keyword{models} 
\keyword{smooth} 
\keyword{regression}


\examples{
## see examples for tobit1 and tobit2
}
