\name{print.matchedwake}
\alias{print.matchedwake}
\title{
print function for objects of class \code{matchedwake}
}
\description{
Overloads the default \code{print()} for objects of class \code{matchedwake}.
}
\usage{
\method{print}{matchedwake}(x, ...)
}
\arguments{
  \item{x}{
  object of class \code{\link{matchedwake}}.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\value{
Returns a \code{data.frame} with all significant results (significance level is \kbd{alpha1} as retrieved from \code{x$parameters}).
}
\references{
Schutte, S., Donnay, K. (2014). ``Matched Wake Analysis: Identifying Causal Relationships in Spatiotemporal Event Data.'' Forthcoming in \emph{Political Geography}.
}
\author{
  Sebastian Schutte and Karsten Donnay.
}
\seealso{
  \code{\link{matchedwake}}
}