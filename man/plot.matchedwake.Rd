\name{plot.matchedwake}
\alias{plot.matchedwake}
\title{
plot function for objects of class \code{matchedwake}
}
\description{
Overloads the default \code{plot()} for objects of class \code{matchedwake}. Returns a contour plot: The lighter the color the larger the estimated treatment effect. The corresponding standard errors are indicated by shading out some of the estimates: No shading corresponds to \eqn{p<alpha1} for the treatment effect in the diff-in-diffs analysis. Dotted lines indicate p-values between \code{alpha1} and \code{alpha2} and full lines indicate \eqn{p>alpha2}. The cells indicating effect size and significance level are arranged in a table where each field corresponds to one specific combination of spatial and temporal sizes.
}
\usage{
\method{plot}{matchedwake}(x, ...)
}
\arguments{
  \item{x}{
  object of class \code{\link{matchedwake}}.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\references{
Schutte, S., Donnay, K. (2014). ``Matched wake analysis: Finding causal relationships in spatiotemporal event data.'' \emph{Political Geography} 41:1-10.
}
\author{
  Sebastian Schutte and Karsten Donnay.
}
\seealso{
  \code{\link{matchedwake}}
}