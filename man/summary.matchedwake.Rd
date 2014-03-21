\name{summary.matchedwake}
\alias{summary.matchedwake}
\title{
summary function for objects of class \code{matchedwake}
}
\description{
Overloads the default \code{summary()} for objects of class \code{matchedwake}.
}
\usage{
\method{summary}{matchedwake}(object, detailed = FALSE, ...)
}
\arguments{
  \item{object}{
  object of class \code{\link{matchedwake}}.
}
  \item{detailed}{
  Boolean indicating whether or not a detailed summary should be returned. \code{Default = TRUE}.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\value{
Returns a \code{data.frame} with an overview of all significant results (significance level is \kbd{alpha1} as retrieved from \code{x$parameters}). If \code{detailed = TRUE} this overview includes a number of matching statistics and statistics on overlaps of the spatiotemporal cylinders.
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