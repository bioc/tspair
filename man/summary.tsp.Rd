\name{summary.tsp}
\alias{summary.tsp}
\alias{summary}
\title{Summary of a tsp object}
\description{
This function can be used to summarize a TSP object.   
}
\usage{
  \method{summary}{tsp}(object,select=NULL,printall=FALSE,...)
}
\arguments{
  \item{object}{A tsp object}  	  
  \item{select}{A numerical indicator of which TSP to summarize, if NULL, all TSPs are summarized one at a time}
  \item{printall}{If printall=TRUE all TSPs are summarized at once}  
  \item{...}{Summary arguments (ignored)} 
}
\details{
  summary() accepts a tsp object calculated on an expression set
  or gene expression matrix. The result is a table of frequencies for
  the TSP indicator versus the group indicator. 
}
\value{
  A table. 
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.
}
\author{Jeffrey T. Leek \email{jtleek@jhu.edu}}
\seealso{\code{\link{tspplot}}, \code{\link{ts.pair}}, \code{\link{tspcalc}},\code{\link{tspsig}}, \code{\link{predict.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 

  ## Run tspcalc() on a data matrix and grp vector
  tsp1 <- tspcalc(dat,grp)

  ## Get the summary for the tsp object. 
  summary(tsp1)
  summary.tsp(tsp1)
 }
}

\keyword{misc}
