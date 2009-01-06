\name{tspsig}
\alias{tspsig}
\title{Significance calculation for top scoring pairs}
\description{
This function calculates the significance of a top-scoring pair. It
can be run after tspcalc() to calculate how strong a TSP is. 
}
\usage{
  tspsig(dat,grp,B=50,seed=NULL)
}
\arguments{
  \item{dat}{Can take two values: (a) an m genes by n arrays matrix of expression data or (b) an eSet object}
  \item{grp}{Can take one of two values: (a) A group indicator incharacter or numeric form, (b) an integer indicating the column of
  pData(dat) to use as the group indicator}
  \item{B}{The number of permutations to perform in calculation of the
  p-value, default is 50.}
  \item{seed}{If this is a numeric argument, the seed will be set for
  reproducible p-values.}
}
\details{
  tspsig() only works for two group classification. The computation time grows rapidly in the number of 
  genes, so for large gene expression matrices one should be prepared
  to wait or do a pre-filtering step. A progress bar is shown which
  gives some indication of the time until the calculation is complete.
  The top scoring pairs methodology was originally described in Geman et al. (2004).
}
\value{
  \item{p}{A p-value for testing the null hypothesis that there is no TSP
  for the data set dat.}
  \item{nullscores}{The null TSP scores from the permutation test.}
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.
}

\author{Jeffrey T. Leek \email{jtleek@jhu.edu}}
\seealso{\code{\link{tspplot}}, \code{\link{ts.pair}}, \code{\link{tspcalc}},\code{\link{predict.tsp}}, \code{\link{summary.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 

  ## Run tspcalc() on a data matrix and grp vector
  tsp1 <- tspcalc(dat,grp)

  ## Run tspsig() to get a p-value
  p <- tspsig(dat,grp)
  p
 }
}

\keyword{misc}
