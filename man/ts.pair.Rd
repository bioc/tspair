\name{ts.pair}
\alias{ts.pair}
\title{Calculation of top scoring pairs.}
\description{
This function calculates the pair of genes that show the maximum
difference in ranking between two user specified groups. The function
tspcalc() is the general version of this function that accepts eSets.
}
\usage{
  ts.pair(dat,grp)
}
\arguments{
  \item{dat}{An m genes by n arrays matrix of expression data.}
  \item{grp}{A group indicator with values 0 and 1 for the two groups.}
}
\details{
  ts.pair() only works for two group classification. All the caveats
  for tspcalc() apply here, but this function only works for matrix
  and vector arguments, for a more general function see tspcalc().
  The top scoring pairs methodology was originally described in Geman et al. (2004).
  and the unique TSP idea was described in Tan et al. (2005).
}

\value{
  A tsp object with elements:
  \item{index}{A two-column matrix where each row is a pair of indices from the data matrix achieving the top score}
  \item{tspscore}{The rank based score described in Geman et al. (2004), essentially the empirical average of sensitivity and specificity for the pair.} 
  \item{score}{The tie-breaking score described in Tan et al. (2005).}
  \item{grp}{The binary group indicator.}
  \item{tspdat}{Row i and row (i + Number of total TSPs) of this data matrix represent the expression data for the ith TSP pair.}
  \item{labels}{The group labels for the two groups as defined by the grp variable.}
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.

  A.C. Tan, D.Q. Naiman, L. Xu, R.L. Winslow, D. Geman, "Simple decision rules for classifying human cancers from gene expression
  profiles," Bioinformatics, 21: 3896-3904, 2005.
}

\author{Jeffrey T. Leek \email{jtleek@jhu.edu}}
\seealso{\code{\link{tspplot}}, \code{\link{tspcalc}}, \code{\link{tspsig}},\code{\link{predict.tsp}}, \code{\link{summary.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 
  tsp1 <- ts.pair(dat,grp)
  tsp1$index
 }
}

\keyword{misc}
