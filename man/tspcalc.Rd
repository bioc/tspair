\name{tspcalc}
\alias{tspcalc}
\title{Estimate top scoring pairs from a gene expression matrix}
\description{
This function calculates the pair of genes that show the maximum difference in ranking between two user specified groups. The "top scoring pair" (TSP) maximizes the average of sensitivity and specificity over all rank based classifiers using a pair of genes in the data set.
}
\usage{
  tspcalc(dat,grp)
}
\arguments{
  \item{dat}{Can take two values: (a) an m genes by n arrays matrix of expression data or (b) an eSet object}
  \item{grp}{Can take one of two values: (a) A group indicator in character or numeric form, (b) an integer indicating the column of pData(dat) to use as the group indicator.}
}
\details{
  tspcalc only works for two group classification. The computation time grows rapidly in the number of 
  genes, so for large gene expression matrices one should be prepared to wait or do a pre-filtering step.
  The top scoring pairs methodology was originally described in Geman et al. (2004) and the unique TSP
  idea was described in Tan et al. (2005). 
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
\seealso{\code{\link{tspplot}}, \code{\link{ts.pair}}, \code{\link{tspsig}},\code{\link{predict.tsp}}, \code{\link{summary.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 

  ## Run tspcalc() on a data matrix and grp vector
  tsp1 <- tspcalc(dat,grp)
  tsp1$index

  ## Run tspcalc() on an expression set and a column of the pData matrix
  tsp2 <- tspcalc(eSet1,1)
  tsp2$index
 }
}

\keyword{misc}
