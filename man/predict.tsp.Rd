\name{predict.tsp}
\alias{predict.tsp}
\alias{predict}
\title{Prediction based on a tsp object}
\description{
This function can be used to predict outcome values for a data set
based on a tsp object.  
}
\usage{
  \method{predict}{tsp}(object,dat=NULL,select=NULL,...)
}
\arguments{
  \item{object}{A tsp object}
  \item{dat}{Can take two values: (a) an m genes by n arrays matrix of
  expression data or (b) an eSet object}
  \item{select}{An indicator of which TSP to use,defaults to the first
  TSP.}
  \item{...}{Plotting arguments (ignored)}	     	     
}
\details{
  predict() accepts a tsp object calculated on an expression set
  or gene expression matrix. If no other data set is included, the tsp
  predictions for the original data set are produced. If a second gene
  expression matrix or expression set is included, predict() looks
  for the gene names of the TSP in tspobj and attempts to match them
  in the rownames or featureNames of the gene expression matrix. If
  rownames or featureNames are not available, the prediction is based
  on the row numbers.  If a match is identified, predict() makes a
  prediction for each gene based on the output.
}
\value{
  \item{predict}{A class prediction for each array of dat based on the
  TSP from tspobj}
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.

  A.C. Tan, D.Q. Naiman, L. Xu, R.L. Winslow, D. Geman, "Simple decision rules for classifying human cancers from gene expression profiles,"
  Bioinformatics, 21:3896-3904, 2005. 
}
\author{Jeffrey T. Leek \email{jtleek@jhu.edu}}
\seealso{\code{\link{tspplot}},\code{\link{ts.pair}}, \code{\link{tspcalc}},\code{\link{tspsig}},\code{\link{summary.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 

  ## Run tspcalc() on a data matrix and grp vector
  tsp1 <- tspcalc(dat,grp)

  ## Get predictions for a new eSet or data matrix
  predict.tsp(tsp1,dat2,1)
  predict(tsp1,eSet2,1)
 }
}

\keyword{misc}
