\name{tspplot}
\alias{tspplot}
\title{Graphical display of tsp objects}
\description{
  Graphical display of tsp objects
}
\usage{
  tspplot(tspobj)
}
\arguments{
  \item{tspobj}{A tsp object.}
}
\details{
  tspplot() creates a plot of each top scoring pair.The figure plots 
  the expression for the first gene in the TSP pair versus the expression 
  for the second gene in the TSP pair across arrays. The user defined 
  groups are plotted in the colors red and blue. The score for the pair is
  shown across the top of each plot. If there is more than one TSP, 
  hitting return will cycle from one TSP to the next. 
}
\value{
  Nothing of interest. 
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene
  expression profiles from pairwise mRNA comparisons, Statist. Appl. in Genetics and Molecular Biology, 3, 2004.
}
\author{Jeffrey T. Leek \email{jtleek@jhu.edu}}
\seealso{\code{\link{tspcalc}}, \code{\link{ts.pair}}, \code{\link{tspsig}},\code{\link{predict.tsp}}, \code{\link{summary.tsp}}}
\examples{
  \dontrun{
  ## Load data
  data(tspdata) 
  tsp1 <- tspcalc(dat,grp)
  tspplot(tsp1)
 }
}
\keyword{misc}
