\name{viewDoc}
\alias{viewDoc}
\title{View documentation}
\description{
  View bayespack documentation, which contains a brief description of
  the implemented methods. The viewDoc function is exactly the same as
  the faux.vignette function from Roger Koenker's quantreg package.
}
\usage{
  viewDoc(topic)
}
\arguments{
  \item{topic}{ Which documentation pdf file to show. When missing the
    list of documentation pdfs is shown. }
}
\examples{
  ## Documents available
  ## (currently only original BAYESPACK Fortran manual, which provides
  ## a nice overview of the implemented methods)
  viewDoc()
}
\keyword{ distribution }

