
\name{twentymeas}
\alias{twentymeas}
\docType{data}
\title{Measles incidence data from 20 cities
}
\description{
twentymeas is a list containing 20 dataframes with cases, births, populations. Each dataframe is a 22 year time series at biweekly (i.e. IP=2) intervals.
}
\usage{data("twentymeas")}


\source{
From Bryan Grenfell}

\examples{
names(twentymeas)
london <- twentymeas[["London"]]
plotdata(london)}
\keyword{datasets}
