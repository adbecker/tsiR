# tsiR
tsiR: A package to run time-series Susceptible-Infected-Recovered models in R. Soon to be on CRAN!

Example 1.

require(devtools)
install_github('adbecker/tsiR')
require(tsiR)
require(kernlab)
require(ggplot2)
require(grid)
require(reshape2)

data <- twentymeas[['London']]
plotdata(data)

res <- runtsir(data,method = 'negbin')

plotres(res)

Example 2.

ex <- exampletsiR()


