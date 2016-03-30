## tsiR 
## an R package for running TSIR models

### Version
0.0.1.0

### Installation

You can currently install **tsiR** via devtools which will install all dependencies, except rjags which must be install manually. 
```sh
require(devtools)
install_github('adbecker/tsiR')
```
Load the package and dependencies **kernlab**, **ggplot2**, **reshape2**, **grid**.

```sh
require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
```

### Data 

We include 20 data sets form the UK each with 20 years of measles data (biweekly data, i.e. IP = 2 in the *runtsir*, *estpars*, *simulatetsir*, *mcmctsir*, *mcmcestpars* functions). 

We require that the data must be a data.frame with column names 'time','cases', 'births', and 'pop'. Data can be plotted using *plotdata* or *plotcases*.

### Example 1
First load and plot the data:
```sh
names(twentymeas)
LondonMeas <- twentymeas[["London"]]
head(LondonMeas)
plotdata(LondonMeas)
```

The default settings for *runtsir* (and all functions) can be accessed through *runtsir* in the R console. Based on measles, everything is written in modulo two weeks, however this can be changed by setting *IP=* to the number of weeks between each time step. The minimal necessary input for a biweekly data set is thus simply *data = *. Here we run a simple example using the default options: cumulative cases on the x axis, a Gaussian regression (default) between cumulative cases and births, estimating both *Sbar* and *alpha*, estimating a 26 (52/IP) point contact parameter, and we run the forward simulation completely forward. We specify the option to draw the next time step from a negative binomial distribution as well as specifying a Poisson link function in the GLM. Additionally, we specify to do 100 simulations.

```sh
basictsir <- runtsir(LondonMeas,method='negbin',
                     family='poisson',nsim=100)
```

We plot the full diagnostic using the *plotres* function. 

```sh
plotres(basictsir)
```
Walking through the eight plots: the first gives the regression model *Yhat* between cumulative births and cases, the second gives the change in reporting rate over time, the third gives the residuals and then the fitted susceptible dynamics, the fourth is the profiled susceptibility, the fifth is the repeating seasonal contact parameter with 95% confidence intervals, the sixth is the predicted-observed regression, and the final two are the forward predictions (red) and the data (blue). The top forward prediction shows each of the *nsim* individually, while the bottom shows the mean with 95\% confidence intervals.

### Example 2
ex <- exampletsiR()


