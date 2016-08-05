[![Build Status](https://travis-ci.org/adbecker/tsiR.svg?branch=master)](https://travis-ci.org/adbecker/tsiR)

## tsiR 
## an R package for fitting TSIR models 

## Introduction 

This R package fits and runs the TSIR model using a number of different options, both Bayesian and frequentist. 

### Version
0.0.1.0

### Installation

You can currently install **tsiR** via **devtools** (although the package will soon be put on CRAN) which will install all dependencies, except rjags which must be installed manually (https://cran.r-project.org/web/packages/rjags/index.html). 
```sh
require(devtools)
install_github('adbecker/tsiR')
```
Load the package and dependencies **kernlab**, **ggplot2**, **reshape2**, and **grid**, The first dependency is for the Gaussian Process regression, while the others are for plotting. 

```sh
require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
```

### Data 

We include 20 data sets form the UK each with 20 years of measles data (biweekly data, i.e. infectious period of two weeks (```IP = 2```) in the ```runtsir```, ```estpars```, ```simulatetsir```, ```mcmctsir```, ```mcmcestpars``` functions). 

We require that the data must be a ```data.frame``` with column names 'time','cases', 'births', and 'pop'. You can load these four vectors into the function ```tsiRdata``` (``` example <- tsiRdata(time = , cases = , births = , pop = , IP = )```) where ```IP=``` designates the generation time to interpolate the data on. Data can be plotted using ```plotdata``` or ```plotcases```. If using any of the plot functions, the time column of the data.frame must be recognized by **ggplot2**.

### Main Functions

A table of the main functions can be found below. Individual options can be found via ```?function_name``` in the R console.

| Function | Purpose |
|----------|-----------|
|```estpars``` |  Using simple regression: reconstructs susceptibles and estimates parameters|
|```exampletsiR``` |  Using the London data: walks one through the fitting options of ```runtsir```|
|```maxthreshold``` |  Using simple regression: optimizes the threshold parameter for sparse data|
|```mcmcestpars``` |  Using MCMC: reconstructions susceptibles and estimates parameters|
|```mcmctsir``` |  Using MCMC: reconstructions susceptibles, estimates parameters, and runs the simulation|
|```plotcases``` |  Plots just the cases data|
|```plotcomp``` |  Plots the data versus the simulations|
|```plotdata``` |  Plots the (interpolated) cases, births, and population time series|
|```plotres``` |  Plots the fitted regressions, reporting, susceptible reconstruction, estimated parameters, and the fit diagonstics|
|```runtsir``` |  Using simple regression: reconstructions suceptibles, estimates parameters, and runs the simulation|
|```simulatetsir``` |  Runs the simulation taking in the output from ```estpars``` or ```mcmcestpars```|
|```tsiRdata``` |  Interpolates cases, births, and pop vectors to the generation time of the disease|
|```twentymeas``` |  Complete biweekly (IP=2) data sets from twenty UK cities|

In the fitting functions, a number of options such as fixing alpha (```alpha = 0.97```, for example) and Sbar (```sbar = 0.025```, for example) are provided. Additionally, contact can be seasonal, driven by the school term calendar, or a single non time-varying parameter (```seasonality = "standard", "schoolterm" or "none")``` respectively. These functions work for data with any fixed generation time, although in the examples we will use biweekly data ```(IP = 2)```.

Several examples follow, although they do not completely exhaust the fitting options.

### Example 1

After loading the dependencies, load and plot the data:

```sh
names(twentymeas)
LondonMeas <- twentymeas[["London"]]
head(LondonMeas)
plotdata(LondonMeas)
```

The default settings for ```runtsir``` (and all functions) can be accessed through ```?runtsir``` in the R console. Based on measles, everything is written in modulo two weeks, however this can be changed by setting ```IP=``` to the number of weeks between each time step. The minimal necessary input for a biweekly data set is thus simply ```data = ```. Here we run a simple example using the default options: cumulative cases on the x axis, a Gaussian regression (default) between cumulative cases and births, estimating both ```sbar``` and ```alpha```, estimating a 26 (52/IP) point contact parameter, and we run the forward simulation completely forward. We specify the option to draw the next time step from a negative binomial distribution and we specificy a Gaussian family (default) with an identity link (default) function in the GLM. Any family and link can be inputted, however, the options are essentially ```quasipoisson, poisson, gaussian``` where ```poisson/ quassipoisson``` take ```link=log``` and ```family=gaussian``` takes either a ```log``` or ```identity``` link.

Options for the GLM family 

Additionally, we specify to do 100 simulations. Other regression types can be specified using ```regtype=``` where the options are ```lm, lowess, loess, spline``` for a linear regression, a lowess, a loess, and a spline regression with 2.5 degrees freedom.

```sh
LondonRes <- runtsir(data=LondonMeas,IP=2,method='negbin',regtype='gaussian',
                     family='gaussian',link='identity',nsim=100)
```

We plot the full diagnostic using the ```plotres``` function. 

```sh
plotres(LondonRes)
```
Walking through the eight plots: the first gives the regression model ```Yhat``` between cumulative births and cases, the second gives the reporting rate over time (for ```regtype = 'lm'``` this is a constant), the third gives the residuals and then the reconstructed susceptible dynamics, the fourth is the profiled average susceptibility, the fifth is the repeating seasonal contact parameter with 95% confidence / credible intervals, the sixth is the predicted-observed regression, and the final two are the forward predictions (red) and the data (blue). The top forward prediction shows each of the ```nsim``` individually, while the bottom shows the mean with 95% confidence intervals.

### Example 2

We can run this same model in a Bayesian framework with schooltime forcing. We can specify almost all the same functions as the ```estpars``` / ```runtsir``` functions, but also number of iterations, chains, burnins, adapts, etc. This example takes about 30-45 minutes to run. Run time scales with the number of parameters, thus ```seasonality = 'standard'``` will take longer. We can plot the traces and posteriors.

```sh
require(rjags)
parms <- mcmcestpars(data=LondonMeas,IP=2,seasonality='schoolterm',update.iter = 1e4,n.iter=1e5,n.chains=3)
names(parms)
plot(parms$mcmcsamples)
```

We can now enter these parameters into the ```simulatetsir``` function. Here we will use a Poisson distribution to draw the next generation and do twenty simulations. Note that this model gives similar results, and 95% credible intervals are constructed for the parameter estimates themselves. 

```sh
sim <- simulatetsir(data=LondonMeas,parms=parms,method='pois',nsim=20)
plotres(sim)
```

Alternatively, we could run the parameter estimation and forward simulations all in one, by calling the ```mcmctsir``` function.

### Example 3

Often times however, the population size is not large enough to support a persistent epidemic like in London. These epidemics are often sporadic and may or may not be spatially coupled. An example of this is the city of Mold.

```sh
MoldMeas <- twentymeas[["Mold"]]
plotdata(MoldMeas)
```

Note that while there are some small case counts, there are essentially only outbreaks every year or two. These epidemics cannot be predicted fully forward, but rather must be done epidemic-ahead (Caudron et al 2015). In order to do this, a threshold parameter must be set which indicts when a new epidemic has been sparked. This threshold parameter can be set arbitrarily, i.e. ```runtsir(data = MoldMeas, epidemics = 'break', threshold = 10)``` for example. Alternatively, the threshold can be fit by maximizing Rsquared predictions using the ```maxthreshold``` function. Epidemic-ahead predictions can then be accomplished using the ```epidemics = break``` command in the ```runtsir``` function. 

First, parameter estimates must be acquired using either ```estpars``` or ```mcmcestpars```. While not necessary, here we fix ```alpha = 0.97```. The average number of susceptibles could be fixed the same way, i.e. ```sbar = ```.

```sh
parms <- estpars(data=MoldMeas,alpha=0.97,IP=2)
tau <- maxthreshold(data=MoldMeas,parms=parms,IP=2)
MoldRes <- simulatetsir(data=MoldMeas,parms=parms,epidemics='break',threshold=tau,
                            nsim=10,method='negbin')
plotres(MoldRes)
```
### Final note

Note that all of these examples used biweekly data, ```IP = 2```. If the data is aggregated by a different number of weeks, you must specificy this globally and consistently.

### Issues

Please email me at adbecker at princeton dot edu with any issues or suggestions. If you find the package useful, please cite us. The bibtex file can be called by ```citation('tsiR')``` in the R console.








