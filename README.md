## tsiR 
## an R package for fitting TSIR models 

## Introduction 

Currently fits and runs the TSIR model using a number of different fitting options, both Bayesian and frequentist. Soon will be updated with TSIR(S) functions.

### Version
0.0.1.0

### Installation

You can currently install **tsiR** via devtools (will soon be put on CRAN) which will install all dependencies, except rjags which must be install manually. 
```sh
require(devtools)
install_github('adbecker/tsiR')
```
Load the package and dependencies **kernlab** **ggplot2** **reshape2** **grid** The first dependency is for the Gaussian Process regression, while the others are just for plotting. 

```sh
require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
```

### Data 

We include 20 data sets form the UK each with 20 years of measles data (biweekly data, i.e. ```IP = 2``` in the ```runtsir```, ```estpars```, ```simulatetsir```, ```mcmctsir```, ```mcmcestpars``` functions). 

We require that the data must be a data.frame with column names 'time','cases', 'births', and 'pop'. You can load these four vectors into the function ```tsiRdata``` (``` example <- tsiRdata(time = , cases = , births = , pop = , IP = )```) where ```IP=``` designates the generation time to interpolate the data on.

Data can be plotted using ```plotdata``` or ```plotcases```.

### Main Functions

A quick table of the main functions can be found below. Individual options can be found via ?function_name in R.

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

In the fitting functions, number of options such as fixing alpha (```alpha = 0.97```, for example) and Sbar (```sbar = 0.025```, for example) are provided. Additionally, contact can be seasonal or driven by the school term calendar (```seasonality = "standard", "schoolterm" or "none")``` for the null model. These functions work for data with any fixed generation time, although in the examples we will use biweekly data.

A couple examples (but not using all functions or options) follow.

### Example 1

After loading the dependencies, load and plot the data:
```sh
names(twentymeas)
LondonMeas <- twentymeas[["London"]]
head(LondonMeas)
plotdata(LondonMeas)
```

The default settings for ```runtsir``` (and all functions) can be accessed through ```runtsir``` in the R console. Based on measles, everything is written in modulo two weeks, however this can be changed by setting ```IP=``` to the number of weeks between each time step. The minimal necessary input for a biweekly data set is thus simply ```data = ```. Here we run a simple example using the default options: cumulative cases on the x axis, a Gaussian regression (default) between cumulative cases and births, estimating both ```sbar``` and ```alpha```, estimating a 26 (52/IP) point contact parameter, and we run the forward simulation completely forward. We specify the option to draw the next time step from a negative binomial distribution as well as specifying a Poisson link function in the GLM. Additionally, we specify to do 100 simulations. Other regression types can be specified using ```regtype=``` where the options are ```lm, lowess, loess, spline``` for a linear regression, a lowess, a loess, and a spline regression with 2.5 degrees freedom.

```sh
basictsir <- runtsir(data=LondonMeas,IP=2,method='negbin',
                     family='poisson',nsim=100)
```

We plot the full diagnostic using the ```plotres``` function. 

```sh
plotres(basictsir)
```
Walking through the eight plots: the first gives the regression model ```Yhat``` between cumulative births and cases, the second gives the change in reporting rate over time, the third gives the residuals and then the fitted susceptible dynamics, the fourth is the profiled susceptibility, the fifth is the repeating seasonal contact parameter with 95% confidence intervals, the sixth is the predicted-observed regression, and the final two are the forward predictions (red) and the data (blue). The top forward prediction shows each of the ```nsim``` individually, while the bottom shows the mean with 95% confidence intervals.

### Example 2

We can run this same model in a Bayesian framework with schooltime forcing. We can specify pretty much all the same functions as the ```estpars``` / ```runtsir``` function, but also number of iterations, chains, burnins, adapts, etc. This takes about 30-45 minutes to run.

```sh
require(rjags)
parms <- mcmcestpars(data=LondonMeas,IP=2,seasonality='schoolterm',update.iter = 1e4,n.iter=1e5,n.chains=3)
plot(parms$mcmcsamples)
```

We can now enter these parameters into the ```simulatetsir``` function. Here we will use a Poisson distribution to draw the next generation and do twenty simulations. Note that this model gives similar results, and 95% credible intervals are constructed for the paramter estimates themselves.

```sh
sim <- simulatetsir(data=LondonMeas,parms=parms,method='pois',nsim=20)
plotres(sim)
```

### Example 3

Often times however, the population size is not large enough to support a persistent epidemic like in London. These epidemics are often sporadic and may or may not be spatially coupled. An example of this is the city of Mold.

```sh
MoldMeas <- twentymeas[["Mold"]]
plotdata(MoldMeas)
```

Note that while there are some small case counts, there are essentially only outbreaks every year or two. These epidemics cannot be predicted fully forward, but rather must be done epidemic ahead (Caudron et al 2015). In order to do this, a threshold parameter must be set which indicts when a new epidemic has been sparked. This can be accomplished using the ```maxthreshold``` function, which maximizes the Rsq fit. Epidemic ahead predictions can then be accomplished using the ```break``` command in the ```runtsir``` function. Alternatively, the threshold paramter can be set to any integer without being maximized.

First, parameter estimates must be acquired using either ```estpars``` or ```mcmcestpars```. While not necessary, here we fix ```alpha = 0.97```. The average number of susceptibles can be fixed the same way, i.e. ```sbar = ```.

```sh
parms <- estpars(data=MoldMeas,alpha=0.97,IP=2)
tau <- maxthreshold(data=MoldMeas,parms=parms,IP=2)
MoldRes <- simulatetsir(data=MoldMeas,parms=parms,epidemics='break',threshold=tau,
                            nsim=10,method='negbin')
plotres(MoldRes)
```
### Final note

Remember all of these examples used ```IP = 2```, if your data is aggregated by one week, three weeks, etc, you have to specify this consistently or globally.

### Issues

Please email me at adbecker at princeton dot edu with any issues. If you find the package helpful or use it please cite us ```citation('tsiR')```.








