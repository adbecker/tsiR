[![Build Status](https://travis-ci.org/adbecker/tsiR.svg?branch=master)](https://travis-ci.org/adbecker/tsiR)

## tsiR: An R package for time-series Susceptible-Infected-Recovered models of epidemics

## Information regarding this package and a short tutorial can be found here: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185528

## If you've found the package useful in your research, I ask that you cite the above PLOS ONE paper.

## This package can be installed via CRAN. Updates to the package post the PLOS ONE paper will follow below as they occur.

## Note : 03/15/2019
#### This warning message was added to V.0.4.1 but worth adding here as well -- if you find very unreasonable reporting rates along either endpoints or in highly variable epidemic regions using a gaussian regression, it may be worth reducing 'sigmamax' in either runtsir or estpars away from the default of 3 close to 0.5 or so. 

## Update V.0.4.1 : Minor update  01/29/2019
#### The tsiR package has been updated to include further warning messages, bug fixes, and further annotations.

## Update V.0.4.0 : Lyapunov Analysis 08/20/2018
#### The tsiR package has been updated to include Global and Local Lyapunov Exponents. You can learn more and find examples of how to use this function for the London data by typing ?TSIR_LE ?TSIR_LLE and ?plotLLE in the R console. This updateis now on CRAN.
