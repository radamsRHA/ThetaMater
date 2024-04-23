# ThetaMater: Rapid and scalable Bayesian estimation of the population size parameter theta=4Nu from diverse genetic data
**NOTE See the file https://github.com/radamsRHA/ThetaMater/tree/master/vignettes for detailed instructions**

Steps for Bayesian estimation of theta and alpha (see tutorial for more details)

1. Read input files and convert to infinite-sites data used by ThetaMater 
2. Choose a specific Bayesian model for estimating theta along (`ThetaMater.M1`), theta with a fixed alpha shape parameter (`ThetaMater.M2`), or the joint posterior distribution of theta and alpha (`ThetaMater.M3`).
3. (Optional) Convert the results of ThetaMater into estimates of effective population size by multiplying the posterior distribution of theta by a given mutation rate.
4. (Optional) Filter out loci with high mutation counts that likely represent spurious loci using the posterior predictive simulator function (`ThetaMater.PPS`).
5. (Optional) Identifying the maximum bounds of mutation counts from the posterior distribution, remove loci with outlier mutation counts that are too high beyond the posterior distribution, and re-estimate Theta using filtered dataset

## Installing R package GppFst from github
*** IMPORTANT: Please make sure that the most recent version of Rcpp is installed. You can use this command to install it: `install_github("https://github.com/RcppCore/Rcpp")`
***

*** ThetaMater require the Boost c++ libraries to be installed for the fast calculation of the underlying likelihood functions. I recommend install the Boost libraries using brew for Mac systems (brew install boost). The latest version of the Boost c++ libraries here:
https://www.boost.org. You can also install boost libraries by installing brew and using the command `brew update` and then `brew install boost` 
*** 

The R package ThetaMater is freely available to download and distribute from github <https://github.com/radamsRHA/ThetaMater/>. To install and load ThetaMater, you must first install the R packages `devtools`, `Rcpp`, and `MCMCpack`. Additionally, make sure the most updated version of R version 3.3.3 is installed (see above warning). Download R version 3.3.3 here: https://cran.r-project.org/bin/macosx/.
```
install.packages("devtools")
install.packages("MCMCpack")
install_github("https://github.com/RcppCore/Rcpp")
```
Now using devtools we can install `ThetaMater` from github:

```
library(devtools)
install_github("radamsRHA/ThetaMater")
library(ThetaMater) # Load package ThetaMater
library(MCMCpack) # Load dependancy phybase
```


To begin using `ThetaMater` try using our vignette with example files provided with this package. 


