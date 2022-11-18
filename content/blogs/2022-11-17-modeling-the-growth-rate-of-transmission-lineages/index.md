---
title: Bayesian change-point modeling using Nimble
author: Magnus Nygård Osnes
date: '2022-11-17'
slug: []
categories:
  - TLs
  - LineageHomology
  - Changepoint modeling
tags: []
image: ~
description: ~
toc: ~
---

## Modeling the growth rate of transmission lineages

We are working on a project where we try to illuminate some aspect of lineage dynamics of *Neisseria gonorrhoeae*. We are focusing on four locations divided into small scale: Norway and Australia, and large scale: Europe and the USA. For each of the four geographical locations we have used [LineageHomology](https://github.com/magnusnosnes/LineageHomology) to define transmission lineages wrt. each of those locations. A “transmission lineage” is group of tips on a phylogeny that are are estimated to stem from the same importation event. Hence such groups may be used to retrospectively assess how important a single importation events in terms of onward transmission of the disease. \## **Research question** In this post I look at the relation between the growth rate of the lineages and the observed gender distribution of the infected hosts in the transmission lineages. The reason we wanted to investigate this question is that the the observed gender distribution of the hosts likely reflects the sexual network structure (such as heterosexual communities, men who have sex with men or mixed networks). We wanted to see if there were differences in the growth rate of the lineages depending on the networks they “inhabit”.

## A look at the data

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-1.png" width="672" />

Looking at the data we thought that there may be some threshold in the observed fraction of men that are indicative of sexual network structure and hence could imply a different growth rate. Visual inspection of the data suggests that above 70% men there tends to be higher growth rates. This may indicate that we have a change-point around 0.7 like this:

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="672" />

Instead of defining this threshold ourselves from visual inspection, we wanted to formally estimate it from using the data, along with the and expected value of the growth rate above and below the threshhold.

<!--  ### Size of the transmission lineages


 
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="768" />
-->

## Bayesian change-point modeling

To describe the “threshold” / “change-point” model we use the r-package [NIMBLE](https://r-nimble.org/) for statistical programming.

## Change-point model

### Defining data for NIMBLE

Our data looks like:

``` r
library(nimble)

head(small_scale[,c("GrowthRate","sexDistribution", "locations")])
##    GrowthRate sexDistribution locations
## 1 0.161981919       0.9944134    Norway
## 2 0.370304650       0.9666667    Norway
## 3 1.070982756       0.5795455    Norway
## 4 0.350432551       0.9759036    Norway
## 5 2.829087867       0.9565217    Norway
## 6 0.009947276       0.8596491    Norway
summary(small_scale[,c("GrowthRate","sexDistribution", "locations")])
##    GrowthRate       sexDistribution      locations 
##  Min.   :-0.06797   Min.   :0.2857   Australia:30  
##  1st Qu.: 0.17727   1st Qu.:0.7500   Norway   :20  
##  Median : 0.60928   Median :0.8750                 
##  Mean   : 0.93092   Mean   :0.8350                 
##  3rd Qu.: 1.38518   3rd Qu.:0.9286                 
##  Max.   : 4.50239   Max.   :1.0000
```

NIMBLE requires the data and explanatory variables to in a named lists. Additionally we need to supply other constants that will be used in the model definition such as the length of the data “n”.

``` r
data <- list(y=small_scale[,"GrowthRate"]) #The growthrate that we are trying to explain.
n = length(small_scale[,"GrowthRate"]) #Number of data points
constants = list(n=n, genderDistribution=small_scale$sexDistribution, #Constants to be used in the Nimblemodel
                 location=as.numeric(as.factor(small_scale$locations)))
```

### Defining the change-point model in NIMBLE

We specify an threshold parameter “alpha” as a continuous value on the scale \[0,1\] fraction of men. We assume that the mean growth rate in is different before and after, but similar for all values. We assume the the errors follow a normal distribution around a mean. For the change-point we assume a prior uniform distribution on the interval \[0.5, 1\].

``` r
changepoint <- nimbleCode({
  
  #Prior distributions 
  alpha ~ dunif(0.5,1) #Changepoint
  sigma1 ~ dunif(0, 5) # prior for variance components based on Gelman (2006)
  for(k in 1:2){
    mu1[k]~ dnorm(0, sd = 10) # Expected value before the change point
    mu2[k]~ dnorm(0, sd = 10) # Expected value before the after the change point
  }
  
  #Model
  for(i in 1:n) {
    turn_on[i] <- step(genderDistribution[i]-alpha) #0 if y[i] is smaller than and 1 otherwise.
    mug[i] <- mu1[location[i]]*turn_on[i]+mu2[location[i]]*(1-turn_on[i]) #Switch means depending on if we are above or below the changepoint.
    
    #Likelihood of the data (growthrate)
    y[i] ~ dnorm(mug[i],sigma1)
  }
})
```

### Setting initial values

Before running the model we need to supply initial values to all the parameters to initialize the [MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo).

``` r
inits <- list(mu1 =c(1,1),mu2=c(1,1),
              sigma1 = 1,
              #sigma2 = 1,
              alpha=0.8)
```

### Running the change-point model

We run the model on the small scale data from Australia and Norway ” and summarize the results on

``` r
#Run the MCMC
samples_small_scale = nimbleMCMC(code = changepoint, constants = constants,
                       data = data, inits = inits,
                       nchains = 2, niter = 1000000,
                       summary = TRUE, WAIC = TRUE,thin = 100,
                       monitors = c('mu1','mu2','sigma1',
                                    #'sigma2',
                                    'alpha'))
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
##   [Warning] There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.
```

Next, we summarize results using the [MCMCvis](https://cran.r-project.org/web/packages/MCMCvis/index.html) r-package.

### Australia and Norway

``` r
library(MCMCvis)
MCMCsummary(object = samples_small_scale$samples)
##             mean        sd         2.5%       50%      97.5% Rhat n.eff
## alpha  0.6677818 0.1791179   0.50327949 0.5663644  0.9959115    1  6684
## mu1[1] 0.9824472 3.8958175  -9.59835442 1.1865059  9.5462415    1 20000
## mu1[2] 0.5530838 0.2959798  -0.03316428 0.5512741  1.1554907    1 20000
## mu2[1] 0.5716374 7.8573234 -17.32340834 1.2563820 17.2369773    1 20115
## mu2[2] 0.2767973 0.5462392  -0.89603444 0.3190279  1.3150664    1 19152
## sigma1 1.1223159 0.2313403   0.71645324 1.1062569  1.6264014    1 20325
MCMCplot(object = samples_small_scale$samples,params=c("alpha"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="672" />

We ensure that the chains are mixing properly by inspecting the effective sample size (last column)

``` r

MCMCplot(object = samples_small_scale$samples,params=c("mu1", "mu2"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-1.png" width="672" />

``` r
#MCMCtrace(object = samples_small_scale$samples,params=c("alpha")) #Look at the traceplots
```

mu1\[1\] and mu2\[1\] is the expected value above and below the change-point Australia, and mu1\[2\], mu1\[2\] for Norway. From the width

### Europe and the USA

``` r
MCMCsummary(object = samples_large_scale$samples)
##             mean         sd        2.5%        50%     97.5% Rhat n.eff
## alpha  0.9314393 0.08110011  0.71532038 0.96714021 0.9987211    1 10315
## mu1[1] 0.1338956 0.28843610 -0.51585054 0.14071660 0.7336394    1 10000
## mu1[2] 0.5182609 0.35910112 -0.03367521 0.49459113 1.2640045    1 10476
## mu2[1] 0.1999654 0.13375860 -0.05993548 0.19689025 0.4765697    1 10213
## mu2[2] 0.1017999 0.10303781 -0.09277714 0.09853546 0.3138090    1 10000
## sigma1 4.7388575 0.24570241  4.09715258 4.81082907 4.9929444    1 10000
MCMCplot(object = samples_large_scale$samples,params=c("alpha"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" width="672" />

``` r
MCMCplot(object = samples_large_scale$samples,params=c("mu1", "mu2"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-2.png" width="672" />

``` r
#MCMCtrace(object = samples_large_scale$samples,params=c("alpha")) #Look at the traceplots
```

## Regression model

We use the same dataset but change to a regression model.

### Defining the model

``` r
regression <- nimbleCode({
  
  sigma1 ~ dunif(0, 5)        # prior for variance components based on Gelman (2006)
  
  #Priors for regression slope
  for(k in 1:2){
    mu1[k]~ dnorm(0, sd = 10) 
    beta[k] ~ dnorm(0, sd = 10) #Regression slope
  }
  
  
  
  for(i in 1:n) {
    mug[i] <- mu1[location[i]]+beta[location[i]]*genderDistribution[i] #Define the expected value
    
    #Likelihood of the growth rate
    y[i] ~ dnorm(mug[i],sigma1)
  }
})
```

### Norway and Australia

``` r
data <- list(y=small_scale[,1])
constants = list(n=n, genderDistribution=small_scale$sexDistribution,
                 location=as.numeric(as.factor(small_scale$locations)))
inits <- list(mu1 =c(1,1),
              sigma1 = 1,
              #sigma2 = 1,
              beta=c(0.8,0.8))

regression_samples_small_scale = nimbleMCMC(code = regression, constants = constants,
                       data = data, inits = inits,
                       nchains = 2, niter = 1000000,nburnin = 1000000/2,
                       summary = TRUE, WAIC = TRUE,thin = 100,
                       monitors = c('mu1',
                                    'sigma1',
                                    'beta'))
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
##   [Warning] There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.
par(mfrow=c(1,1))
MCMCsummary(object = regression_samples_small_scale$samples)
##                mean        sd       2.5%        50%    97.5% Rhat n.eff
## beta[1] -2.21144862 2.0655486 -6.2728213 -2.2330787 1.808282    1  4693
## beta[2]  0.56143496 1.0625851 -1.5663978  0.5804113 2.651546    1 10299
## mu1[1]   3.09452913 1.7571339 -0.3373730  3.1101156 6.538150    1  4689
## mu1[2]   0.03882466 0.8947434 -1.7045187  0.0328491 1.824354    1 10337
## sigma1   1.12868494 0.2313801  0.7213114  1.1101052 1.621900    1 11386
MCMCplot(object = regression_samples_small_scale$samples,params=c("mu1","beta"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" width="672" />

``` r
#MCMCtrace(object = regression_samples_small_scale$samples,params=c("beta"))
```

### Europe and the USA

``` r
data <- list(y=large_scale[,1])
constants = list(n=n, genderDistribution=large_scale$sexDistribution,
                 location=as.numeric(as.factor(large_scale$locations)))
inits <- list(mu1 =c(1,1),
              sigma1 = 1,
              #sigma2 = 1,
              beta=c(0.8,0.8))

regression_samples_large_scale = nimbleMCMC(code = regression, constants = constants,
                       data = data, inits = inits,
                       nchains = 2, niter = 1000000,nburnin = 1000000/2,
                       summary = TRUE, WAIC = TRUE,thin = 100,
                       monitors = c('mu1',
                                    'sigma1',
                                    'beta'))
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
##   [Warning] There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.
par(mfrow=c(1,1))
MCMCsummary(object = regression_samples_large_scale$samples)
##                mean        sd       2.5%         50%     97.5% Rhat n.eff
## beta[1] -0.18467150 0.6040629 -1.3913003 -0.17793571 0.9790354    1  9445
## beta[2]  0.03587699 0.5229691 -0.9904779  0.03742878 1.0572973    1 10003
## mu1[1]   0.33209881 0.5029887 -0.6484135  0.32995573 1.3335549    1  9415
## mu1[2]   0.10457476 0.4209626 -0.7157290  0.10495709 0.9319498    1 10040
## sigma1   4.72455469 0.2564046  4.0559068  4.79786131 4.9915727    1 10000
MCMCplot(object = regression_samples_large_scale$samples,params=c("mu1","beta"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" width="672" />

``` r
#MCMCtrace(object = regression_samples_large_scale$samples,params=c("beta"))
```

## Summing up

Although visually inspection indicated a change-point in the data. We are unable to estimate such a point with a narrow credibility interval based on the data. The regression model find any effect of sex-distribution (regression slope 97.5% credibility interval covered zero for all locations. Hence we remain uncertain that a change point exists and if the observed gender distribution in a transmission lineage can be used to determine its growth rate.
