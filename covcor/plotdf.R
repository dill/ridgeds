## Simulation -- correlated covariates in the detection function

# in this situation we have 1 covariate that's actually useful, the other
# is correlated with the first

# generic simulation stuff

library(Distance)
library(mmds)

options(stingsAsFactors=FALSE)
source("../generic_sim_fns/sim.hn.R")
source("../generic_sim_fns/make.data.R")

## options
# number of simulations to do per sample size
n.sims <- 200
# number of distances observed
n.samps <- c(30,60,120,240,480)
# correlations
corrs <- c(0.7,0.8,0.9)
# truncation
width <- 1
# formula used to *generate* the data
formula <- ~cov1

### BETA distributed covariates
### simulation setting 1
b <- 1.1
a <- 4
pars <- c(0.1,-0.2)

### beta testing
betamean <- a/(a+b)
betavar <- (a*b)/((a+b)^2 * (a+b+1))
### testing parameters
corr <- 0.8
sample.size <- 1000
#covdata <- data.frame(cov1=rnorm(sample.size))
covdata <- data.frame(cov1=rbeta(sample.size,a,b))
covdata$cov2 <- corr*covdata$cov1 + sqrt(1-corr^2)*rnorm(sample.size,betamean,sqrt(betavar))


### T distributed covariates
### semi-stolen from Marques et al (2007)
#library(tmvtnorm)
#pars <- c(log(0.9),-0.8)
#corr <- 0.8
#sample.size <- 1000
#covdata <- data.frame(cov1=rtmvt(sample.size,mean=0,lower=0,upper=360,sigma=8000))
###hist(covdata$cov1)
### create pseduo hours
#covdata$cov2 <- round((covdata$cov1/360)*4,0)
### standardise
#covdata$cov1 <- covdata$cov1/sd(covdata$cov1)

# make the distance data frame
this.data <- make.data(sample.size,width,pars,covdata,formula)

#### plot that
#par(mfrow=c(2,3))
#plot(this.data$distance,this.data$cov1,xlab="Distance",ylab="cov1")
#plot(this.data$distance,this.data$cov2,xlab="Distance",ylab="cov2")
#plot(this.data$cov1,this.data$cov2,xlab="cov1",ylab="cov2")
#hist(this.data$distance,xlab="Distance")
#hist(this.data$cov1,xlab="cov1")
#hist(this.data$cov2,xlab="cov2")
##
##
#mod1 <- try(ds(this.data,truncation=width,order=0))
#mod1.cov <- try(ds(this.data,truncation=width,formula=~cov1,order=0))
#mod12.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2,order=0))
#par(mfrow=c(1,3))
#plot(mod1)
#plot(mod1.cov)
#plot(mod12.cov)
#
#dev.new()

fake.object <- list()
fake.object$width <- width
fake.object$data <- this.data
fake.object$pars <- pars
fake.object$mix.terms <- 1
fake.object$zdim <- 1
fake.object$z <- list(as.matrix(covdata))
fake.object$pt <- FALSE
fake.object$model.formula <- formula

mmds:::plot.ds.mixture(fake.object)

