## Simulation -- correlated covariates in the detection function

# in this situation we have 1 covariate that's actually useful, the other
# is correlated with the first

# generic simulation stuff

library(Distance)
library(doMC)
library(foreach)

options(stingsAsFactors=FALSE)
source("../generic_sim_fns/sim.hn.R")
source("../generic_sim_fns/make.data.R")
source("../generic_sim_fns/store_results.R")

options(cores=2)
registerDoMC()

## options
# number of simulations to do per sample size
n.sims <- 200
# number of distances observed
n.samps <- c(30,60,120,240,480)
# truncation
width <- 0.75
#parameters
pars <- c(log(0.3),log(1.5))
# include the covariate name in the formula to have it be non-spurious
formula <- ~cov1

## test plots
corr <- 0.8
#sample.size <- 1000
#covdata <- data.frame(cov1=rnorm(sample.size))
#covdata$cov2 <- corr*covdata$cov1 + sqrt(1-corr^2)*rnorm(sample.size)
# make the distance data frame
#this.data <- make.data(sample.size,width,pars,covdata,formula)

#par(mfrow=c(2,3))
#plot(this.data$distance,this.data$cov1,xlab="Distance",ylab="cov1")
#plot(this.data$distance,this.data$cov2,xlab="Distance",ylab="cov2")
#plot(this.data$cov1,this.data$cov2,xlab="cov1",ylab="cov2")
#hist(this.data$distance,xlab="Distance")
#hist(this.data$cov1,xlab="cov1")
#hist(this.data$cov2,xlab="cov2")

big.res <- c()

for(sample.size in n.samps){
  # loop over sims
  results <- foreach(sim = 1:n.sims, .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {

    res <- c()

    ## generate data
    covdata <- data.frame(cov1=rnorm(sample.size))
    covdata$cov2 <- corr*covdata$cov1 + sqrt(1-corr^2)*rnorm(sample.size)
    # make the distance data frame
    this.data <- make.data(sample.size,width,pars,covdata,formula)

    # fit the "true" model -- model data was generated from
    mod1 <- try(ds(this.data,truncation=width))
    res <- store_results(mod1, "hn", sample.size, sim, res)

    ## fit with "wrong" models...

    # hn with covariate 1
    mod1.cov <- try(ds(this.data,truncation=width,formula=~cov1))
    res <- store_results(mod1.cov, "hn+cov1", sample.size, sim, res)

    # hn with covariate 1 + 2
    mod12.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2))
    res <- store_results(mod12.cov, "hn+cov1+cov2", sample.size, sim, res)

    # hazard rate
    mod1.hr <- try(ds(this.data,truncation=width,key="hr"))
    res <- store_results(mod1.hr, "hr", sample.size, sim, res)

    # hazard rate with covariate 1
    mod1.hr.cov <- try(ds(this.data,truncation=width,formula=~cov1,key="hr"))
    res <- store_results(mod1.hr.cov, "hr+cov1", sample.size, sim, res)

    # hazard rate with covariate 1 + 2
    mod12.hr.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2,key="hr"))
    res <- store_results(mod12.hr.cov, "hr+cov1+cov2", sample.size, sim, res)

    return(res)
  }

  big.res <- rbind(big.res,results)
}

write.csv(big.res, file="covcor.csv")





