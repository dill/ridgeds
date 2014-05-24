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

options(cores=8)
registerDoMC()

## options
# number of simulations to do per sample size
n.sims <- 2#00
# number of distances observed
n.samps <- c(30,60,120,240,480)
# correlations
corrs <- c(0.7,0.8,0.9)
# truncation
width <- 1
#parameters
pars <- c(log(0.75),1.2)
# include the covariate name in the formula to have it be non-spurious
formula <- ~cov1

#### testing parameters
#corr <- 0.8
#sample.size <- 100
##covdata <- data.frame(cov1=rnorm(sample.size))
a<-0.45;b<-0.45
betamean <- a/(a+b)
betavar <- (a*b)/((a+b)^2 * (a+b+1))
#covdata <- data.frame(cov1=2*(rbeta(sample.size,a,b)-1/2))
#covdata$cov2 <- corr*covdata$cov1 + sqrt(1-corr^2)*rnorm(sample.size,betamean,sqrt(betavar))

# make the distance data frame
#this.data <- make.data(sample.size,width,pars,covdata,formula)

### plot that
#par(mfrow=c(2,3))
#plot(this.data$distance,this.data$cov1,xlab="Distance",ylab="cov1")
#plot(this.data$distance,this.data$cov2,xlab="Distance",ylab="cov2")
#plot(this.data$cov1,this.data$cov2,xlab="cov1",ylab="cov2")
#hist(this.data$distance,xlab="Distance")
#hist(this.data$cov1,xlab="cov1")
#hist(this.data$cov2,xlab="cov2")

#mod1 <- try(ds(this.data,truncation=width,order=0))
#mod1.cov <- try(ds(this.data,truncation=width,formula=~cov1,order=0))
#mod12.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2,order=0))
#par(mfrow=c(1,3))
#plot(mod1)
#plot(mod1.cov)
#plot(mod12.cov)

#####################

big.res <- c()

# build the simulation settings
setting.grid <- expand.grid(corr=corrs,
                            sample.size=n.samps)


for(i in 1:nrow(setting.grid)){

  # grab this set of simulation settings
  these.settings <- setting.grid[i,]
  corr <- these.settings$corr
  sample.size <- these.settings$sample.size

  # loop over sims
  results <- foreach(sim = 1:n.sims, .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {

    res <- c()

    ## generate data
    covdata <- data.frame(cov1=2*(rbeta(sample.size,a,b)-1/2))
    covdata$cov2 <- corr*covdata$cov1 +
                      sqrt(1-corr^2)*rnorm(sample.size,betamean,sqrt(betavar))
    # make the distance data frame
    this.data <- make.data(sample.size,width,pars,covdata,formula)

    # fit the "true" model -- model data was generated from
    mod1 <- try(ds(this.data,truncation=width,order=0))
    res <- store_results(mod1, "hn", sample.size, sim, res, corr=corr)

    ## fit with "wrong" models...

    # hn with covariate 1
    mod1.cov <- try(ds(this.data,truncation=width,formula=~cov1,order=0))
    res <- store_results(mod1.cov, "hn+cov1", sample.size, sim, res, corr=corr)

    # hn with covariate 1 + 2
    mod12.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2,order=0))
    res <- store_results(mod12.cov, "hn+cov1+cov2", sample.size, sim, res, corr=corr)

#    # hazard rate
#    mod1.hr <- try(ds(this.data,truncation=width,key="hr",order=0))
#    res <- store_results(mod1.hr, "hr", sample.size, sim, res, corr=corr)
#
#    # hazard rate with covariate 1
#    mod1.hr.cov <- try(ds(this.data,truncation=width,formula=~cov1,key="hr",order=0))
#    res <- store_results(mod1.hr.cov, "hr+cov1", sample.size, sim, res, corr=corr)
#
#    # hazard rate with covariate 1 + 2
#    mod12.hr.cov <- try(ds(this.data,truncation=width,formula=~cov1+cov2,
#                           key="hr",order=0))
#    res <- store_results(mod12.hr.cov, "hr+cov1+cov2", sample.size, sim, res, corr=corr)

    return(res)
  }

  big.res <- rbind(big.res,results)
}

write.csv(big.res, file="covcor.csv")


