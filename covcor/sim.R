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
a<-0.45;b<-0.45
pars <- c(1,-1.2)
betamean <- a/(a+b)
betavar <- (a*b)/((a+b)^2 * (a+b+1))

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
    # BETA
    covdata <- data.frame(cov1=rbeta(sample.size,a,b))
    covdata$cov2 <- corr*covdata$cov1 + sqrt(1-corr^2)*rnorm(sample.size,betamean,sqrt(betavar))

    ## non central T
    #covdata <- data.frame(cov1=rtmvt(sample.size,mean=100,
    #                                 lower=0,upper=360,sigma=8000))
    #covdata$cov2 <- round((covdata$cov1/360)*4,0)
    ## standardise
    #covdata$cov1 <- covdata$cov1/sd(covdata$cov1)

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

write.csv(big.res, file=paste0("covcor-beta-",
                               paste(pars,collapse="-"),".csv"))


