library(wisp)
library(Distance)
library(doMC)
library(foreach)

options(cores=8)
registerDoMC()

options(stingsAsFactors=FALSE)
source("../generic_sim_fns/store_results.R")

set.seed(1234)

t.seed <- NULL


## options
# number of simulations to do per sample size
n.sims <- 200
# correlations
corrs <- c(0.7,0.8,0.9)
# formula used to *generate* the data
formula <- ~cov1
# truncation
width <- 0.9
# population sizes
n.pops <- c(1000,5000)


## generate region and density
myreg <- generate.region(x.length=100, y.width=100)
mydens <- generate.density(myreg, southwest=10,southeast=10,
                           northwest=10, nint.x=10, nint.y=10)

big.res <- c()

for(this.n.pop in n.pops){
  # population
  mypop.pars <- setpars.population(mydens, number.groups=this.n.pop,
                                   size.method="user", size.prob=1, size.values=1,
                                   exposure.method="beta", exposure.min=2,
                                   exposure.max=10, exposure.mean=6,
                                   exposure.shape=1)
  mypop <- generate.population(mypop.pars, seed=t.seed)
  # set up the design
  mydes.pars <- setpars.design.lt(myreg, n.transects=10, n.units=10, visual.range=3.5)
  mydes <- generate.design.lt(mydes.pars, seed=t.seed)
  mysurvey.pars <- setpars.survey.lt(mypop, mydes, disthalf.min=0.3,
                                     disthalf.max=0.6, model = "half.normal")

  # exposure parameters, used below to generate correlated data
  a <- mypop.pars$exposure.alpha
  b <- mypop.pars$exposure.beta
  betamean <- a/(a+b)
  betavar <- (a*b)/((a+b)^2 * (a+b+1))
  scaling <- mypop.pars$exposure.max-mypop.pars$exposure.min

  # loop over correlations
  for(corr in corrs){

    ## simulation
    results <- foreach(sim = 1:n.sims, .combine=rbind,
                       .inorder=FALSE, .init=c()) %dopar% {

      # storage
      res <- c()

      # generate data
      mysamp <- generate.sample.lt(mysurvey.pars, seed=t.seed)
      samp.ind <- mysamp$detected==1 & !is.na(mysamp$detected)
      dat <- data.frame(distance = mysamp$distance[samp.ind],
                        cov1     = mysamp$population$exposure[samp.ind])

      dat$cov2 <- corr*dat$cov1 + scaling*sqrt(1-corr^2)*rnorm(nrow(dat),betamean,sqrt(betavar))


      # fit models
      mm1 <- try(ds(dat, truncation=width,adjustment=NULL))
      mm2 <- try(ds(dat, truncation=width,formula=~cov1,adjustment=NULL))
      mm3 <- try(ds(dat, truncation=width,formula=~cov1+cov2,adjustment=NULL))


      # extract the results
      res <- store_results(mm1, "hn", this.n.pop, sim, res, corr=corr)
      res <- store_results(mm2, "hn+cov1", this.n.pop, sim, res, corr=corr)
      res <- store_results(mm3, "hn+cov1+cov2", this.n.pop, sim, res,corr=corr)

      return(res)

    } # end sim i

    big.res <- rbind(big.res,results)

  } # end loop over correlations

} # end loop over population sizes

write.csv(big.res, file="covcor-wisp.csv")


