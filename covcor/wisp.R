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
n.sims <- 500
# correlations
corrs <- c(0.7,0.8,0.9,0.95,0.99)
# formula used to *generate* the data
formula <- ~cov1
# truncation
width <- 1
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
                                   exposure.method="beta", exposure.min=0,
                                   exposure.max=10, exposure.mean=2,
                                   exposure.shape=2)
  mypop <- generate.population(mypop.pars, seed=t.seed)
  # set up the design
  mydes.pars <- setpars.design.lt(myreg, n.transects=20, n.units=20,
                                  visual.range=1.1)
  mydes <- generate.design.lt(mydes.pars, seed=t.seed)

# first sim was 0.1, 0.6
# second sim was 0.1, 1.2
# second sim was 0.7, 1.2
  mysurvey.pars <- setpars.survey.lt(mypop, mydes, disthalf.min=0.15,
                                     disthalf.max=0.8, model = "half.normal")

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
      dat$cov3 <- corr*dat$cov1 + scaling*sqrt(1-corr^2)*rnorm(nrow(dat),betamean,sqrt(betavar))

# normalise the covariates
dat[,2:4] <- (dat[,2:4]-colMeans(dat[,2:4]))/apply(dat[,2:4],2,sd)

      # fit models
      mm1 <- try(ds(dat, truncation=width,adjustment=NULL))
      mm2 <- try(ds(dat, truncation=width,formula=~cov1,adjustment=NULL))
      mm3 <- try(ds(dat, truncation=width,formula=~cov1+cov2,adjustment=NULL))
      mm5 <- try(ds(dat, truncation=width,formula=~cov1+cov2+cov3,adjustment=NULL))

# now do some PCA
dat2 <- as.matrix(dat[,2:4])
pc <- eigen((t(dat2)%*%dat2))
pc.dat <- dat2%*%pc$vectors
pc.dat <- as.data.frame(pc.dat)
names(pc.dat) <- c("cov1","cov2","cov3")
pc.dat$distance <- dat$distance

# fit PCA models
mm2.pc <- try(ds(pc.dat, truncation=width,formula=~cov1,adjustment=NULL))
mm3.pc <- try(ds(pc.dat, truncation=width,formula=~cov1+cov2,adjustment=NULL))
mm5.pc <- try(ds(pc.dat, truncation=width,formula=~cov1+cov2+cov3,adjustment=NULL))


Ncov <- sum(!is.na(mysamp$detected))

      # extract the results
      res <- store_results(mm1, "hn", this.n.pop, sim, res, corr=corr, Ncov)
      res <- store_results(mm2, "hn+cov1", this.n.pop, sim, res, corr=corr, Ncov)
      res <- store_results(mm3, "hn+cov1+cov2", this.n.pop, sim, res,corr=corr, Ncov)
      res <- store_results(mm5, "hn+cov1+cov2+cov3", this.n.pop, sim, res,corr=corr, Ncov)

# save pca results
      res <- store_results(mm2.pc, "PCAhn+cov1", this.n.pop, sim, res, corr=corr, Ncov)
      res <- store_results(mm3.pc, "PCAhn+cov1+cov2", this.n.pop, sim, res,corr=corr, Ncov)
      res <- store_results(mm5.pc, "PCAhn+cov1+cov2+cov3", this.n.pop, sim, res,corr=corr, Ncov)



      return(res)

    } # end sim i

    big.res <- rbind(big.res,results)

  } # end loop over correlations

} # end loop over population sizes

write.csv(big.res, file="covcor-wisp-smallp-PCA.csv")


