library(wisp)
library(Distance)
library(doMC)
library(foreach)

options(cores=8)
registerDoMC()

options(stingsAsFactors=FALSE)
source("../generic_sim_fns/store_results.R")

#set.seed(1234)

t.seed <- NULL


## options
# number of simulations to do per sample size
n.sims <- 1#400
# correlations
corrs <- c(0.8,0.9,0.95)
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

this.n.pop <- n.pops[1]
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

mysurvey.pars <- setpars.survey.lt(mypop, mydes, disthalf.min=0.15,
                                   disthalf.max=0.8, model = "half.normal")

# exposure parameters, used below to generate correlated data
a <- mypop.pars$exposure.alpha
b <- mypop.pars$exposure.beta
betamean <- a/(a+b)
betavar <- (a*b)/((a+b)^2 * (a+b+1))
scaling <- mypop.pars$exposure.max-mypop.pars$exposure.min

# loop over correlations
corr <- corrs[1]

cov1 <- mypop$exposure
cov2 <- corr*cov1 + scaling*sqrt(1-corr^2)*rnorm(this.n.pop,
                                                 betamean,sqrt(betavar))
cov3 <- corr*cov1 + scaling*sqrt(1-corr^2)*rnorm(this.n.pop,
                                                 betamean,sqrt(betavar))

# generate data
mysamp <- generate.sample.lt(mysurvey.pars, seed=t.seed)
samp.ind <- mysamp$detected==1 & !is.na(mysamp$detected)
dat <- data.frame(distance = mysamp$distance[samp.ind],
                  cov1     = mysamp$population$exposure[samp.ind],
                  cov2     = cov2[samp.ind],
                  cov3     = cov3[samp.ind])

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


par(mfrow=c(2,3))
plot(mm2)
plot(mm3)
plot(mm5)
plot(mm2.pc)
plot(mm3.pc)
plot(mm5.pc)

