# generic simulation stuff

library(Distance)
library(doMC)
library(foreach)

options(stingsAsFactors=FALSE)
source("sim.hn.R")

options(cores=2)
registerDoMC()

# function to generate data
make.data <- function(sample.size,width,pars,covdata=NULL,formula=~1){

  if(formula != "~1"){
    z <- model.matrix(formula,covdata)
  }else{
    z <- rep(1,sample.size)
  }

  # simulate some data
  ddat <- sim.hn(sample.size,z, pars, width)

  # build a data.frame
  if(is.null(covdata)){
    dat <- data.frame(distance=ddat)
  }else{
    dat <- cbind(data.frame(distance=ddat),
                 as.data.frame(covdata))
  }
  return(dat)
}

# shortcut function to store results
store_results <- function(mod, model.name, sample.size, sim, res){
  # mod         - fitted model object (call was wrapped in a try())
  # model.name  - string to identify model
  # sample.size - number of observations
  # sim         - simulation id
  # results     - data.frame to append to
  if(class(mod)!="try-error" | !is.null(mod$model)){
    p <- sample.size/mod$ddf$Nhat
  }else{
    p <- NA
  }
  res <- rbind(res,
               c(sim, sample.size, model.name, sample.size/mod$ddf$Nhat))
  return(res)
}

## options
# number of simulations to do per sample size
n.sims <- 200
# number of distances observed
n.samps <- c(30,60,120,240,480)
# truncation
width <- 0.75
#parameters
pars <- log(0.3)
# include the covariate name in the formula to have it be non-spurious
formula <- ~1


#plot(this.data$distance,this.data$cov1,xlab="Distance")

big.res <- c()

for(sample.size in n.samps){
  # loop over sims
  results <- foreach(sim = 1:n.sims, .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {

    res <- c()

    ## generate data
    # what is the covariate we're including?
    covdata <- data.frame(cov1=runif(sample.size,-1,1))
    # make the distance data frame
    this.data <- make.data(sample.size,width,pars,covdata,formula)

    # fit the "true" model -- model data was generated from
    mod1 <- try(ds(this.data,truncation=width))
    res <- store_results(mod1, "hn", sample.size, sim, res)

    ## fit with "wrong" models...

    # hn with covariate
    mod1.cov <- try(ds(this.data,truncation=width,formula=~cov1))
    res <- store_results(mod1.cov, "hn+cov1", sample.size, sim, res)

    # hazard rate
    mod1.hr <- try(ds(this.data,truncation=width,key="hr"))
    res <- store_results(mod1.hr, "hr", sample.size, sim, res)

    # hazard rate with covariate
    mod1.hr.cov <- try(ds(this.data,truncation=width,formula=~cov1,key="hr"))
    res <- store_results(mod1.hr.cov, "hr+cov1", sample.size, sim, res)

    return(res)
  }

  big.res <- rbind(big.res,results)
}

write.csv(big.res, file="junk.csv")


## quick plot
# load data here 
results <- as.data.frame(big.res)
names(results) <- c("sim","n","model","p")
results$p <- as.numeric(as.character(results$p))
results$n <- as.factor(as.numeric(as.character(results$n)))
library(ggplot2)
p <- ggplot(results)
p <- p + geom_boxplot(aes(x=n,y=p))
p <- p + facet_wrap(~model)
p <- p + scale_y_continuous(limits=c(0,1))
print(p)




