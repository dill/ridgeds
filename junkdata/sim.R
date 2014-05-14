# generic simulation stuff

library(Distance)
library(doMC)
library(foreach)

options(stingsAsFactors=FALSE)
source("sim.hn.R")

options(cores=2)
registerDoMC()



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

## options
# number of simulations to do per sample size
n.sims <- 1000
# number of distances observed
n.samps <- c(30,60,120,240,480)
# truncation
width <- 0.75
#parameters
pars <- log(0.3)
# what is the covariate we're including?
covdata <- data.frame(cov1=runif(sample.size,-1,1))
# include the covariate name in the formula to have it be non-spurious
formula <- ~1


#plot(this.data$distance,this.data$cov1,xlab="Distance")

big.res <- c()

for(sample.size in n.samps){
  # loop over sims
  results <- foreach(sim = 1:n.sims, .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {

    # generate data
    this.data <- make.data(sample.size,width,pars,covdata,formula)

    # fit the "true" model -- model data was generated from
    mod1 <- try(ds(this.data,truncation=width))

    if(class(mod1)!="try-error"){
      res <- c(sim, sample.size, "hn", sample.size/mod1$ddf$Nhat)
    }else{
      res <- c(sim, sample.size, "hn", NA)
    }

    # fit with "wrong" models...

    # hn with covariate
    mod1.cov <- ds(this.data,truncation=width,formula=~cov1)
    if(class(mod1.cov)!="try-error"){
      res <- rbind(res,
                   c(sim, "hn+cov1", sample.size/mod1.cov$ddf$Nhat))
    }else{
      res <- rbind(res,
                   c(sim, "hn+cov1", NA))
    }

    # hazard rate
    mod1.hr <- ds(this.data,truncation=width,key="hr")
    if(class(mod1.hr)!="try-error"){
      res <- rbind(res,
                   c(sim, "hr", sample.size/mod1.hr$ddf$Nhat))
    }else{
      res <- rbind(res,
                   c(sim, "hr", NA))
    }

    # hazard rate with covariate
    mod1.hr.cov <- ds(this.data,truncation=width,formula=~cov1,key="hr")
    if(class(mod1.hr.cov)!="try-error"){
      res <- rbind(res,
                   c(sim, "hr+cov1", sample.size/mod1.hr.cov$ddf$Nhat))
    }else{
      res <- rbind(res,
                   c(sim, "hr+cov1", NA))
    }

    return(res)

  }

  big.res <- rbind(big.res,results)
}

write.csv(big.res, file="junk.csv")


## quick plot
#results <- as.data.frame(results)
#names(results) <- c("sim","model","p")
#library(ggplot2)
#p <- ggplot(results)
#p <- p + geom_boxplot(aes(x=model,y=p))
#print(p)




