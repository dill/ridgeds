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
