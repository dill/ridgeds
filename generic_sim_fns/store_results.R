# shortcut function to store results
store_results <- function(mod, model.name, sample.size, sim, res, corr=0){
  # mod         - fitted model object (call was wrapped in a try())
  # model.name  - string to identify model
  # sample.size - number of observations
  # sim         - simulation id
  # results     - data.frame to append to

  # dummy
  p <- NA
  varp <- NA
  cov1 <- NA
  cov2 <- NA
  b0 <- NA

  # other variables
  other.vars <- c(sim, sample.size, model.name, corr)

  # if fitting went okay, then store some values
  if(!is.null(mod)){
    if(class(mod)!="try-error"){
      if(!is.null(mod$ddf)){
        p <- sample.size/mod$ddf$Nhat
        varp <- summary(mod$ddf)$average.p.se

        # pull out parameters
        b0 <- mod$ddf$par[1]
        if(!is.null(mod$ddf$par["cov1"])){
          cov1 <- mod$ddf$par["cov1"]
        }
        if(!is.null(mod$ddf$par["cov2"])){
          cov2 <- mod$ddf$par["cov2"]
        }
      }
    }
  }
  res <- rbind(res,
               c(other.vars, "p", p),
               c(other.vars, "b0", b0),
               c(other.vars, "cov1", cov1),
               c(other.vars, "cov2", cov2),
               c(other.vars, "varp", varp))
  return(res)
}
