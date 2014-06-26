# shortcut function to store results
store_results <- function(mod, pop.size, sim, corr=0, Ncov){
  # mod         - fitted model object (call was wrapped in a try())
  # model.name  - string to identify model
  # pop.size    - population size
  # sim         - simulation id

  # dummy
  p <- NA
  varp <- NA
  cov1 <- NA
  cov2 <- NA
  b0 <- NA
  aic <- NA
  n <- NA

  # other variables
  other.vars <- c(sim, pop.size, mod$model.name, corr,Ncov)

  # if fitting went okay, then store some values
  if(!is.null(mod)){
    if(class(mod)!="try-error"){
      if(!is.null(mod$ddf)){

        # average p is n/Nhat
        p <- summary(mod$ddf)$average.p

        # variance of average p and Nhat
        varp <- summary(mod$ddf)$average.p.se
        varN <- summary(mod$ddf)$Nhat.se

        # pull out parameters
        b0 <- mod$ddf$par[1]
        if(!is.null(mod$ddf$par["cov1"])){
          cov1 <- mod$ddf$par["cov1"]
        }
        if(!is.null(mod$ddf$par["cov2"])){
          cov2 <- mod$ddf$par["cov2"]
        }

        aic <- mod$ddf$criterion
        n <- summary(mod$ddf)$n

      }
    }
  }
  res <- rbind(
               c(other.vars, "p", p),
               c(other.vars, "n", n),
               c(other.vars, "b0", b0),
               c(other.vars, "cov1", cov1),
               c(other.vars, "cov2", cov2),
               c(other.vars, "aic", aic),
               c(other.vars, "varp", varp))
               c(other.vars, "varN", varN))
  return(res)
}
