# shortcut function to store results
store_results <- function(mod, model.name, sample.size, sim, res){
  # mod         - fitted model object (call was wrapped in a try())
  # model.name  - string to identify model
  # sample.size - number of observations
  # sim         - simulation id
  # results     - data.frame to append to

  # dummy
  p <- NA
  varp <- NA

  # if fitting went okay, then store some values
  if(!is.null(mod)){
    if(class(mod)!="try-error"){
      if(!is.null(mod$ddf)){
        p <- sample.size/mod$ddf$Nhat
        varp <- summary(mod$ddf)$average.p.se
      }
    }
  }
  res <- rbind(res, c(sim, sample.size, model.name, p, varp))
  return(res)
}
