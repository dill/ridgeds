sim.hn <- function(n,pars,width){
  sigma <- function(z,beta){
    # z has n rows (# obs)
    #       p cols (length beta)
    exp(as.matrix(z)%*%beta)
  }

  # need to integrate these into the rest of the code?
  detfct.tmp<-function(x,pars,z=NULL){
    # evaluate the total detection function for either hn or hr keys - no error
    # checking on inputs, assumes that's been done earlier
    sigmas<-sigma(z,pars)
    return(mrds:::keyfct.hn(x,sigmas))
  }

  dat <- matrix(NA,n,2)

counter<-0
  while(counter<n){

a<-0.45;b<-0.45

      # for covariate models, set z for each observation
      z <- t(c(1,rbeta(1,a,b)))

      proposal<-runif(1,0,width)
      U<-runif(1)

      # accept/reject
      if(U<=detfct.tmp(proposal,pars,z)){
         counter<-counter+1
         dat[counter,]<-c(proposal,z[2])
      }

   }
  return(dat)
}


