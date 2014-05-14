sim.hn <- function(n,z,beta,width){
  sigma <- function(z,beta){
    # z has n rows (# obs)
    #       p cols (length beta)
    exp(as.matrix(z)%*%beta)
  }

  d <- abs(rnorm(2*n, 0, sigma(z,beta)))
  d <- d[d <= width] # only obsns within truncation
  d[1:n] # we generated 2n samples, just take the first n
}


