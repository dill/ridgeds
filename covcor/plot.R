## summarise and plot results

options(stingsAsFactors=FALSE)
library(ggplot2)
library(plyr)
big.res <- read.csv("covcor.csv")


## quick plot
# load data here
results <- as.data.frame(big.res)
results <- results[,-1]
names(results) <- c("sim","n","model","corr","parameter","value")

# drop some of the results
results <- results[!grepl("hr",results$model),]
#results <- results[results$n!=30,]

## true p lines
#int1 <- integrate(function(x,pars){
#                    exp(-x^2/(2*exp(pars)^2))
#                  },
#                  upper=0.75, lower=0, pars=log(0.3))
#pa.lines <- data.frame(true.p = rep(int1$value/0.75,length(unique(results$n))),
#                       metric = rep("average p",length(unique(results$n))),
#                       n      = unique(results$n))

presults <- results[results$parameter=="p",]
p <- ggplot(presults)
p <- p + geom_boxplot(aes(x=model,y=value))
p <- p + facet_grid(corr~n,scales="free_y")
#p <- p + geom_hline(aes(yintercept=true.p),data=pa.lines)
#p <- p + scale_y_continuous(limits=c(0,1))
p <- p + coord_cartesian(ylim=c(0.25,1))
print(p)

# how do estimates of beta_0 change?
# want to compare b_0 estimates between model with cov1 and cov1+cov2

# true b0?
pars <- c(log(0.75),1.2)
# only want the b0 estimates
b0res <- results[results$parameter=="b0",]
# boxplot of b0
boxb0 <- ggplot(b0res)
boxb0 <- boxb0 + geom_boxplot(aes(x=model,y=value))
boxb0 <- boxb0 + facet_wrap(corr~n,scales="free",nrow=3)
boxb0 <- boxb0 + geom_hline(aes(yintercept=pars[1]))
print(boxb0)


# how do coefficient estimates change per sim between models
# only for cov1 and cov1+cov2 models
b0res <- b0res[b0res$model%in%c("hn+cov1","hn+cov1+cov2"),]
pb <- ggplot(b0res)
pb <- pb + geom_point(aes(x=model,y=value))
pb <- pb + geom_line(aes(x=model,y=value,group=sim),alpha=0.5)
pb <- pb + facet_grid(corr~n,scales="free_y")
print(pb)




#fdiff <- function(x) diff(x$value)
#
#b0dif <- ddply(b0res,.(sim,corr,n),fdiff)
#
## performing a sign test on these differences give non-significant result
#pd <- ggplot(b0dif)
#pd <- pd + geom_histogram(aes(V1))
#pd <- pd + facet_grid(corr~n,scales="free")
#pd <- pd + geom_vline(xintercept=0,colour="red")
#print(pd)



# only want the cov1 beta estimates
cov1res <- results[results$parameter=="cov1",]

# boxplot of b0
boxcov1 <- ggplot(cov1res)
boxcov1 <- boxcov1 + geom_boxplot(aes(x=model,y=value))
boxcov1 <- boxcov1 + facet_wrap(corr~n,scales="free",nrow=3)
boxcov1 <- boxcov1 + geom_hline(aes(yintercept=pars[1]),colour="red")
print(boxcov1)

# only for cov1 and cov1+cov2 models
cov1res <- cov1res[cov1res$model%in%c("hn+cov1","hn+cov1+cov2"),]

pb <- ggplot(cov1res)
pb <- pb + geom_point(aes(x=model,y=value))
pb <- pb + geom_line(aes(x=model,y=value,group=sim),alpha=0.5)
pb <- pb + facet_grid(corr~n,scales="free_y")
print(pb)


## performing a sign test on these differences give non-significant result
#cov1dif <- ddply(cov1res,.(sim,corr,n),fdiff)
#
#pd <- ggplot(cov1dif)
#pd <- pd + geom_histogram(aes(V1),binwidth=0.01)
#pd <- pd + facet_grid(corr~n,scales="free")
#pd <- pd + geom_vline(xintercept=0,colour="red")
#print(pd)




# which was the winner?
aicres <- results[results$parameter=="aic",]

aic.winner <- function(x) as.character(x$model[which.min(x$value)])

aicw <- ddply(aicres,.(sim,corr,n),aic.winner)

pd <- ggplot(aicw)
pd <- pd + geom_histogram(aes(V1))#,binwidth=0.01)
pd <- pd + facet_grid(corr~n,scales="free")
print(pd)

# how many AIC winners were better by more than 3 points than BOTH
# models -- i.e. unambiguous winner
unambig.aic.winner <- function(x){
  this.winner <- which.min(x$value)
  vals <- x$value[-this.winner]-x$value[this.winner]
  #return(all(vals > 3))
  if(any(is.na(vals))){
    return(NA)
  }else if(all(vals > 3)){
    return(as.character(x$model[which.min(x$value)]))
  }else{
    return(NA)
  }
}
unambig.aicw <- ddply(aicres,.(sim,corr,n),unambig.aic.winner)
unambig.aicw <- unambig.aicw[!is.na(unambig.aicw$V1),]

pd <- ggplot(unambig.aicw)
pd <- pd + geom_histogram(aes(V1))#,binwidth=0.01)
pd <- pd + facet_grid(corr~n,scales="free")
print(pd)


# what about variance issues?
pvarresults <- results[results$parameter=="varp",]
p <- ggplot(pvarresults)
p <- p + geom_boxplot(aes(x=model,y=value))
p <- p + facet_grid(corr~n,scales="free_y")
#p <- p + geom_hline(aes(yintercept=true.p),data=pa.lines)
p <- p + scale_y_continuous(limits=c(0,1))
#p <- p + coord_cartesian(ylim=c(0.25,1))
print(p)



