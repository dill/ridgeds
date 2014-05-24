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
results <- results[results$n!=30,]

## true p lines
#int1 <- integrate(function(x,pars){
#                    exp(-x^2/(2*exp(pars)^2))
#                  },
#                  upper=0.75, lower=0, pars=log(0.3))
#pa.lines <- data.frame(true.p = rep(int1$value/0.75,length(unique(results$n))),
#                       metric = rep("average p",length(unique(results$n))),
#                       n      = unique(results$n))


p <- ggplot(results)
p <- p + geom_boxplot(aes(x=model,y=value))
p <- p + facet_grid(parameter~n,scales="free_y")
#p <- p + geom_hline(aes(yintercept=true.p),data=pa.lines)
#p <- p + scale_y_continuous(limits=c(0,1))
print(p)


# how do estimates of beta_0 change?
# want to compare b_0 estimates between model with cov1 and cov1+cov2

# only want the b0 estimates
b0res <- results[results$parameter=="b0",]
# only for cov1 and cov1+cov2 models
b0res <- b0res[b0res$model%in%c("hn+cov1","hn+cov1+cov2"),]

pb <- ggplot(b0res)
pb <- pb + geom_point(aes(x=model,y=value))
pb <- pb + geom_line(aes(x=model,y=value,group=sim),alpha=0.5)
pb <- pb + facet_grid(corr~n,scales="free_y")
print(pb)



fdiff <- function(x) diff(x$value)

b0dif <- ddply(b0res,.(sim,corr,n),fdiff)

# performing a sign test on these differences give non-significant result
pd <- ggplot(b0dif)
pd <- pd + geom_histogram(aes(V1))
pd <- pd + facet_grid(corr~n,scales="free")
pd <- pd + geom_vline(xintercept=0,colour="red")
print(pd)



# only want the cov1 beta estimates
cov1res <- results[results$parameter=="cov1",]
# only for cov1 and cov1+cov2 models
cov1res <- cov1res[cov1res$model%in%c("hn+cov1","hn+cov1+cov2"),]

pb <- ggplot(cov1res)
pb <- pb + geom_point(aes(x=model,y=value))
pb <- pb + geom_line(aes(x=model,y=value,group=sim),alpha=0.5)
pb <- pb + facet_grid(corr~n,scales="free_y")
print(pb)


# performing a sign test on these differences give non-significant result
cov1dif <- ddply(cov1res,.(sim,corr,n),fdiff)

pd <- ggplot(cov1dif)
pd <- pd + geom_histogram(aes(V1),binwidth=0.01)
pd <- pd + facet_grid(corr~n,scales="free")
pd <- pd + geom_vline(xintercept=0,colour="red")
print(pd)




# which was the winner?
aicres <- results[results$parameter=="aic",]

aic.winner <- function(x) as.character(x$model[which.min(x$value)])

aicw <- ddply(aicres,.(sim,corr,n),aic.winner)

pd <- ggplot(aicw)
pd <- pd + geom_histogram(aes(V1),binwidth=0.01)
pd <- pd + facet_grid(corr~n,scales="free")
print(pd)

