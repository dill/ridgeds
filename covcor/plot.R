options(stingsAsFactors=FALSE)
library(ggplot2)
big.res <- read.csv("covcor.csv")


## quick plot
# load data here
results <- as.data.frame(big.res)
results <- results[,-1]
names(results) <- c("sim","n","model","p","varp")

# append the variances to the bottom of the results data.frame
results <- rbind(results,results)
results$p[(nrow(results)/2+1):nrow(results)] <- results$varp[(nrow(results)/2+1):nrow(results)]
results$varp <- NULL
results$metric <- c(rep("average p",nrow(results)/2),
                    rep("se average p",nrow(results)/2))

# drop some of the results
results <- results[!grepl("hr",results$model),]
results <- results[results$n!=30,]

results$p <- as.numeric(as.character(results$p))
#results$varp <- as.numeric(as.character(results$varp))
results$n <- as.factor(as.numeric(as.character(results$n)))

# true p lines
int1 <- integrate(function(x,pars){
                    exp(-x^2/(2*exp(pars)^2))
                  },
                  upper=0.75, lower=0, pars=log(0.3))
pa.lines <- data.frame(true.p = rep(int1$value/0.75,length(unique(results$n))),
                       metric = rep("average p",length(unique(results$n))),
                       n      = unique(results$n))


p <- ggplot(results)
p <- p + geom_boxplot(aes(x=model,y=p))
p <- p + facet_grid(metric~n,scales="free_y")
p <- p + geom_hline(aes(yintercept=true.p),data=pa.lines)
#p <- p + scale_y_continuous(limits=c(0,1))
print(p)



