# GMU 13E Black Bear data
# Adapted from a script provided by Earl Becker
# Alaska Department of Fish and Game

# file creates the R data file

library(optimx)
library(mrds)
library(stats)

## load data
BlackLTdata13E <- read.csv("GMU13E_BlackLTdata.csv", header=TRUE,row.names=1)

## recoding
# percentage snow
BlackLTdata13E$PchSnow <- ifelse(BlackLTdata13E$PCTSNOW>50, (100-BlackLTdata13E$PCTSNOW), BlackLTdata13E$PCTSNOW)
# group size
BlackLTdata13E$AdjGroupSize <- BlackLTdata13E$GROUPSIZE
BlackLTdata13E$AdjGroupSize[BlackLTdata13E$GROUPTYPE == "female w/ coy(s)"] <- 1
BlackLTdata13E$AdjGroupSize[BlackLTdata13E$GROUPTYPE == "female w/ yearling cub(s)"] <-
  (0.5+(0.5*BlackLTdata13E$GROUPSIZE[BlackLTdata13E$GROUPTYPE == "female w/ yearling cub(s)"]))

#  table(BlackLTdata13E$ACTIVITY)  # 108 bedded down black bears
#  No sitting black bears observered!
BlackLTdata13E$Bedded <- ifelse(BlackLTdata13E$ACTIVITY =="bedded dow", 1, 0)

BlackLTdata13E$PatchSnow <- ifelse(BlackLTdata13E$PCTSNOW>=30 & BlackLTdata13E$PCTSNOW <= 70, 1, 0)
BlackLTdata13E$HiSnowCov <-  ifelse(BlackLTdata13E$PCTSNOW>=70, 1, 0)

cat("BlackLTdata13E contains",dim(BlackLTdata13E)[1],"bears observations\n")
cat("   ",sum(BlackLTdata13E$distance < 450),"/",dim(BlackLTdata13E)[1],
    "(", 100*sum(BlackLTdata13E$distance < 450)/dim(BlackLTdata13E)[1],"%) are in the 22-450m range \n")


## now drop some unused columns

BlackLTdata13E$BEARID <- NULL
BlackLTdata13E$ACTIVITY <- NULL
BlackLTdata13E$GROUPTYPE <- NULL

par(mfrow=c(1,2))
plot(BlackLTdata13E$distance,BlackLTdata13E$EFF2TRANS,xlab="Distance",ylab="Search distance")

# throw a loess through that?
ll <- loess(EFF2TRANS~distance,data=BlackLTdata13E)
pr <- predict(ll,newdata=data.frame(distance=seq(0,700,len=1500)))
lines(seq(0,700,len=1500),pr,col="red",lwd=2)

plot(BlackLTdata13E$distance,log(BlackLTdata13E$EFF2TRANS),xlab="Distance",ylab="log(Search distance)")
ll <- loess(log(EFF2TRANS)~distance,data=BlackLTdata13E)
pr <- predict(ll,newdata=data.frame(distance=seq(0,700,len=1500)))
lines(seq(0,700,len=1500),pr,col="red",lwd=2)

# save that data
save(BlackLTdata13E,file="dave-13e.RData")


