# Try a simple analysis, ignoring MR component

library(Distance)

load("dave-13e.RData")

## this fit is horrible as we have a spike away from zero (22)
#bl.hn.nosd <- ds(BlackLTdata13E, truncation=list(left=22,right=450))
## anything with adjustments results in functions with P_a>1
#bl.hr.nosd <- ds(BlackLTdata13E, truncation=list(left=22,right=450),key="hr")
## try with no adjustments
#bl.hr.nosd <- ds(BlackLTdata13E, truncation=list(left=22,right=450),key="hr",order=0)
# So: let's truncate heavily and standardise the search distance

BlackLTdata13E$effscale <- BlackLTdata13E$EFF2TRANS/sd(BlackLTdata13E$EFF2TRANS)

# truncate but manually on the data
bl.dat <- BlackLTdata13E[BlackLTdata13E$distance>99,]
bl.dat$distance <- bl.dat$distance-99
bl.dat <- bl.dat[,c("distance","effscale","PCTCOVER","PCTSNOW")]

# standardise the covariates
bl.dat[,2:4] <- apply(bl.dat[,2:4],2,scale)

bl.nosd <- ds(bl.dat, truncation=list(left=0,right=450))
bl.allcov <- ds(bl.dat, truncation=list(left=0,right=450),formula=~effscale+PCTCOVER+PCTSNOW)

bl.PCA <- as.matrix(bl.dat[,2:4])
pc <- eigen((t(bl.PCA)%*%bl.PCA))
pc.dat <- bl.PCA%*%pc$vectors
pc.dat <- as.data.frame(pc.dat)
names(pc.dat) <- paste0("cov",1:3)
pc.dat$distance <- bl.dat$distance


bl.pca1 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1)
bl.pca12 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1+cov2)
bl.pca123 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1+cov2+cov3)

#par(mfrow=c(1,4))
#plot(bl.allcov,main="all covariates")
#plot(bl.pca1,main="PCA1")
#plot(bl.pca12,main="PCA12")
#plot(bl.pca123,main="PCA123")
#
#
#summary(bl.allcov)
#summary(bl.pca1)
#summary(bl.pca12)
#summary(bl.pca123)


# print results
res <- data.frame(p = c(summary(bl.allcov)$ds$average.p,
                        summary(bl.pca1)$ds$average.p,
                        summary(bl.pca12)$ds$average.p,
                        summary(bl.pca123)$ds$average.p),
                  p.se = c(summary(bl.allcov)$ds$average.p.se,
                           summary(bl.pca1)$ds$average.p.se,
                           summary(bl.pca12)$ds$average.p.se,
                           summary(bl.pca123)$ds$average.p.se),
                  AIC = c(bl.allcov$ddf$criterion,
                          bl.pca1$ddf$criterion,
                          bl.pca12$ddf$criterion,
                          bl.pca123$ddf$criterion))
rownames(res) <- c("allcov",paste0("pca",1:3))

print(res)

#summary(bl.allcov)$ds$average.p.se/summary(bl.pca1)$ds$average.p.se




#### very correlated
bl.dat <- BlackLTdata13E[BlackLTdata13E$distance>99,]
bl.dat$distance <- bl.dat$distance-99

## add an extra covariate 0.9 correlated w eff
corr <- 0.9
covcor <- corr*bl.dat$effscale + sqrt(1-corr^2)*rnorm(nrow(bl.dat),mean(bl.dat$effscale,sd(bl.dat$effscale)))

#print(cor(bl.dat$effscale,covcor))

bl.dat <- bl.dat[,c("distance","effscale","PCTCOVER","PCTSNOW")]
bl.dat$covcor <- covcor

# standardise the covariates
bl.dat[,2:5] <- apply(bl.dat[,2:5],2,scale)

#bl.nosd <- ds(bl.dat, truncation=list(left=0,right=450))
bl.allcov <- ds(bl.dat, truncation=list(left=0,right=450),formula=~effscale+PCTCOVER+PCTSNOW+covcor)

bl.allcov.t <- ds(bl.dat, truncation=list(left=0,right=450),formula=~effscale+PCTCOVER+PCTSNOW)

bl.PCA <- as.matrix(bl.dat[,2:5])
pc <- eigen((t(bl.PCA)%*%bl.PCA))
pc.dat <- bl.PCA%*%pc$vectors
pc.dat <- as.data.frame(pc.dat)
names(pc.dat) <- paste0("cov",1:4)
pc.dat$distance <- bl.dat$distance


bl.pca1 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1)
bl.pca12 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1+cov2)
bl.pca123 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1+cov2+cov3)
bl.pca1234 <- ds(pc.dat, truncation=list(left=0,right=450),formula=~cov1+cov2+cov3+cov4)
#
#par(mfrow=c(1,4))
#plot(bl.allcov,main="all covariates")
#plot(bl.pca1,main="PCA1")
#plot(bl.pca12,main="PCA12")
#plot(bl.pca123,main="PCA123")
#
#
#summary(bl.allcov)
#summary(bl.pca1)
#summary(bl.pca12)
#summary(bl.pca123)
#
#
#summary(bl.allcov)$ds$average.p.se/summary(bl.pca1)$ds$average.p.se
# print results
res <- data.frame(p = c(summary(bl.allcov)$ds$average.p,
                        summary(bl.allcov.t)$ds$average.p,
                        summary(bl.pca1)$ds$average.p,
                        summary(bl.pca12)$ds$average.p,
                        summary(bl.pca123)$ds$average.p,
                        summary(bl.pca1234)$ds$average.p),
                  p.se = c(summary(bl.allcov)$ds$average.p.se,
                           summary(bl.allcov.t)$ds$average.p.se,
                           summary(bl.pca1)$ds$average.p.se,
                           summary(bl.pca12)$ds$average.p.se,
                           summary(bl.pca123)$ds$average.p.se,
                           summary(bl.pca1234)$ds$average.p.se),
                  AIC = c(bl.allcov$ddf$criterion,
                          bl.allcov.t$ddf$criterion,
                          bl.pca1$ddf$criterion,
                          bl.pca12$ddf$criterion,
                          bl.pca123$ddf$criterion,
                          bl.pca1234$ddf$criterion))
rownames(res) <- c("allcov","allcov.t",paste0("pca",1:4))

print(res)
