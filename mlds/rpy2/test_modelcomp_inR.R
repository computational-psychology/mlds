library(MLDS)
library(psyphy)

d.df <- read.table('first.csv', sep=, header=TRUE)  
d.df2 <- read.table('second.csv', sep=, header=TRUE)  

stim <- sort(unique(c(d.df$s1, d.df$s2, d.df$s3)))

d1 <- with(d.df, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))
d1 <- as.mlbs.df(d1, st = stim)

d2 <- with(d.df2, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))
d2 <- as.mlbs.df(d2, st = stim)


# not same results when sigma is very small
m1.glm <- mlds(d1, lnk=probit.2asym(g = 0.0, lam = 0.0))
m2.glm <- mlds(d2, lnk=probit.2asym(g = 0.0, lam = 0.0))

m1.glm2 <- mlds(d1, lnk = "probit")
m2.glm2 <- mlds(d2, lnk = "probit")


bias1 <- m1.glm$pscale - m1.glm2$pscale
bias2 <- m2.glm$pscale - m2.glm2$pscale


anova(m1.glm$obj, m2.glm$obj)
anova(m1.glm2$obj, m2.glm2$obj)


#########
library(mgcv)

dfst <- with(d.df, data.frame(resp = Response, S1 = s1, S2 = s2, S3 = s3))
dfst <- as.mlbs.df(dfst, st = stim)  # not nec. except to attach invord attr
### from mlds.mlbs.df to effect SwapOrder on triads
dfst[attr(dfst, "invord"), -1] <- dfst[attr(dfst, "invord"), 4:2]
dfst[attr(dfst, "invord"), 1] <- 1 - dfst[attr(dfst, "invord"), 1]

k = 4

by.mat <- matrix(c(1, -2, 1), nc = 3, nr = nrow(dfst), byrow = TRUE)
S.mat <- with(dfst, cbind(S1 = S1, S2 = S2, S3 = S3))
m.gam <- gam(resp ~ s(S.mat, k = k, by = by.mat), family = binomial(probit), data = dfst)

plot(m.gam, shade = TRUE, 
     ylim = c(-1, 5) - mean(coef(m1.glm)),
     shift = mean(coef(m1.glm)), 
     xlab = "Slant (deg)", ylab = "Difference Scale")

