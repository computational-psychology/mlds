library(MLDS)

d.df <- read.table('first.csv', sep=, header=TRUE)  
d.df2 <- read.table('second.csv', sep=, header=TRUE)  
d.all <- rbind(d.df, d.df2)

stim <- sort(unique(c(d.df$s1, d.df$s2, d.df$s3)))

d1 <- with(d.df, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))
d1 <- as.mlbs.df(d1, st = stim)

d2 <- with(d.df2, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))
d2 <- as.mlbs.df(d2, st = stim)

dall <- with(d.all, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))
dall <- as.mlbs.df(dall, st = stim)

# not same results when sigma is very small
m1.glm <- mlds(d1, lnk=probit.2asym(g = 0.0, lam = 0.0))
m2.glm <- mlds(d2, lnk=probit.2asym(g = 0.0, lam = 0.0))

m1.glm2 <- mlds(d1, lnk = "probit")
m2.glm2 <- mlds(d2, lnk = "probit")


bias1 <- m1.glm$pscale - m1.glm2$pscale
bias2 <- m2.glm$pscale - m2.glm2$pscale


anova(m1.glm$obj, m2.glm$obj)
anova(m1.glm2$obj, m2.glm2$obj)

m.glm <- mlds(dall, lnk = "probit")

#########
library(mgcv)


dfst <- with(d.all, data.frame(resp = Response, S1 = s1, S2 = s2, S3 = s3))
dfst <- as.mlbs.df(dfst, st = stim)  # not nec. except to attach invord attr
dfst[attr(dfst, "invord"), -1] <- dfst[attr(dfst, "invord"), 4:2]
dfst[attr(dfst, "invord"), 1] <- 1 - dfst[attr(dfst, "invord"), 1]

k <- 4

by.mat <- matrix(c(1, -2, 1), nc = 3, nr = nrow(dfst), byrow = TRUE)
S.mat <- with(dfst, cbind(S1 = S1, S2 = S2, S3 = S3))
m.gam <- gam(resp ~ s(S.mat, k = k, by = by.mat), family = binomial(probit), data = dfst)

pang <- seq(0, 1, len = 100)  # stimulus vector
nd <-  list(S.mat = cbind(S1 = pang, S2 = 0, S3 = 0), by.mat = matrix(1, nc = 3, nr = length(pang))) # new data object
m.pred <- predict(m.gam, newdata = nd, type = "link") # predict from m.gam object, using new data
m.pred <- m.pred - mean(m.pred)  # centering at zero 

plot(m.gam, shade = TRUE, 
     ylim = c(-1, 4) - mean(coef(m.glm)),
     shift = mean(coef(m.glm)), 
     xlab = "Stimulus", ylab = "Difference Scale")

# + mean(coef(m1.glm))   --> puts minimum at zero
lines(pang, m.pred + mean(c(0, coef(m.glm))), lty = 2, lwd = 4, col = "orange")  
points(m.glm, pch = 21, bg = "white")
# 
# plot(pang, m.pred + mean(coef(m.glm))) # this command to be used

# splitting into two
first <- rbind( by.mat[1:nrow(d1),] , by.mat[(nrow(d1)+1):nrow(d.all), ] * 0)
second <- rbind( by.mat[1:nrow(d1),]*0 , by.mat[(nrow(d1)+1):nrow(d.all), ] )

m.gam2 <- gam(resp ~ s(S.mat, k = k, by = first) + s(S.mat, k = k, by = second), 
              family = binomial(probit), data = dfst)
summary(m.gam2)
anova(m.gam, m.gam2, test = "Chisq")    # ANOVA correctly detect a better fit if separated.

plot(m.gam2)


nd1 <- list(S.mat = cbind(S1 = pang, S2 = 0, S3 = 0), 
            first = matrix(1, nc = 3, nr = length(pang)),
            second = matrix(0, nc = 3, nr = length(pang))
)
nd2 <- list(S.mat = cbind(S1 = pang, S2 = 0, S3 = 0), 
            first = matrix(0, nc = 3, nr = length(pang)),
            second = matrix(1, nc = 3, nr = length(pang))
)
pr1 <- predict(m.gam2, newdata = nd1, type = "link")
pr2 <- predict(m.gam2, newdata = nd2, type = "link")
pr1 <- pr1 - mean(pr1)
pr2 <- pr2 - mean(pr2)


par(mfrow = c(1, 2), pty = "s")
plot(m.gam2, select = 1, shift = -pr1[1])
lines(pang, pr1-pr1[1], col = "blue")
plot(m.gam2, select = 2, shift = -pr2[1])
lines(pang, pr2-pr2[1], col = "blue")


#
par(mfrow = c(1, 2), pty = "s")
plot(pang, pr1-pr1[1], col = "blue", ylim = c(0, 3) )
par(new=T)
plot(m1.glm, ylim = c(0, 3))

#
plot(pang, pr2-pr2[1], col = "red", ylim = c(0, 3))
par(new=T)
plot(m2.glm, ylim = c(0, 3))

anova(m.gam, m.gam2, test = "Chisq")    # ANOVA correctly detect a better fit if separated.
