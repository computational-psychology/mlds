library(MLDS)
library(psyphy)
d.df <- read.table('SS_ARu.csv', sep=, header=TRUE)
stim <- sort(unique(c(d.df$s1, d.df$s2, d.df$s3)))
results <- with(d.df, data.frame(resp = Response, S1= match(s1, stim), S2=match(s2, stim), S3=match(s3, stim)))
results <- as.mlbs.df(results, st=stim)

obs.mlds <- mlds(results)
obs.mlds2 <- mlds(results, lnk=probit.2asym(g = 0.0, lam = 0.0))

bias  <- obs.mlds$pscale - obs.mlds2$pscale


print(obs.mlds)
print(obs.mlds2)

resp <- obs.mlds$obj$data$resp
ps <- psyfun.2asym( cbind(resp, 1-resp) ~ . -1 , data= obs.mlds$obj$data, link= probit.2asym)


obs.mlds3 <- mlds(results, lnk=probit.2asym(g = ps$gam, lam = ps$lambda))
print(obs.mlds3)

plot(obs.mlds3)

bias2 <- obs.mlds$pscale - obs.mlds3$pscale

pmc(obs.mlds)
pmc(obs.mlds2)
pmc(obs.mlds3)

par(mfrow = c(1, 3), pty = "s")
plot(obs.mlds)
plot(obs.mlds2)
plot(obs.mlds3)

