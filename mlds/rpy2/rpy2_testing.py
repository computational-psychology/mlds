# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:56:35 2014

@author: guille
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import mlds

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
MLDS = importr("MLDS")
mgcv = importr("mgcv")

r = robjects.r

file1 = 'first.csv'
file2 = 'second.csv'
file3 = 'third.csv'

obs1 = mlds.MLDSObject( file1, standardscale =False, boot=False)
obs2 = mlds.MLDSObject( file2, standardscale =False, boot=False)
obs3 = mlds.MLDSObject( file2, standardscale =False, boot=False)
obs1.run()
obs2.run()
obs3.run()


##

d_df1 = r('d.df1 <- read.table("%s", sep=, header=TRUE)' % file1)
d_df2 = r('d.df2 <- read.table("%s", sep=, header=TRUE)' % file2)
d_all = r('d.all <- rbind(d.df1, d.df2)')

stim = r('stim <- sort(unique(c(d.all$s1, d.all$s2, d.all$s3)))')

r('d1 <- with(d.df1, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))')
r('d2 <- with(d.df2, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))')
r('dall <- with(d.all, data.frame(resp = Response, S1 = i1, S2 = i2, S3 = i3))')

d1 = r('d1 <- as.mlbs.df(d1, st = stim)')
d2 = r('d2 <- as.mlbs.df(d2, st = stim)')
dall = r('dall <- as.mlbs.df(dall, st = stim)')

# MLDS calculated via GLM
m1_glm = r('m1.glm <- mlds(d1, lnk="probit")')
m2_glm = r('m2.glm <- mlds(d2, lnk="probit")')
m_glm = r('m.glm <- mlds(dall, lnk="probit")')


########## GAM

k = r('k <- 4') # smoothing parameter

# preparing dataframe

dfst = r('dfst <- with(d.all, data.frame(resp = Response, S1 = s1, S2 = s2, S3 = s3))')
dfst = r('dfst <- as.mlbs.df(dfst, st = stim) ')
dfst = r('dfst[attr(dfst, "invord"), -1] <- dfst[attr(dfst, "invord"), 4:2]') # reverting the stimulus order according to invord attribute
dfst = r('dfst[attr(dfst, "invord"), 1] <- 1 - dfst[attr(dfst, "invord"), 1]') # reverting responses according to invord attribute


# gam call
by_mat = r('by.mat <- matrix(c(1, -2, 1), nc = 3, nr = nrow(dfst), byrow = TRUE)')
S_mat =  r('S.mat <- with(dfst, cbind(S1 = S1, S2 = S2, S3 = S3))')
m_gam =  r('m.gam <- gam(resp ~ s(S.mat, k = k, by = by.mat), family = binomial(probit), data = dfst)')
print m_gam

# 
svec = r('svec <- seq(0, 1, len = 100)')  # stimulus vector
nd = r('nd <-  list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), by.mat = matrix(1, nc = 3, nr = length(svec)))') # new data object
m_pred = r('m.pred <- predict(m.gam, newdata = nd, type = "link")') # predict from m.gam object, using new data
m_pred = r('m.pred <- m.pred - mean(m.pred)')  # centering at zero 
#m_zero = r('m.pred + mean(c(0, coef(m.glm)))')  # zero as anchor at stim = zero
m_zero = r('m.pred - m.pred[1]')  # zero as anchor at stim = zero



# evaluating separatedely
first = r('first <- rbind( by.mat[1:nrow(d1),] , by.mat[(nrow(d1)+1):nrow(d.all), ] * 0)')
second = r('second <- rbind( by.mat[1:nrow(d1),]*0 , by.mat[(nrow(d1)+1):nrow(d.all), ] )')

m_gam2 = r('m.gam2 <- gam(resp ~ s(S.mat, k = k, by = first) + s(S.mat, k = k, by = second), family = binomial(probit), data = dfst)')

rsummary = r('summary')
ranova = r('anova')

print rsummary(m_gam2)
anovares = ranova(m_gam, m_gam2, test = "Chisq") 
print anovares   # ANOVA correctly detect a better fit if separated.
print anovares[4]  # p- value


nd1 = r('nd1 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), first = matrix(1, nc = 3, nr = length(svec)), second = matrix(0, nc = 3, nr = length(svec)))')
nd2 = r('nd2 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), first = matrix(0, nc = 3, nr = length(svec)), second = matrix(1, nc = 3, nr = length(svec)))')


m1_pred = r('m1.pred <- predict(m.gam2, newdata = nd1, type = "link")') # predict from m.gam object, using new data
m1_pred = r('m1.pred <- m1.pred - mean(m1.pred)')  # centering at zero 
#m1_zero = r('m1.pred + mean(c(0, coef(m1.glm)))')  # zero as anchor
m1_zero = r('m1.pred - m1.pred[1]')                # zero as anchor 


m2_pred = r('m2.pred <- predict(m.gam2, newdata = nd2, type = "link")') # predict from m.gam object, using new data
m2_pred = r('m2.pred <- m2.pred - mean(m2.pred)')  # centering at zero 
#m2_zero = r('m2.pred + mean(c(0, coef(m2.glm)))')  # zero as anchor a
m2_zero = r('m2.pred - m2.pred[1]')  # zero as anchor 


##
plt.figure()
plt.plot(np.array(svec), np.array(m_zero), label="GAM all", linewidth=2)
plt.plot(np.array(svec), np.array(m1_zero), label="GAM first", linewidth=2)
plt.plot(np.array(svec), np.array(m2_zero), label="GAM second", linewidth=2)
mlds.plotscale(obs1, observer="GLM first", color='green', marker='o', linewidth=0)
mlds.plotscale(obs2, observer="GLM second", color='red', marker='o', linewidth=0)
plt.legend(loc=2)
plt.show()

