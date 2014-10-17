# -*- coding: utf-8 -*-
"""
Utilities for MLDS estimation using GAM. Interface with R via rpy2 module.

@author: G. Aguilar, Oct 2014.

TODO: improve data preparation for multiple files.

"""

import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
MLDS = importr("MLDS")
mgcv = importr("mgcv")
r = robjects.r

class MLDSGAMCompare:
    """
    MLDSGAMCompare object reads files of triad experiment, calculate difference
    scales using GAM method, and compares them statistically.
    
    * Usage
    
    * Methods

    * Attributes
    
    * References:
        Wood 
        
    """
    
    def __init__(self, files):
        
        self.files = files
        self.scales    = []
        self.scale_all = []
        self.stim = []
        self.anovares = []
        
    
    def run(self):
        """
        Calculate difference scales using GAM for each file, and
        run an analysis of deviance to establish signifficant differences.
        
        """
        
        ##
        l = []
        nrows = []
        for i,file in enumerate(self.files):
            # read data from file to dataframes
            df= r('df%d <- read.table("%s", sep=, header=TRUE)' % (i, file))
            nrows.append( list(r['nrow'](df))[0] )
            l.append('df%d' % i)
            
        d_all = r('dfall <- rbind(%s)' % ",".join(l))
        stim = r('stim <- sort(unique(c(dfall$s1, dfall$s2, dfall$s3)))')    
        nrowall = list(r['nrow'](d_all))[0]
        assert(nrowall == sum(nrows))
        
        
        ########## GAM
        k = r('k <- 4') # smoothing parameter
        
        # preparing dataframe
        dfst = r('dfst <- with(dfall, data.frame(resp = Response, S1 = s1, S2 = s2, S3 = s3))')
        dfst = r('dfst <- as.mlbs.df(dfst, st = stim) ')
        dfst = r('dfst[attr(dfst, "invord"), -1] <- dfst[attr(dfst, "invord"), 4:2]') # reverting the stimulus order according to invord attribute
        dfst = r('dfst[attr(dfst, "invord"), 1] <- 1 - dfst[attr(dfst, "invord"), 1]') # reverting responses according to invord attribute
        
        
        # gam call for all data together
        by_mat = r('by.mat <- matrix(c(1, -2, 1), nc = 3, nr = nrow(dfst), byrow = TRUE)')
        S_mat =  r('S.mat <- with(dfst, cbind(S1 = S1, S2 = S2, S3 = S3))')
        m_gam =  r('m.gam <- gam(resp ~ s(S.mat, k = k, by = by.mat), family = binomial(probit), data = dfst)')
        print m_gam
        
        # 
        svec = r('svec <- seq(0, 1, len = 100)')  # stimulus vector
        nd = r('nd <-  list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), by.mat = matrix(1, nc = 3, nr = length(svec)))') # new data object
        
        m_pred = r('m.pred <- predict(m.gam, newdata = nd, type = "link")') # predict from m.gam object, using new data
        m_pred = r('m.pred <- m.pred - mean(m.pred)')  # centering at zero 
        m_zero = np.array(r('m.pred - m.pred[1]'))  # zero as anchor at stim = zero


        ### evaluating separatedely. This preparation code could be done better,
        # but I cannot find a way to  do it cleaner and dynamic.. 
        # instead I'm  writing all the conditions. not wasting time on this, maybe in the future
        s = np.concatenate(([0], np.cumsum(nrows[:-1])))
        e = s+nrows
        s=s+1
        
        block= 'by.mat[%d:%d,]'
        noblock = 'by.mat[%d:%d, ] * 0'
        
        if len(self.files)==4:
            
            sec0 = r('sec0 <- rbind( %s , %s, %s, %s )' % ( block % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]), 
                                                        noblock % (s[2], e[2]),
                                                        noblock % (s[3], e[3])  ) )
                                                        
            sec1 = r('sec1 <- rbind( %s , %s, %s, %s )' % ( noblock % (s[0], e[0]), 
                                                        block % (s[1], e[1]), 
                                                        noblock % (s[2], e[2]),
                                                        noblock % (s[3], e[3])  ) )
                                                        
            sec2 = r('sec2 <- rbind( %s , %s, %s, %s )' % ( noblock % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]), 
                                                        block % (s[2], e[2]),
                                                        noblock % (s[3], e[3])  ) )
                                                        
            sec3 = r('sec3 <- rbind( %s , %s, %s, %s )' % ( noblock % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]), 
                                                        noblock % (s[2], e[2]),
                                                        block % (s[3], e[3])  ) )                                            
                                                        
            gamcall = 'm.gam2 <- gam(resp ~ s(S.mat, k = k, by = sec0) + s(S.mat, k = k, by = sec1)+ s(S.mat, k = k, by = sec2) + s(S.mat, k = k, by = sec3), family = binomial(probit), data = dfst)'                                           
               
            nd0 = r('nd0 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(1, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)), sec2 = matrix(0, nc = 3, nr = length(svec)), sec3 = matrix(0, nc = 3, nr = length(svec)))')
            nd1 = r('nd1 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(1, nc = 3, nr = length(svec)), sec2 = matrix(0, nc = 3, nr = length(svec)), sec3 = matrix(0, nc = 3, nr = length(svec)))')
            nd2 = r('nd2 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)), sec2 = matrix(1, nc = 3, nr = length(svec)), sec3 = matrix(0, nc = 3, nr = length(svec)))')
            nd3 = r('nd3 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)), sec2 = matrix(0, nc = 3, nr = length(svec)), sec3 = matrix(1, nc = 3, nr = length(svec)))')
                  
               
        elif len(self.files)==3:
            
            sec0 = r('sec0 <- rbind( %s , %s, %s )' % ( block % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]), 
                                                        noblock % (s[2], e[2])  ) )
                                                        
            sec1 = r('sec1 <- rbind( %s , %s, %s )' % ( noblock % (s[0], e[0]), 
                                                        block % (s[1], e[1]), 
                                                        noblock % (s[2], e[2])  ) )
                                                        
            sec2 = r('sec2 <- rbind( %s , %s, %s )' % ( noblock % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]), 
                                                        block % (s[2], e[2])  ) )
            gamcall = 'm.gam2 <- gam(resp ~ s(S.mat, k = k, by = sec0) + s(S.mat, k = k, by = sec1)+ s(S.mat, k = k, by = sec2), family = binomial(probit), data = dfst)'                                           
             
            nd0 = r('nd0 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(1, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)), sec2 = matrix(0, nc = 3, nr = length(svec)))')
            nd1 = r('nd1 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(1, nc = 3, nr = length(svec)), sec2 = matrix(0, nc = 3, nr = length(svec)))')
            nd2 = r('nd2 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)), sec2 = matrix(1, nc = 3, nr = length(svec)))')
                                                     
        elif len(self.files)==2:
            
            sec0 = r('sec0 <- rbind( %s , %s )' % ( block % (s[0], e[0]), 
                                                        noblock % (s[1], e[1]) ) )
                                                        
            sec1 = r('sec1 <- rbind( %s , %s )' % ( noblock % (s[0], e[0]), 
                                                        block % (s[1], e[1])  ) )
                                                        
            gamcall = 'm.gam2 <- gam(resp ~ s(S.mat, k = k, by = sec0) + s(S.mat, k = k, by = sec1), family = binomial(probit), data = dfst)'                                           
                
            
            nd0 = r('nd0 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(1, nc = 3, nr = length(svec)), sec1 = matrix(0, nc = 3, nr = length(svec)))')
            nd1 = r('nd1 <- list(S.mat = cbind(S1 = svec, S2 = 0, S3 = 0), sec0 = matrix(0, nc = 3, nr = length(svec)), sec1 = matrix(1, nc = 3, nr = length(svec)))')
        
        
        else:
            raise('not implemented for more than 4 files (yet)')
        
        
        # calling GAM                         
        m_gam2 = r(gamcall)
        
        ranova = r('anova')
        
        print r['summary'](m_gam2)
        anovares = ranova(m_gam, m_gam2, test = "Chisq") 
        print anovares   # ANOVA correctly detect a better fit if separated.
        print anovares[4]  # p- value
        
        def predict(i):
            r('m0.pred <- predict(m.gam2, newdata = nd%d, type = "link")' % i) # predict from m.gam object, using new data
            r('m0.pred <- m0.pred - mean(m0.pred)')  # centering at zero 
            #m_zero = r('m0.pred + mean(c(0, coef(m0.glm)))')  # zero as anchor
            m_zero = r('m0.pred - m0.pred[1]')                # zero as anchor 
            return np.array(m_zero)
            
        
        mzeros = np.zeros((len(svec), len(self.files)))
        for i in range(len(self.files)):
            mzeros[:,i] = predict(i)
        
        # saving to object self
        self.stim = svec
        self.scale_all = m_zero
        self.scales = mzeros
        self.anovares = anovares
        
    
