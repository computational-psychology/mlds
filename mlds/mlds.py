# -*- coding: utf-8 -*-
"""
Utilities for MLDS estimation, interfacing with R.

@author: G. Aguilar, Nov 2013, rev. Apr 2014
"""

import numpy as np
import os
import csv
import itertools
import random
import subprocess
import uuid


class MLDSObject:
    """
    MLDS Object allows the quick analysis of perceptual scales from triads or quadruples experiments
    using the maximum-likelihood difference scaling method (Maloney & Yang, 2003). 
    
    It runs the analysis on R, using the MLDS package.
    It requires, therefore, a running R version with the package installed.
    
    Usage
    ----------
    
    obs = MLDSObject( datafile, *flags  )
   
    boolean flags: 
    
    boot -> calculate bootstrap estimates (default : False) \n
    keepfiles -> keep temporary R and csv files (default: False) \n
    standardscale -> calculate perceptual scale in range [0,1] (default: True)\n
    getlinearscale -> calculate linear scale as if differencing perfecly (default=False)\n
   
    Methods
    ----------
    
    obs.run()  --> run MLDS analysis


    Attributes
    ----------
    
    obs.stim  --> stimulus scale\n
    obs.scale --> estimated perceptual scale\n
    obs.sigma --> estimated noise parameter\n
    obs.lscale --> linear scale (assuming perfect differentiating)\n
    
    when boot = True, \n
    obs.mns and obs.ci95  --> mean and 95% confidence interval for scale estimates after bootstrap.         
    obs.sigmamns and obs.sigmaci95 -->   "                     for sigma estimate after bootstrap. 
    """
    
    def __init__(self, filename, boot=False, keepfiles=False, standardscale=True, getlinearscale=False, verbose=False):
        
        self.status=0 # 0: just initialized, 1: mlds computed, 2: bootstrapped, -1: error
        self.verbose = verbose
        
        # flags
        self.boot = boot
        self.keepfiles = keepfiles
        self.standardscale = standardscale
        self.getlinearscale = getlinearscale
        
        self.filename = filename  # csv datafile containing observer responses
        self.scale = None
        self.lscale = None
        self.stim = None
        self.sigma = None
        
        self.mns = None
        self.ci95 = None
        self.sigmamns = None
        self.sigmaci95= None
        
        self.seq=[]  # sequence of commands in R
        self.mldsfile=''
        self.returncode=-1
        
        # initialize commands for execution in R
        self.initcommands()
        
    ###################################################################################################  
    def initcommands(self):
        
        seq1 = ["library(MLDS)\n", 
               "df <- read.table('%s', sep=" ", header=TRUE)\n" % self.filename,
                "stim <- sort(unique(df$s1))\n",
                "results <- data.frame(resp = as.integer(df$Response), S1= match(df$s1, stim), S2=match(df$s2, stim), S3=match(df$s3, stim))\n",
                "attr(results, \"stimulus\") <- stim\n",
                "attr(results, \"invord\") <- as.logical( df$invord )\n",
                "class(results) <- c(\"mlbs.df\", \"data.frame\")\n",
                 "obs.mlds <- mlds(results)\n"]
                 
        # writing perceptual scale calculation
        if self.standardscale:
            seq1a = ["ll<- length(obs.mlds$pscale)\n",
                     "pscale <- obs.mlds$pscale/obs.mlds$pscale[ll]\n",
                     "sigma  <- 1/obs.mlds$pscale[ll]\n"]           
        else: 
            seq1a = ["pscale <- obs.mlds$pscale\n",
                     "sigma  <- obs.mlds$sigma\n"]  
            
        seq1 = seq1 + seq1a
            
        # getting or not linear scale
        if self.getlinearscale or self.boot:
            seq1b = ["df$Response2 <- as.integer(abs(df$s1 - df$s2) < abs(df$s2 - df$s3))\n",
                "results2 <- data.frame(resp = as.integer(df$Response2),  S1= match(df$s1, stim), S2=match(df$s2, stim), S3=match(df$s3, stim))\n",
                "attr(results2, \"stimulus\") <- stim\n",
                "attr(results2, \"invord\") <- as.logical( df$invord )\n",
                "class(results2) <- c(\"mlbs.df\", \"data.frame\")\n",
                "obs.mlds2 <- mlds(results2)\n"]
                
            if self.standardscale:
                seq1c = ["ll<- length(obs.mlds2$pscale)\n",
                         "pscale2 <- obs.mlds2$pscale/obs.mlds2$pscale[ll]\n",
                         "sigma2  <- 1/obs.mlds2$pscale[ll]\n"]
            else: 
                seq1c = ["pscale2 <- obs.mlds2$pscale\n",
                         "sigma2  <- obs.mlds2$sigma\n"]     
                                
            seq1 = seq1 + seq1b + seq1c
        
        # if we want confidence intervals, we need to bootstrap
        if self.boot:
            seq2a = ["obs.bt <- boot.mlds(obs.mlds, 10000)\n"]
            
            if self.standardscale:  # still to add
                seq2b = ["samples <- obs.bt$boot.samp\n"]
            else:
                seq2b = ["n <- nrow(obs.bt$boot.samp)\n",
                         "samples <- apply(obs.bt$boot.samp, 2, function(x) x/x[n])\n"]
            # manually computing mean and sd. from the bootstrap samples. R routine for unconstrained scales does not work properly        
            seq2c = ["obs.mns <- c(0, apply(samples, 1, mean))\n",
                     "obs.sd  <- c(0, apply(samples, 1, sd))\n",
                     "obs.95ci <- qnorm(0.975) * obs.sd\n"
                     "dd <- data.frame( row.names = c( obs.mlds$stimulus, 'sigma'), pscale_obs = c( pscale, sigma), pscale_diff = c( pscale2, sigma2), mns = obs.mns, ci95=obs.95ci)\n"]
            seq2 = seq2a + seq2b + seq2c
        
        # if not, we just save the obtained scale, with or without the linear one
        elif self.getlinearscale:
            seq2 = ["dd <- data.frame( row.names = c( obs.mlds$stimulus, 'sigma'), pscale_obs = c( pscale, sigma), pscale_diff = c( pscale2, sigma2))\n"]
        
        else: 
            seq2 = ["dd <- data.frame( row.names = c( obs.mlds$stimulus, 'sigma'), pscale_obs = c( pscale, sigma))\n"]
        
        # writing file
        self.mldsfile = str(uuid.uuid4()) + '.csv' # csv file with mlds results
        seq3 = ["write.csv(dd, file=\"%s\")\n" % self.mldsfile]
        self.seq = seq1 + seq2 + seq3
    
    
    ###################################################################################################    
    def run(self):
       
        # run mlds analysis on R
        if self.verbose:
            print "executing in R..."     
            
        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in self.seq:
            proc.stdin.write( line )
        proc.communicate()
        
        self.returncode = proc.returncode
        
        self.readresults()
   
    ################################################################################################### 
    def readresults(self):
        
        ## reading results
        if self.returncode==0:
            if self.verbose:
                print "reading MLDS results"          
            data=[]
            csvfile = open(self.mldsfile, 'rb')
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            reader.next()
            for row in reader:
                data.append( np.asarray(row) )
            csvfile.close()
            
            arr=np.asarray(data[:-1], dtype=float)
            
            # parsing to attributes  
            self.stim = arr.T[0]
            self.scale= arr.T[1]
            self.sigma = float(data[-1][1])
            self.status=1
            
            if self.getlinearscale or self.boot:
                self.lscale = arr.T[2]
                
                if self.boot:
                    self.mns = arr.T[3]
                    self.ci95 = arr.T[4]  
                    self.sigmamns = float(data[-1][3])
                    self.sigmaci95= float(data[-1][4])
                    self.status=2
                
        else:
            print "error in execution"      
            self.status= -1
            
        if not self.keepfiles:
            os.remove(self.mldsfile)


        
###############################################################################
########################## utilities for this class  #########################   
     
def plotscale(s, observer="", color='blue', offset=0, linewidth=1, elinewidth=1):
    
    import matplotlib.pyplot as plt

    if s.boot:
        if s.standardscale:
            label = "%s, $\hat{\sigma}=%.3f \pm %.3f$" % (observer, s.sigmamns, s.sigmaci95)
        else:
            label = "%s" % observer
            
        #plt.errorbar(s.stim, s.mns, yerr=s.ci95, label=label)
        stimoff= np.hstack((0, s.stim[1:]+ offset))
        plt.errorbar(stimoff, s.mns, yerr= s.ci95, fmt=None, ecolor= color, elinewidth=elinewidth)
        plt.plot(s.stim, s.mns, color= color, label=label, linewidth=linewidth)
        
    else:
        if s.standardscale:
            label = "%s, $\hat{\sigma}=%.3f$" % (observer, s.sigma)
        else:
            label = "%s" % observer
        
        plt.plot(s.stim, s.scale, color= color, label=label, linewidth=linewidth)

        
###############################################################################
########################## utilities for experiments #########################
## Utilities for MLDS  experiments
# @author: G. Aguilar, Nov 2013

def generate_quadruples(stim):
    """
    Generate quadruples for difference scaling experiment.
    
    Usage:   quads, indx, invord = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: quads:  list of all possible non-overlapping quadruples (stimulus values)
            indx:    (idem)   but indices values
            invord:  0 or 1, if non-inverted or inverted order of stimulus values
    
    Ported from routine "runQuadExperiment.R" from the R package "MLDS"
    
    Reference: 
        Knoblauch, K. and Maloney, L. T. (2008) MLDS: Maximum likelihood difference scaling in R.
        Journal of Statistical Software, 25:2, 1–26, http://www.jstatsoft.org/v25/i02.

    """
    numstim = len(stim)

    # getting all possible combinations in a list, C(n,k)
    allTrials = [list(s) for s in itertools.combinations(range(numstim), 4)]
    
    # randomize order of trials
    random.shuffle(allTrials)
    
    # randomize order or pairs 
    topbot = np.random.binomial(1, 0.5, len(allTrials) )
    
    for t in range(len(allTrials)):
        if topbot[t]:
            allTrials[t] = [ allTrials[t][2], allTrials[t][3], allTrials[t][0], allTrials[t][1] ]
            
    # create list of stimulus quadruples
    quadruples = [stim[t]   for t in allTrials]

    # returns the quadruples, and the indices
    return (quadruples, allTrials, topbot)
    
    
    
def generate_triads(stim):
    """
    Generate triads for difference scaling experiment.
    
    Usage:   traids, indx, invord = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: triads:  list of all possible non-overlapping triads (stimulus values)
            indx:    (idem)   but indices values
            invord:  0 or 1, if non-inverted or inverted order of stimulus values


    Ported from routine "runTriadExperiment.R" from the R package "MLDS"
    
    Reference: 
        Knoblauch, K. and Maloney, L. T. (2008) MLDS: Maximum likelihood difference scaling in R.
        Journal of Statistical Software, 25:2, 1–26, http://www.jstatsoft.org/v25/i02.
        
    """
    numstim = len(stim)

    # getting all possible combinations in a list, C(n,k)
    allTrials = [list(s) for s in itertools.combinations(range(numstim), 3)]
    
    # randomize order of trials
    random.shuffle(allTrials)
    
    # randomize order or pairs 
    topbot = np.random.binomial(1, 0.5, len(allTrials) )
    
    for t in range(len(allTrials)):
        if topbot[t]:
            allTrials[t] = [ allTrials[t][2], allTrials[t][1], allTrials[t][0] ]
            # careful here to make non-overlapping triads.. when inverting, the stim 1 must stay in the middle.
     
    # create list of stimulus triads
    triads = [stim[t]   for t in allTrials]

    # returns the triads, and the indices
    return (triads, allTrials, topbot)
 

###############################################################################
############################# Testing routine ##################################

if __name__ == "__main__":
    
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    
    ## running a test example
    # 'test.csv' was created simulating a power function, with exponent=2, 5 blocks and noise parameter sigma=0.1
    
    # Case 1: Simple scale, unconstrained 
    obs = MLDSObject( 'test.csv', boot=False, standardscale=False) # initialize object
    obs.run()                                                   # runs Mlds in R
     
    fig=plt.figure()
    plt.plot(obs.stim, obs.scale)
    plt.xlabel('Stimulus')
    plt.ylabel('Difference scale')
    plt.title('$\sigma = %.2f$' % obs.sigma)
    #plt.show()
    fig.savefig('Fig1.png')

    
    
    # Case 2: Standard scale
    obs = MLDSObject( 'test.csv', boot=False, standardscale=True) 
    obs.run()                                                   

    fig=plt.figure()
    plt.plot(obs.stim, obs.scale)
    plt.xlabel('Stimulus')
    plt.ylabel('Difference scale')
    plt.title('$\sigma = %.2f$' % obs.sigma)
    #plt.show()
    fig.savefig('Fig2.png')

    
    
    # Case 3: Standard scale and bootstrap
    obs = MLDSObject( 'test.csv', boot=True, standardscale=True) 
    obs.run()                                                   

    fig=plt.figure()
    plt.errorbar(obs.stim, obs.mns, yerr=obs.ci95)
    plt.xlabel('Stimulus')
    plt.ylabel('Difference scale')
    plt.title('$\sigma = %.2f \pm %.2f$' % (obs.sigmamns, obs.sigmaci95))
    plt.xlim(0, 1.05)
    #plt.show()
    fig.savefig('Fig3.png')

    
    
    # Case 4: Unconstrained scale and bootstrap
    obs = MLDSObject( 'test.csv', boot=True, standardscale=False) 
    obs.run()                                                  

    fig=plt.figure()
    plt.errorbar(obs.stim, obs.mns, yerr=obs.ci95)
    plt.xlabel('Stimulus')
    plt.ylabel('Difference scale')
    plt.title('$\sigma = %.2f \pm %.2f$' % (obs.sigmamns, obs.sigmaci95))
    plt.xlim(0, 1.05)
    #plt.show()
    fig.savefig('Fig4.png')
    
    
#EOF    
