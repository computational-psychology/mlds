# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 14:27:11 2014

@author: G. Aguilar, Apr 2014
"""

import numpy as np
import os
import csv
import itertools
import random
import subprocess


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
    
    obs.do()  --> run MLDS analysis


    Attributes
    ----------
    
    obs.stim  --> stimulus scale\n
    obs.scale --> perceptual scale\n
    obs.lscale --> linear scale (assuming perfect differentiating)\n
    
    when boot = True, \n
    obs.mns and obs.ci95  --> mean and 95% confidence interval of estimates after bootstrap.         
    
    """
    def __init__(self, filename, boot=False, keepfiles=False, standardscale=True, getlinearscale=False):
        
        self.status=0 # 0: just initialized, 1: mlds computed, 2: bootstrapped, -1: error
        
        # flags
        self.boot = boot
        self.keepfiles = keepfiles
        self.standardscale = standardscale
        self.getlinearscale = getlinearscale
        
        self.filename = filename  # csv datafile containing observer responses
        self.scale = None
        self.lscale = None
        self.stim = None
        
        self.mns = None
        self.ci95 = None
    
    def run(self):
        
        # write R file         
        Rfile    = "tmp.R"         # r script
        mldsfile = "tmp.csv"       # csv file with mlds results
        #fid = open(Rfile, "w+")
        
        seq1 = ["library(MLDS)\n", 
               "df <- read.table('%s', sep=" ", header=TRUE)\n" % self.filename,
                "stim <- seq(from=0, to=60, by=5)\n",  # this must be changed to accomodate other stimulus configurations
                "results <- data.frame(resp = as.integer(df$Response), S1= df$i1, S2=df$i2, S3=df$i3)\n",
                "attr(results, \"stimulus\") <- stim\n",
                "attr(results, \"invord\") <- as.logical( df$invord )\n",
                "class(results) <- c(\"mlbs.df\", \"data.frame\")\n",
                 "obs.mlds <- mlds(results)\n"]
                 
        # writing perceptual scale calculation
        if self.standardscale:
            seq1a = ["ll<- length(obs.mlds$pscale)\n",
                     "pscale <- obs.mlds$pscale/obs.mlds$pscale[ll]\n"]           
        else: 
            seq1a = ["pscale <- obs.mlds$pscale\n"]
            
        seq1 = seq1 + seq1a
            
        # getting or not linear scale
        if self.getlinearscale or self.boot:
            seq1b = ["df$Response2 <- as.integer(abs(df$s1 - df$s2) < abs(df$s2 - df$s3))\n",
                "results2 <- data.frame(resp = as.integer(df$Response2), S1= df$i1, S2=df$i2, S3=df$i3)\n",
                "attr(results2, \"stimulus\") <- stim\n",
                "attr(results2, \"invord\") <- as.logical( df$invord )\n",
                "class(results2) <- c(\"mlbs.df\", \"data.frame\")\n",
                "obs.mlds2 <- mlds(results2)\n"]
                
            if self.standardscale:
                seq1c = ["ll<- length(obs.mlds2$pscale)\n",
                         "pscale2 <- obs.mlds2$pscale/obs.mlds2$pscale[ll]\n"]
            else: 
                seq1c = ["pscale2 <- obs.mlds2$pscale\n"]     
                                
            seq1 = seq1 + seq1b + seq1c
        
        # if we want confidence intervals, we need to bootstrap
        if self.boot:
            seq2a = ["obs.bt <- boot.mlds(obs.mlds, 10000)\n"]
            
            if self.standardscale:
                seq2b = ["obs.bt.res <- summary(obs.bt)\n",
                        "obs.mns <- obs.bt.res[,1]\n",
                        "obs.95ci <- qnorm(0.975) * obs.bt.res[,2]\n"]
            else:
                seq2b = ["n <- nrow(obs.bt$boot.samp)\n",
                         "obs.mns <- c(0, obs.bt$bt.mean[-n] / obs.bt$bt.mean[n])\n",
                        "obs.sd <- c(0, obs.bt$bt.sd[-n] / obs.bt$bt.mean[n])\n",
                        "obs.95ci <- qnorm(0.975) * obs.sd\n"]
            seq2c = ["dd <- data.frame( row.names = obs.mlds$stimulus, pscale_obs = pscale, pscale_diff = pscale2, mns = obs.mns, ci95=obs.95ci)\n"]
            seq2 = seq2a + seq2b + seq2c
        
        # if not, we just save the obtained scale, with or without the linear one
        elif self.getlinearscale:
            seq2 = ["dd <- data.frame( row.names = obs.mlds$stimulus, pscale_obs = pscale, pscale_diff = pscale2)\n"]
        
        else: 
            seq2 = ["dd <- data.frame( row.names = obs.mlds$stimulus, pscale_obs = pscale)\n"]
        
        # writing file
        seq3 = ["write.csv(dd, file=\"%s\")\n" % mldsfile]
        seq = seq1 + seq2 + seq3
        #fid.writelines(seq)
        #fid.close()  
                  
        # run mlds analysis on R
        print "executing in R..."     
        # option 1
        #st = os.system("R --no-save < %s" % Rfile)
        # option 2
        #proc = subprocess.Popen("R --no-save < %s" % Rfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #(out,err) = proc.communicate()
        # option 3 
        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in seq:
            proc.stdin.write( line )        
        proc.communicate()
    
        ## reading results
        if proc.returncode==0:
            print "reading MLDS results"          
            data=[]
            with open(mldsfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                header = reader.next()
                for row in reader:
                    data.append( np.asarray(row, dtype=float) )
                csvfile.close()
                arr = np.array(data).T    
                
                # parsing to attributes  
                self.stim = arr[0]
                self.scale= arr[1]
                self.status=1
                
                if self.getlinearscale or self.boot:
                    self.lscale = arr[2]
                    
                    if self.boot:
                        self.mns = arr[3]
                        self.ci95 = arr[4]           
                        self.status=2
                
        else:
            print "error in execution"      
            self.status= -1
            
        if not self.keepfiles:
            os.remove(mldsfile)
            #os.remove(Rfile)
        
        
###############################################################################
########################## utilities for experiments #########################
## Utilities for MLDS  experiments
# @author: G. Aguilar, Nov 2013

def generate_quadruples(stim):
    """
    Generate quadruples for difference scaling experiment.
    
    Usage:   quads, indx = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: quads:  list of all possible non-overlapping quadruples (stimulus values)
            indx:    (idem)   but indices values
    
    Ported from routines runQuadExperiment.R   from R package "MLDS"
    
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
    
    Usage:   traids, indx = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: triads:  list of all possible non-overlapping triads (stimulus values)
            indx:    (idem)   but indices values


    Ported from routines runTriadExperiment.R   from R package "MLDS"
    
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
 
      
#EOF            
    
