# -*- coding: utf-8 -*-
"""
Utilities for MLDS estimation, interfacing with R.

@author: G. Aguilar, Nov 2013, rev. Apr 2014.
rev. Aug 2014: diagnostics added, saving to R datafile.

TO ADD: save a MLDSObject to a file. Better to save it with pickle, and
do it loading the csvfile to a data array or so. So the object can be run 
from csvfile or from data array. This would imply, if data is non-empty, to 
save the data on a csvfile and that is passed to R.

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
    MLDS Object allows the quick analysis of perceptual scales from triad experiments
    using the maximum-likelihood difference scaling method (Maloney & Yang, 2003). 
    
    It runs the analysis on R, using the MLDS package.
    It requires, therefore, a running R version with packages MLDS and psyphy installed.
    
    Usage
    ----------
    
    obs = MLDSObject( datafile, *flags  )
   
    boolean flags: 
    
    boot -> calculate bootstrap estimates (default : False) \n
    keepfiles -> keep temporary R and csv files (default: False) \n
    standardscale -> calculate perceptual scale in range [0,1] (default: True)\n
    getlinearscale -> calculate linear scale as if differencing perfecly (default=False)\n
    save -> save results into an R datafile\n
   
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
    
    def __init__(self, filename, boot=False, keepfiles=False, standardscale=True, getlinearscale=False, verbose=False, save=True):
        
        self.status=0 # 0: just initialized, 1: mlds computed, 2: bootstrapped, -1: error
        self.verbose = verbose
        
        # flags
        self.boot = boot
        self.keepfiles = keepfiles
        self.standardscale = standardscale
        self.getlinearscale = getlinearscale
        self.saveRobj = save
        
        self.linktype ="probit"
        self.linkgam  = 0.0
        self.linklam  = 0.0
        
        self.filename = filename  # csv datafile containing observer responses
        
        self.getRdatafilename()
         
        # scale and noise param, stimulus vector
        self.scale = None
        self.lscale = None
        self.stim = None
        self.sigma = None
        
        # bootstrap 
        self.mns = None
        self.ci95 = None
        self.sigmamns = None
        self.sigmaci95= None
        
        # diagnostic measures
        self.AIC = None
        self.DAF = None
        self.prob = None
                
        
        self.seq=[]  # sequence of commands in R
        self.mldsfile=''
        self.returncode=-1
        
        # parallel execution of bootsrap
        self.parallel= False
        self.workers = ['"localhost"', '"localhost"']
        self.master = '"localhost"'
        
        # initialize commands for execution in R
        self.initcommands()
    
    def getRdatafilename(self, force_refit=False):
        rootname = self.filename.split('.')[0]
        
        if self.linkgam==0.0 and self.linklam==0.0:
            self.Rdatafile = rootname + '_' + self.linktype + '.MLDS'
            
        else:
            self.Rdatafile = rootname + '_' + self.linktype + '_refit' + '.MLDS'
            
        if force_refit:
            self.Rdatafile = rootname + '_' + self.linktype + '_refit' + '.MLDS'
       
    ###################################################################################################  
    def initcommands(self):
        
        self.getRdatafilename()
        
        seq1 = ["library(MLDS)\n",
                "library(psyphy)\n",
               "df <- read.table('%s', sep=" ", header=TRUE)\n" % self.filename,
                "stim <- sort(unique(df$s1))\n",
                "results <- data.frame(resp = as.integer(df$Response), S1= match(df$s1, stim), S2=match(df$s2, stim), S3=match(df$s3, stim))\n",
                "attr(results, \"stimulus\") <- stim\n",
                "attr(results, \"invord\") <- as.logical( df$invord )\n",
                "class(results) <- c(\"mlbs.df\", \"data.frame\")\n",
                 "obs.mlds <- mlds(results, lnk=%s.2asym(g = %f, lam = %f))\n" % (self.linktype, self.linkgam, self.linklam) ]
                 
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
            if self.parallel:
                seq2a = ["library(snow)\n",
                         "source('~/git/slantfromtex/mlds/pboot.mlds.R')\n",
                         "workers <- c(%s)\n" % ",".join(self.workers),
                         "master <- %s\n" % self.master,
                         "obs.bt <- pboot.mlds(obs.mlds, 10000, workers = workers, master=master )\n"]
            else:
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
        
        # saving R objects into datafile
        if self.boot and self.saveRobj:
            seq4 = ["save(results, obs.mlds, obs.bt, obs.mns, obs.95ci, file='%s')\n" % self.Rdatafile]
        elif self.saveRobj:
            seq4 = ["save(results, obs.mlds, file='%s')\n" % self.Rdatafile]
        else:
            seq4 = ["\n"]
            
        self.seq = seq1 + seq2 + seq3 + seq4
    
    
    ###################################################################################################    
    def run(self):
        
        self.initcommands()
       
        # run MLDS analysis in R
        if self.verbose:
            print "executing in R..."     
            
        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in self.seq:
            proc.stdin.write( line )
        (out, err) = proc.communicate()
        
        self.returncode = proc.returncode
                
        if self.verbose:
            print out                    
        
        if self.returncode==0:
            self.readresults()
        else:
            print err
            raise RuntimeError("Error in execution within R (see error output above)")
        
        
            
   
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
            
    #################################################################### 
    def rundiagnostics(self): 
        
        self.getRdatafilename()
        
        seqdiag = ["library(MLDS)\n", 
        "load('%s')\n" % self.Rdatafile, 
        "library(snow)\n",
        "source('~/git/slantfromtex/mlds/pbinom.diagnostics.R')\n",
        "workers <- c(%s)\n" % ",".join(self.workers),
        "master <- %s\n" % self.master,
        "obs.diag.prob <- pbinom.diagnostics (obs.mlds, 10000, workers=workers, master=master)\n"
        "save(results, obs.mlds, obs.bt, obs.mns, obs.95ci, obs.diag.prob, file='%s')\n" % self.Rdatafile ]
        
       
        # run MLDS analysis in R
        if self.verbose:
            print "executing in R..."     
            
        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in seqdiag:
            proc.stdin.write( line )
        (out, err) = proc.communicate()
        
        self.returncode = proc.returncode
                
        if self.verbose:
            print out                    
        
        if self.returncode==0:
            self.readdiags()
        else:
            print err
            raise RuntimeError("Error in execution within R (see error output above)")
            
            
    #####
    def readdiags(self):
        
        import rpy2.robjects as robjects

        ### loading file
        objs = robjects.r['load']("%s" % self.Rdatafile)
        objl = list(objs)
        
        # loading R objects into python variables
        obsmlds = robjects.r['obs.mlds']
                
        ### Akaike information criterion
        self.AIC = list(robjects.r['AIC'](obsmlds))[0]
        
        #### DAF
        # function definition
        robjects.r('''
                DAF <- function(g) {
                    (g$obj$null.deviance - deviance(g$obj)) / g$obj$null.deviance  
                }
                ''')
        daf = robjects.globalenv['DAF']        
        
        self.DAF = list( daf(obsmlds) )[0]
        
        # prob
        if 'obs.diag.prob' in objl:
            diagprob = robjects.r['obs.diag.prob']
            self.prob = list( diagprob[4])[0]
        else:
            print "bootstrap diagnostics are not yet calculated"


    ##########################################################################    
    def estgamlam(self):
        
        gamlamfile = str(uuid.uuid4()) + '.csv' 

        # as in the book chapter
        seq = ["library(MLDS)\n",
               "library(psyphy)\n",
               "load('%s')\n" % self.Rdatafile,
               "resp <- obs.mlds$obj$data$resp\n",
               "ps <- psyfun.2asym( cbind(resp, 1-resp) ~ . -1 , data= obs.mlds$obj$data, link= %s.2asym)\n" % self.linktype,
               "c(ps$lambda, plogis(qlogis(ps$lambda) + c(-ps$SElambda, ps$SElambda)))\n",   # as in psyfun.2asym source code
               "c(ps$gam, plogis(qlogis(ps$gam) + c(-ps$SEgam, ps$SEgam)))\n",
               "dd <- data.frame( gamma=c(ps$gam, plogis(qlogis(ps$gam) + c(-ps$SEgam, ps$SEgam))),   lambda = c(ps$lambda, plogis(qlogis(ps$lambda) + c(-ps$SElambda, ps$SElambda))) )\n",
               "write.csv(dd, file=\"%s\")\n" % gamlamfile]
        
        
        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
        for line in seq:
            proc.stdin.write( line )
        (out, err) = proc.communicate()
        
        if self.verbose:
            print out                    
        
        print self.returncode
        
        try:
            self.readgamlam( gamlamfile )
        except:
            print "Error"
            print out
            print err
            raise RuntimeError("Error in execution within R (see error output above)")
            
            
    ##########################################################################            
    def readgamlam(self, gamlamfile):
        
        def tofloat(s):
            try:
                return round ( float(s), 5 )
            except ValueError:
                return np.nan
            
        
        data=[]
        csvfile = open(gamlamfile, 'rb')
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        reader.next()
        for row in reader:
            data.append( np.asarray(row) )
        csvfile.close()
        
                
        self.gam = tofloat(data[0][1])
        self.lam = tofloat(data[0][2])
        
        
        if self.gam > 0.0:
            self.SEgam = np.array([tofloat(data[1][1]), tofloat(data[2][1])])
        else:
            self.SEgam = None
            
        if self.lam > 0.0:
            self.SElam = np.array([tofloat(data[1][2]), tofloat(data[2][2])])
        else:
            self.SElam = None
        
        os.remove(gamlamfile)
        
    ##########################################################################                
    def load(self, force_refit=False):
        
        self.getRdatafilename( force_refit)
               
        # scale and bootstrap 
        self.readobjectresults()
        
        # diagnostic measures
        self.readdiags()
        
        
        
    ##########################################################################    
    def readobjectresults(self):
        
        import rpy2.robjects as robjects

        ### loading file
        robjects.r['load']("%s" % self.Rdatafile)

        
        # loading R objects into python variables
        obsmlds = robjects.r['obs.mlds']        
        
        if self.standardscale:
            self.sigma = 1.0/ np.array(obsmlds[0])[-1]
            self.scale = np.array(obsmlds[0]) / np.array(obsmlds[0])[-1]
        
        else:
            self.sigma = list(obsmlds[2])[0]
            self.scale = np.array(obsmlds[0])
        
        self.stim = np.array(obsmlds[1])
        self.status=1
        
        # if bootstrapped
        if self.boot:
            obsmns  = robjects.r['obs.mns']
            obs95ci  = robjects.r['obs.95ci']
            
            self.mns = np.array(obsmns)[:-1]
            self.ci95 =  np.array(obs95ci)[:-1]
            
            self.sigmamns = np.array(obsmns )[-1]
            self.sigmaci95= np.array(obs95ci)[-1]
            self.status=2
        
   
    
        
###############################################################################
########################## utilities for this class  #########################      
    
def plotscale(s, observer="", color='blue', offset=0, linewidth=1, elinewidth=1):
    
    import matplotlib.pyplot as plt

    if s.boot:
        if s.standardscale:
            label = "%s, $\hat{\sigma}=%.3f \, [%.3f, \, %.3f]$" % (observer, s.sigmamns, s.sigmamns - s.sigmaci95, s.sigmamns + s.sigmaci95)
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
    import time
    
    simplecases=False
    bootstrap=False
    benchmark=True
    
    
    
    ## running a test example
    # 'test.csv' was created simulating a power function, with exponent=2, 5 blocks and noise parameter sigma=0.1
    
    if simplecases:
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

    
    if bootstrap:
        
        # Case 3: Standard scale and bootstrap
        obs = MLDSObject( 'test.csv', boot=True, standardscale=True) 
        obs.run()    
    
    
        fig=plt.figure()
        plt.errorbar(obs.stim, obs.mns, yerr=obs.ci95)
        plt.xlabel('Stimulus')
        plt.ylabel('Difference scale')
        plt.title('$\sigma = %.2f \pm %.2f$' % (obs.sigmamns, obs.sigmaci95))
        plt.xlim(0, 1.05)
        plt.show()
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
        plt.show()
        fig.savefig('Fig4.png')
    
    
    # Benchmark parallelization
    if benchmark:
        
        obs = MLDSObject( 'test.csv', boot=True, standardscale=True)
        
        obs.parallel= False
        
        startime=time.time()
        obs.run()
        print 'without parallel: %d sec ' % (time.time() - startime)
        
        
        # parallel but in duo core computer
        obs.parallel= True                           
        
        startime=time.time()
        obs.run()
        print 'with parallel and 2 CPUs: %d sec ' % (time.time() - startime)
        
        
        # parallel but in quad computer
        obs.workers=['"localhost"'] *4                          
        
        startime=time.time()
        obs.run()
        print 'with parallel and 4 CPUs: %d sec' % (time.time() - startime)
        
        
        # in the lab with 8 CPUS
        obs.workers=['"localhost"', '"localhost"', '"localhost"', '"localhost"',  '"130.149.57.124"' , '"130.149.57.124"', '"130.149.57.124"', '"130.149.57.124"']
        obs.master = '"130.149.57.105"'            
        
        startime=time.time()
        obs.run()
        print 'with parallel and 8 CPUs: %d sec' % (time.time() - startime)

     
#EOF    
