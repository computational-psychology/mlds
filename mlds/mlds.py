# -*- coding: utf-8 -*-
"""
Utilities for MLDS estimation, interfacing with R.

@author: G. Aguilar, Nov 2013, rev. Apr 2014.
rev. Aug 2014: diagnostics added, saving to R datafile.
rev. Jan 2015: bootstrap samples can be accessed.
rev. Feb 2015: adds plot diagnostics, as in R.
"""

import numpy as np
import os, sys
import csv
import itertools
import random
import subprocess
import uuid
import multiprocessing


def tofloat(s):
    try:
        return round(float(s), 5)
    except ValueError:
        return np.nan

class MLDSObject:
    """
    MLDS Object allows the quick analysis of perceptual scales from triad
    experiments using the maximum-likelihood difference scaling method
    (Maloney & Yang, 2003).

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
    obs.mns and obs.ci95  --> mean and 95% confidence interval for
                                scale estimates after bootstrap.
    obs.sigmamns and obs.sigmaci95 -->   mean and 95% confidence interval for
                                sigma estimate after bootstrap.

    """

    def __init__(self, filename, boot=False, dimension_unit = False, keepfiles=False, standardscale=True, getlinearscale=False, verbose=False, save=True):

        # 0: just initialized, 1: mlds computed, 2: bootstrapped, -1: error
        self.status = 0
        self.verbose = verbose

        # flags
        self.boot = boot
        self.nsamples = 10000
        if not dimension_unit:
            self.dimension_unit = ''
        else:
            self.dimension_unit = dimension_unit
            
        self.keepfiles = keepfiles
        self.standardscale = standardscale
        self.getlinearscale = getlinearscale  # I will deprecate this
        self.saveRobj = save

        self.linktype = "probit"
        self.linkgam = 0.0
        self.linklam = 0.0

        self.filename = filename  # csv datafile containing observer responses
        self.rootname = self.filename.split('.')[0]

        self.getRdatafilename()
        
        # default column names in the text file to be read
        self.colnames = {'stim' : ['s1', 's2', 's3'], 
                         'response': 'Response'}
        
        
        # scale and noise param, stimulus vector
        self.scale = None
        self.lscale = None  # deprecating this
        self.stim = None
        self.sigma = None

        # bootstrap
        self.correctedCI = False
        self.mns = None
        self.ci95 = None
        self.sd = None
        self.sigmamns = None
        self.sigmaci95 = None
        self.scalesbt = None
        self.sigmabt = None

        # diagnostic measures
        self.AIC = None
        self.DAF = None
        self.prob = None
        self.pmc = None
        self.residuals = None

        ##
        self.seq = []  # sequence of commands in R
        self.mldsfile = ''
        self.returncode = -1

        # parallel execution of bootsrap
        ncpus = multiprocessing.cpu_count()
        if ncpus == 1:
            self.parallel= False
        else:
            self.parallel = True
            self.workers = ['"localhost"']* ncpus
            self.master = '"localhost"'

        # initialize commands for execution in R
        self.initcommands()


    def printinfo(self):
        print "MLDSObject based on file %s " % self.filename
        if self.status==0:
            print "not yet run"
        elif self.status==1 or self.status==2:
            print "   stimulus: ", self.stim
            print "   scale: ", self.scale
            print "   sigma: ", self.sigma

        if self.status==2:
            print "   ci 2.5%: ", self.ci95[0]
            print "   ci 97.5%: ", self.ci95[1]
            if self.correctedCI:
                print "corrected CIs"

    def getRdatafilename(self, force_refit=False):
        
        s = []
        s.append(self.rootname)
        if not self.dimension_unit=='':
            s.append(self.dimension_unit)

        if self.standardscale:
            s.append('norm')
        else:
            s.append('unnorm')

        s.append(self.linktype)
        
        if not (self.linkgam == 0.0 and self.linklam == 0.0) or force_refit:
            s.append('refit')

        self.Rdatafile = "_".join(s) + '.MLDS'
         

    ###################################################################################################
    def initcommands(self):
        """
        Generate list of commands to be executed in R
        
        """
        self.getRdatafilename()

        seq = ["library(MLDS)\n",
                "library(psyphy)\n",
                "d.df <- read.table('%s', sep=" ", header=TRUE)\n" % self.filename,
                "stim <- sort(unique(c(d.df$%s, d.df$%s, d.df$%s)))\n" % (self.colnames['stim'][0], self.colnames['stim'][1], self.colnames['stim'][2]),
                "results <- with(d.df, data.frame(resp = %s, S1= match(%s, stim), S2=match(%s, stim), S3=match(%s, stim)))\n" % (self.colnames['response'], self.colnames['stim'][0], self.colnames['stim'][1], self.colnames['stim'][2]),
                "results <- as.mlbs.df(results, st=stim)\n"]

        # still to be debug is the fact that passing a link function with gamma = lamda = 0 gives different results as passing the link function as a word.
        if self.linkgam == 0.0 and self.linklam == 0:
            seq.append("obs.mlds <- mlds(results, lnk='%s')\n" % self.linktype)
        else:
            seq.append("obs.mlds <- mlds(results, lnk=%s.2asym(g = %f, lam = %f))\n" % (self.linktype, self.linkgam, self.linklam))

        # writing perceptual scale calculation
        if self.standardscale:
            seq.extend(["ll<- length(obs.mlds$pscale)\n",
                     "pscale <- obs.mlds$pscale/obs.mlds$pscale[ll]\n",
                     "sigma  <- 1/obs.mlds$pscale[ll]\n"])
        else:
            seq.extend(["pscale <- obs.mlds$pscale\n",
                     "sigma  <- obs.mlds$sigma\n"])


        # if we want confidence intervals, we need to bootstrap
        if self.boot:
            if self.parallel:
                seq.extend(["library(snow)\n",
                         "source(paste('%s', '/pboot.mlds.R', sep=''))\n" % os.path.dirname(sys.modules[__name__].__file__),
                         "workers <- c(%s)\n" % ",".join(self.workers),
                         "master <- %s\n" % self.master,
                         "obs.bt <- pboot.mlds(obs.mlds, %d, workers = workers, master=master )\n" % self.nsamples])
            else:
                seq.extend(["obs.bt <- boot.mlds(obs.mlds, %d)\n" % self.nsamples])

            if self.standardscale:  # still to add
                seq.extend(["samples <- obs.bt$boot.samp\n"])
            else:
                seq.extend(["n <- nrow(obs.bt$boot.samp)\n",
                         "samples <- apply(obs.bt$boot.samp, 2, function(x) x/x[n])\n"])

            # manually computing mean and sd. from the bootstrap samples. R routine for unconstrained scales does not work properly
            seq.extend(['obs.mns <- c("0" = 0, apply(samples, 1, mean))\n',
                     'obs.sd  <- c("0"= 0, apply(samples, 1, sd))\n',
                     'obs.low <- c("0"=0, apply(samples, 1, quantile, probs = 0.025))\n',
                     'obs.high <- c("0"=0, apply(samples, 1, quantile, probs = 0.975))\n'])

            # Efron's corrected CI, suggested by Ken. Check this in the book, it reverts the low and high CI order.
            if self.correctedCI:
                seq.extend(['obs.low <- 2*obs.mns - obs.high\n',
                            'obs.high <- 2*obs.mns - obs.low\n'])

            seq.append('dd <- data.frame( row.names = c( obs.mlds$stimulus, "sigma"), pscale_obs = c( pscale, sigma), mns = obs.mns, low=obs.low, high=obs.high)\n')

        else:
            seq.extend(["dd <- data.frame( row.names = c( obs.mlds$stimulus, 'sigma'), pscale_obs = c( pscale, sigma))\n"])

        # writing file
        self.mldsfile = str(uuid.uuid4()) + '.csv'  # csv file with mldsresults
        seq.append("write.csv(dd, file=\"%s\")\n" % self.mldsfile)

        # saving R objects into datafile
        if self.boot and self.saveRobj:
            seq.append("save(results, obs.mlds, obs.bt, obs.mns, obs.low, obs.high, obs.sd, samples, file='%s')\n" % self.Rdatafile)
        elif self.saveRobj:
            seq.append("save(results, obs.mlds, file='%s')\n" % self.Rdatafile)
        else:
            seq.append("\n")

        self.seq = seq


    ###################################################################################################
    def run(self):

        self.initcommands()

        # run MLDS analysis in R
        if self.verbose:
            print "executing in R..."

        proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in self.seq:
            proc.stdin.write(line)
        (out, err) = proc.communicate()

        self.returncode = proc.returncode

        if self.verbose:
            print "return code: %d" % self.returncode
            print "output: "
            print out

        if self.returncode==0:
            self.readresults()
        else:
            print err
            raise RuntimeError("Error in execution within R (see error output above)")




    ###################################################################################################
    def readresults(self):

        ## reading results
        if self.returncode == 0:
            if self.verbose:
                print "reading MLDS results"
            data = []
            csvfile = open(self.mldsfile, 'rb')
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            reader.next()
            for row in reader:
                data.append(np.asarray(row))
            csvfile.close()

            arr = np.asarray(data[:-1], dtype=float)

            # parsing to attributes
            self.stim = arr.T[0]
            self.scale = arr.T[1]
            self.sigma = float(data[-1][1])
            self.status = 1

            if self.boot:

                self.mns = arr.T[2]
                self.ci95 = arr.T[3:]
                self.sigmamns = float(data[-1][2])
                self.sigmaci95 = np.array(data[-1][3:], dtype=float)
                self.status = 2

        else:
            print "error in execution"
            self.status = -1

        if not self.keepfiles:
            os.remove(self.mldsfile)

    ###################################################################################################
    def rundiagnostics(self, saveresiduals=False):

        import rpy2.robjects as robjects

        ### loading file
        objs = robjects.r['load']("%s" % self.Rdatafile)
        objl = list(objs)

        if 'obs.diag.prob' in objl:
            self.readdiags()

        else:

            self.getRdatafilename()

            seqdiag = ["library(MLDS)\n",
            "load('%s')\n" % self.Rdatafile]
            
            
            if self.parallel:
                seqdiag.extend(["library(snow)\n",
                "source(paste('%s', '/pbinom.diagnostics.R', sep=''))\n" % os.path.dirname(sys.modules[__name__].__file__),
                "workers <- c(%s)\n" % ",".join(self.workers),
                "master <- %s\n" % self.master,
                "obs.diag.prob <- pbinom.diagnostics (obs.mlds, %d, workers=workers, master=master)\n" % self.nsamples])
                
            else:
                seqdiag.extend(["obs.diag.prob <- binom.diagnostics (obs.mlds, %d)\n" % self.nsamples])
            
            
            if not saveresiduals:
                # I take 95% CI envelopes and necesary variables to plot GoF,
                # and I get rid of all residuals data that is very heavy.
                seqdiag.extend(["alpha <- 0.025\n",
                                "nsim <- dim(obs.diag.prob$resid)[1]\n",
                                "n <- dim(obs.diag.prob$resid)[2]\n",
                                "obs.diag.prob$lowc <- obs.diag.prob$resid[alpha * nsim, ]\n",
                                "obs.diag.prob$highc <- obs.diag.prob$resid[(1 - alpha) * nsim, ]\n",
                                "obs.diag.prob$resid <- NULL\n",
                                "obs.diag.prob$alpha <- alpha\n",
                                "obs.diag.prob$nsim <- nsim\n",
                                "obs.diag.prob$n <- n\n"])
            seqdiag.append('save(results, obs.mlds, obs.bt, obs.mns, obs.low, obs.high, obs.sd, samples, obs.diag.prob, file="%s")\n' % self.Rdatafile)
           
            

            # run MLDS analysis in R
            if self.verbose:
                print "executing in R..."

            proc = subprocess.Popen(["R", "--no-save"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for line in seqdiag:
                proc.stdin.write(line)
            (out, err) = proc.communicate()

            self.returncode = proc.returncode

            if self.verbose:
                print out

            if self.returncode == 0:
                self.readdiags()
            else:
                print err
                raise RuntimeError("Error in execution within R (see error output above)")


    ###################################################################################################
    def readdiags(self):

        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr

        ### loading library
        rmlds = importr("MLDS")

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

        self.DAF = list(daf(obsmlds))[0]

        ### pmc
        self.pmc = robjects.r['pmc'](obsmlds)[0]

        ### residuals
        obj = obsmlds[obsmlds.names.index('obj')]
        self.residuals = np.array(obj[obj.names.index('residuals')])


        #### prob
        if 'obs.diag.prob' in objl:
            self.diagnostics = robjects.r['obs.diag.prob']
            self.prob = list(self.diagnostics[self.diagnostics.names.index('p')])[0]

        #else:
        #   print "bootstrap diagnostics are not yet calculated"

    ###################################################################################################
    def plotdiags(self, width=10, height=5):

        if self.diagnostics is not None:

            import matplotlib.pyplot as plt
            import matplotlib.ticker as ticker


            NumRuns = np.array(self.diagnostics[self.diagnostics.names.index('NumRuns')])
            #resid = np.array(self.diagnostics[1])
            Obs_resid = np.array(self.diagnostics[self.diagnostics.names.index('Obs.resid')])
            ObsRuns = np.array(self.diagnostics[self.diagnostics.names.index('ObsRuns')])[0]

            #nsim = resid.shape[0]
            #n = resid.shape[1]
            n = int(self.diagnostics[self.diagnostics.names.index('n')][0])
            
            lowc = self.diagnostics[self.diagnostics.names.index('lowc')]
            highc = self.diagnostics[self.diagnostics.names.index('highc')]
            
            cdfx = (np.arange(1, n + 1, 1) - 0.5) / n

            fig = plt.figure(figsize=(width, height))
            ax = plt.subplot(1, 2, 1)
            ax.set_xlabel("Deviance residual")
            ax.set_ylabel("Cumulative density function")
            ax.plot(np.sort(Obs_resid), cdfx, 'o', markersize=1,
                     markeredgecolor='k', markerfacecolor='k')
            ax.plot(lowc, cdfx, '-', color='#4C72B0', linewidth=1)
            ax.plot(highc, cdfx, '-', color='#4C72B0', linewidth=1)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.tick_params(right=False, top=False)
            plt.locator_params(axis = 'y', nbins=4)

            ax = plt.subplot(1, 2, 2)
            ax.hist(NumRuns, bins=25, normed=True)
            ax.axvline(ObsRuns, color='k', linewidth=1)
            ax.set_xlabel("Number of Runs")
            ax.set_ylabel("Frequency")
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.tick_params(right=False, top=False)
            plt.locator_params(axis = 'y', nbins=3)

            return fig

#   ###################################################################################################
#    def closeRplot(self):
#        import rpy2.robjects as robjects
#        robjects.r('dev.off()')

    ###################################################################################################
    def setsubset(self, cuts, write=False):
        if self.residuals==None:
            print "subset aborted: you should run diagnostics first"
            return
        else:
            invalid = np.logical_or(self.residuals < cuts[0], self.residuals > cuts[1])

        import pandas as pd


        rootname = self.filename.split('.')[0]
        fname = "%s_subset.csv" % rootname
        
        if write:
            orig = pd.read_csv( self.filename , sep=" ")
            dest = orig[np.logical_not(invalid)]
        
            dest.to_csv(fname,  sep=' ', index=False)
            
            print "subset created, new filename: %s" % fname
            print "you must now run() or load() to update results"
            
        else:
            print "object filename updated"
            print "run() or load() to update results"
        

        self.filename = fname
        self.rootname = self.filename.split('.')[0]


    ###################################################################################################
    def estgamlam(self, writetofile=False):

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

        #print self.returncode
        if writetofile:
            f = open('gamlam.R','w')
            f.writelines(seq)
            f.close()

        try:
            self.readgamlam( gamlamfile )
        except:
            print "Error"
            print out
            print err
            raise RuntimeError("Error in execution within R (see error output above)")


    ##########################################################################
    def readgamlam(self, gamlamfile):

        data = []
        csvfile = open(gamlamfile, 'rb')
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        reader.next()
        for row in reader:
            data.append(np.asarray(row))
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

        self.getRdatafilename(force_refit)

        if not os.path.isfile(self.Rdatafile):
            self.saveRobj = True
            if self.verbose:
                print ".. running analysis.."
            self.run()
        else:
            if self.verbose:
                print "reading results from MLDS file"

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
        
        # GLM object
        self.GLMobject = obsmlds[obsmlds.names.index('obj')]


        # if bootstrapped
        if self.boot:
            obsmns  = robjects.r['obs.mns']
            obslow  = robjects.r['obs.low']
            obshigh = robjects.r['obs.high']
            obssd   = robjects.r['obs.sd']
            bt = robjects.r['samples']

            self.mns = np.array(obsmns)[:-1]
            self.ci95 =  np.vstack( (obslow[:-1], obshigh[:-1]))
            self.sd = np.array(obssd)[:-1]

            self.sigmamns = np.array(obsmns )[-1]
            self.sigmaci95= np.array([obslow[-1], obshigh[-1]])

            arr = np.array(bt)[:-1,:]
            anchor = np.zeros((1, arr.shape[1]))

            self.scalesbt = np.vstack((anchor, arr))
            self.sigmabt = np.array(bt)[-1,:]

            self.status=2

    #########################################################################
    def saveRcommands(self):
        """ Save the commands run by R in each analysis in an R file """
        Rfname = self.Rdatafile.split('.')[0] + '.R'
        f = open(Rfname, 'w')
        f.writelines(self.seq)
        f.close()



###############################################################################
########################## utilities for this class  #########################

def plotscale(s, observer="", color='blue', offset=0, linewidth=1, elinewidth=1, **kargs):

    import matplotlib.pyplot as plt

    if s.boot:
        if s.standardscale:
            label = "%s, $\hat{\sigma}=%.3f \, [%.3f, \, %.3f]$" % (observer, s.sigmamns, s.sigmaci95[0], s.sigmaci95[1])
        else:
            label = "%s" % observer

        #plt.errorbar(s.stim, s.mns, yerr=s.ci95, label=label)
        stimoff= np.hstack((0, s.stim[1:]+ offset))
        if s.ci95.shape[0]==2:
            yerr = abs(s.ci95 - s.mns)
        elif s.ci95.shape[0]==1:
            yerr = s.ci95
        plt.errorbar(stimoff, s.mns, yerr= yerr, fmt="none", ecolor= color, elinewidth=elinewidth)
        plt.plot(s.stim, s.mns, color= color, label=label, linewidth=linewidth, **kargs)

    else:
        if s.standardscale:
            label = "%s, $\hat{\sigma}=%.3f$" % (observer, s.sigma)
        else:
            label = "%s" % observer

        plt.plot(s.stim, s.scale, color= color, label=label, linewidth=linewidth, **kargs)


###############################################################################
###############################################################################
########################## Simulation utilities   #############################

def simulateobserver(sensoryrep, stim, nblocks=1, decisionrule='diff',
                     noisetype='sensory', sigma=0.0, secondstim='indep',
                     debug=False):

    """ Simulate an observer responding to a triads experiment according
    to a given sensory representation function

    Inputs
    ----------

    sensoryrep: callable
                function that returns the value of the sensory representation
                psi at stimulus level x.
                Example: sensoryrep = lambda x: x**2

    stim : vector
            vector of stimulus levels

    nblocks: int. Default : 1
            number of blocks to be simulated

    noisetype: string. Default: 'sensory'
                'sensory' or 'decision', where the gaussian noise originates.
                in case of decision, 'sigma' argument is used to establish
                the decision noise amount

    decisionrule:  string. Default: 'diff'
            'diff', 'absdiff'

    sigma : float.
            Decision noise. Default = 0.0

    secondstim: string. Either 'indep' or 'same'. Default: 'indep'
            tells how to draw the 2nd stimulus, which is the middle stimulus
            in the method of triads, and serves as an anchor. The sensory
            response can be drawn only once for this stimulus ('same') or
            twice ('indep')


    Output
    ----------

    f: string
        filename csv containing data of simulated experiment

    """

    # Simulation triads experiment
    triads, indxt, order = generate_triads(stim)

    # where to save the results
    filename = str(uuid.uuid4()) + '.csv'
    rfl = open(filename, 'wb')
    writer = csv.writer(rfl, delimiter=' ')
    writer.writerow(['Trial', 'Response', 's1', 's2', 's3', 'invord', 'i1', 'i2', 'i3'])

    # going through every trial
    for b in range(nblocks):
        for t in range(len(triads)):

            # stimulus  values of current triad
            if order[t] == 1:
                s1 = triads[t][2]
                s2 = triads[t][1]
                s3 = triads[t][0]
            else:
                s1 = triads[t][0]
                s2 = triads[t][1]
                s3 = triads[t][2]


            # sensory representation variables
            if secondstim=='indep':
                try:
                    f1, f2, f2b, f3 =  sensoryrep([s1, s2, s2, s3])
                    if debug: print "vector form, stimulus 2 drawn twice."

                except:
                    if debug: print "non vector form"
                    f1 = sensoryrep(s1)
                    f2 = sensoryrep(s2)
                    f2b = sensoryrep(s2)
                    f3 = sensoryrep(s3)


            elif secondstim=='same':
                try:
                    f1, f2, f3 = sensoryrep([s1, s2, s3])
                    f2b = f2
                    if debug: print "vector form, stimulus 2 drawn once."
                except:
                    if debug: print "non vector form"
                    f1= sensoryrep(s1)
                    f2 = sensoryrep(s2)
                    f2b = f2
                    f3 = sensoryrep(s3)
            else:
                 raise ValueError('wrong second stimulus rule "%s"' % secondstim)


            ##  check noise sources
            if (sensoryrep.sigmamin!=0 or sensoryrep.sigmamax!=0) and noisetype=='decision':
                print "Note: you have set decision noise AND a sensory function with noise"


            # decision variable
            if decisionrule == 'diff':
                delta = (f3 - f2b) - (f2-f1)
            elif decisionrule == 'absdiff':
                delta = abs(f3 - f2b) - abs(f2 - f1)
            else:
                raise ValueError('wrong decision rule "%s"' % decisionrule)


            if noisetype=='decision':
                delta = delta + random.gauss(0, sigma)
            elif noisetype=='sensory':
                pass
            else:
                raise ValueError('wrong noise type "%s"' % noisetype)

            # binary response
            if delta > 0:
                response = 1
            elif delta < 0:
                response = 0
            else:
                response = np.nan

            if order[t] == 1:
                response = int(not response)

            # saves response
            writer.writerow([t, response, "%.5f" % triads[t][0], "%.5f" % triads[t][1], "%.5f" % triads[t][2], order[t], indxt[t][0]+1, indxt[t][1]+1, indxt[t][2]+1])

    rfl.close()

    if debug and decisionrule == 'diff':
        print "diff"
    elif debug and decisionrule == 'absdiff':
        print "absdiff"

    return filename


###############################################################################
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
#EOF
