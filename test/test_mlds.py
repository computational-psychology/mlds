# -*- coding: utf-8 -*-
"""
Unittests for mlds module

@author: G. Aguilar, Oct 2014
"""

import os
import mlds
import unittest
import numpy as np

base_path = os.path.dirname(__file__)
def abspath(fname):
    return os.path.join(base_path, fname)


# expected results from test.csv
scale = np.array([ 0. , 0.0203168884, 0.0657919702, 0.1426437, 0.1758397487, 0.2587451112, 0.3697451673, 0.4686494742, 0.6344124339, 0.7298032188, 1. ])
mns = np.array([ 0. ,  0.0193694154, 0.0650221787, 0.1420025599, 0.1753727203, 0.2583533223, 0.3692932295, 0.4683024835, 0.6341455306, 0.7301011469, 1.        ])
low = np.array([ 0. , -0.0272724567, 0.0210245681, 0.1014956201, 0.1354316387, 0.2210277348, 0.3336706154, 0.4339024888, 0.5986615941, 0.6890368767,  1.        ])  
high = np.array([ 0., 0.0615440407, 0.1053452521, 0.1796633018, 0.2128018444, 0.2944877629, 0.4028255524, 0.5018942114, 0.6719317338, 0.7750129988,  1.        ])
ci95 = np.vstack((low, high))

sigma = 0.1596638086
sigmamns = 0.1568117077
sigmaci95 = np.array([0.133503277554566, 0.181765764894874])

aic = 570.3924
daf = 0.5187585
prob = 0.93

# corrected CI
low_C= np.array([0.00000000, -0.02149037, 0.02468978, 0.10483483, 0.13793712, 0.22326763, 0.33504963, 0.43554262, 0.59778788, 0.68618777, 1.00000000])
high_C=np.array([0.00000000, 0.06550406, 0.10849741, 0.18226010, 0.21574120, 0.29597577, 0.40533268, 0.50275388, 0.67041595, 0.77105617, 1.00000000])
ci95_corrected = np.vstack((low_C, high_C))

sigmaci95_corrected= np.array([0.13340837, 0.18122606])

#####
d = 2  # decimal positions to compare

########################################################################

class TestMLDSClass(unittest.TestCase):
    
    def compare(self, obs):
        
        np.testing.assert_almost_equal(obs.sigma, sigma,decimal= d)
        np.testing.assert_almost_equal(obs.scale, scale, decimal= d)
        
        np.testing.assert_almost_equal(obs.mns, mns, decimal=d)
        np.testing.assert_almost_equal(obs.ci95, ci95, decimal=d)
        np.testing.assert_almost_equal(obs.sigmamns, sigmamns, decimal=d)
        np.testing.assert_almost_equal(obs.sigmaci95, sigmaci95, decimal=d)
        
        
    def test_simple(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False)
        obs.run()
        
        np.testing.assert_almost_equal(obs.scale, scale, decimal= d)
        
    
    @unittest.skip("skipping bootstrap, saving time")
    def test_bootstrap(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, save=False)
        obs.parallel=False
        obs.nsamples = 1000
        obs.initcommands()
        obs.run()
        
        self.compare(obs)
    
    @unittest.skip("skipping bootstrap, saving time")
    def test_bootstrap_correctedCI(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, save=False)
        obs.correctedCI=True
        obs.parallel=False
        obs.nsamples = 1000
        obs.run()
        
        np.testing.assert_almost_equal(obs.ci95, ci95_corrected, decimal=d)
        np.testing.assert_almost_equal(obs.sigmaci95, sigmaci95_corrected, decimal=d)
        
        
    def test_Rdatafilename(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False)
        print(obs.rootname, obs.filename)
        obs.getRdatafilename()
        assert(obs.Rdatafile == abspath('test_norm_probit.MLDS'))



    def test_Rdatafilename2(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False, standardscale=False)

        obs.getRdatafilename()
        assert(obs.Rdatafile == abspath('test_unnorm_probit.MLDS'))


    def test_Rdatafilename3(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False, standardscale=False)
        obs.linktype = 'cauchit'
        
        obs.getRdatafilename()
        assert(obs.Rdatafile == abspath('test_unnorm_cauchit.MLDS'))


    def test_Rdatafilename4(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False, standardscale=False)
        
        obs.getRdatafilename(force_refit=True)
        assert(obs.Rdatafile == abspath('test_unnorm_probit_refit.MLDS'))
        
        

    def test_readcsv(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, keepfiles=True)
        obs.mldsfile = abspath('output_mldsfile.csv')
        obs.returncode=0

        obs.readresults()                                    
        
        self.compare(obs)
        
    def test_readobjectresults(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, keepfiles=False)
        obs.Rdatafile = abspath('output_test.MLDS')
        
        obs.readobjectresults()
        
        self.compare(obs)

    @unittest.skip("skipping bootstrap diagnostics, saving time")
    def test_rundiags(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, keepfiles=False)
        obs.parallel=False
        obs.nsamples = 250
        obs.run()
        obs.rundiagnostics()
                
        self.assertAlmostEqual(obs.prob, prob, places=1)
        os.remove(obs.Rdatafile)
    
    @unittest.skip("skipping bootstrap diagnostics, saving time")
    def test_rundiags_nosave(self, saveresiduals=False):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, keepfiles=False)
        obs.parallel=False
        obs.nsamples = 250
        obs.run()
        obs.rundiagnostics(saveresiduals=saveresiduals)
                
        self.assertAlmostEqual(obs.prob, prob, places=1)    
        os.remove(obs.Rdatafile)
       
       
    def test_readdiags(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, keepfiles=False)
        obs.Rdatafile = abspath('output_test.MLDS')
        
        obs.readdiags()
        
        self.assertAlmostEqual(obs.AIC, aic, places=d)
        self.assertAlmostEqual(obs.DAF, daf, places=d)
        
    def test_simple_othernamecols(self):
        obs = mlds.MLDSObject(abspath('test_cols.csv'), boot=False, save=False)
        obs.colnames = {'stim' : ['stim1', 'stim2', 'stim3'], 
                        'response': 'resp'}
        obs.run()
        
        np.testing.assert_almost_equal(obs.scale, scale, decimal= d)


    def test_simple_othernamecols2(self):
        obs = mlds.MLDSObject(abspath('test_cols.csv'), boot=False, save=False)
        obs.colnames = {'stim' : ['stim1', 'stim2', 'stim3'], 
                        'response': 'resp'}
        obs.load()
        
        np.testing.assert_almost_equal(obs.scale, scale, decimal= d)
        
        
    def test_dimension_unit(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False,
                              dimension_unit='stim')
        #obs.saveRcommands()
        obs.run()
        
        
    def test_dimension_unit_Rdatafilename(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=False, save=False,
                              dimension_unit='stim')
                
        obs.getRdatafilename()
        assert(obs.Rdatafile == abspath('test_stim_norm_probit.MLDS'))

    @unittest.skip("skipping threshold prediction")
    def test_threshold_prediction(self):
        obs = mlds.MLDSObject(abspath('test.csv'), boot=True, standardscale=False)
        obs.nsamples = 1000
        obs.parallel = False
        obs.load()

        sts = np.array([0.25, 0.99])
        dp = [-2, -0.5, 0.5, 2]

        k = 3
        factor = 2.0
        res = 0.000005
        citype = 'BCa'
        mldsthrs, *_ = mlds.predict_thresholds(obs, sts, dp, k=k, factor=factor,
                                                     citype=citype, tol=0.01, res=res,
                                                     warn=False, debug=False, save=False)

        # comparing point estimate
        expected = np.array([np.nan, 0.17, 0.32, 0.49, 0.91, 0.97, np.nan, np.nan])
        actual = mldsthrs['point_estimate'].to_list()
        np.testing.assert_almost_equal(expected, actual, decimal=d)
        
        # comparing lower bound CI
        expected = np.array([np.nan, 0.14, 0.31, 0.46, 0.90, 0.97, np.nan, np.nan])
        actual = mldsthrs['CIl'].to_list()
        np.testing.assert_almost_equal(expected, actual, decimal=d)
        
        # comparing upper bound CI
        expected = np.array([np.nan, 0.18, 0.335, 0.52, 0.926, 0.9747, np.nan, np.nan])
        actual = mldsthrs['CIh'].to_list()
        np.testing.assert_almost_equal(expected, actual, decimal=d)
        
        
        
# #######################################################################
# class TestMLDSComparison(unittest.TestCase):
#
#
#     def test_gam_run4(self):
#
#         files = ['0.csv', '1.csv', '2.csv', '3.csv']
#         gam = mlds.MLDSGAMCompare(files)
#         gam.run()
#
#     def test_gam_run3(self):
#
#         files = ['1.csv', '2.csv', '3.csv']
#         gam = mlds.MLDSGAMCompare(files)
#         gam.run()
#
#     def test_gam_run2(self):
#
#         files = ['1.csv', '2.csv']
#         gam = mlds.MLDSGAMCompare(files)
#         gam.run()
#
#     def test_gam_run1(self):
#
#         files = ['0.csv']
#         gam = mlds.MLDSGAMCompare(files, dividedby=2)
#         gam.run()
#
#     def test_gam_dividedby(self):
#
#         with self.assertRaises(Exception):
#             mlds.MLDSGAMCompare(['0.csv'])
#
#     def test_gam_argument(self):
#
#         with self.assertRaises(Exception):
#             mlds.MLDSGAMCompare('0.csv')


        

    
if __name__ == '__main__':
    unittest.main()
    
