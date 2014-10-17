# -*- coding: utf-8 -*-
"""
Unittests for mlds module

@author: G. Aguilar, Oct 2014
"""

import sys
sys.path.append('../')
import mlds
import mlds_gam
import unittest
import numpy as np

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
prob = 0.927

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
        obs = mlds.MLDSObject('test.csv', boot=False, save=False)       
        obs.run()
        
        np.testing.assert_almost_equal(obs.scale, scale, decimal= d)
        
    
    @unittest.skip("skipping bootstrap, saving time")
    def test_bootstrap(self):
        obs = mlds.MLDSObject('test.csv', boot=True, save=False)
        obs.parallel=True
        obs.master= '"localhost"'
        obs.workers= ['"localhost"']*4
        
        obs.run()
        
        self.compare(obs)
    
    @unittest.skip("skipping bootstrap, saving time")
    def test_bootstrap_correctedCI(self):
        obs = mlds.MLDSObject('test.csv', boot=True, save=False)
        obs.parallel=True
        obs.master= '"localhost"'
        obs.workers= ['"localhost"']*4
        obs.correctedCI=True
        obs.run()
        
        np.testing.assert_almost_equal(obs.ci95, ci95_corrected, decimal=d)
        np.testing.assert_almost_equal(obs.sigmaci95, sigmaci95_corrected, decimal=d)
        
        
    def test_readcsv(self):
        obs = mlds.MLDSObject('test.csv', boot=True, keepfiles=True)
        obs.mldsfile='output_mldsfile.csv'
        obs.returncode=0

        obs.readresults()                                    
        
        self.compare(obs)
        
    def test_readobjectresults(self):
        obs = mlds.MLDSObject('test.csv', boot=True, keepfiles=False)
        obs.Rdatafile = 'output_test.MLDS'
        
        obs.readobjectresults()
        
        self.compare(obs)


    def test_readdiags(self):
        obs = mlds.MLDSObject('test.csv', boot=False, keepfiles=False)
        obs.Rdatafile = 'output_test.MLDS'
        
        obs.readdiags()
        
        self.assertAlmostEqual(obs.AIC, aic, places=d)
        self.assertAlmostEqual(obs.DAF, daf, places=d)

        
#######################################################################
class TestMLDSComparison(unittest.TestCase):
        

    def test_gam_run4(self):
        
        files = ['0.csv', '1.csv', '2.csv', '3.csv']
        gam = mlds_gam.MLDSGAMCompare(files)
        gam.run()
        
    def test_gam_run3(self):
        
        files = ['1.csv', '2.csv', '3.csv']
        gam = mlds_gam.MLDSGAMCompare(files)
        gam.run()
        
    def test_gam_run2(self):
        
        files = ['1.csv', '2.csv']
        gam = mlds_gam.MLDSGAMCompare(files)
        gam.run()
        
    def test_gam_run1(self):
        
        files = ['0.csv']
        gam = mlds_gam.MLDSGAMCompare(files, dividedby=2)
        gam.run()
        
    def test_gam_dividedby(self):
        
        with self.assertRaises(Exception):
            mlds_gam.MLDSGAMCompare(['0.csv'])
            
    def test_gam_argument(self):
        
        with self.assertRaises(Exception):
            mlds_gam.MLDSGAMCompare('0.csv')


        

    
if __name__ == '__main__':
    unittest.main()
    
