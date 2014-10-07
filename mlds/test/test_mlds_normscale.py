# -*- coding: utf-8 -*-
"""
Unittests for mlds module

@author: G. Aguilar, Oct 2014
"""

import sys
sys.path.append('../')
import mlds
import unittest
import numpy as np

# data
scale = np.array([ 0.        ,  0.02031689,  0.06579197,  0.1426437 ,  0.17583975,
        0.25874511,  0.36974517,  0.46864947,  0.63441243,  0.72980322,  1.        ])
mns = np.array([ 0.        ,  0.01937689,  0.0649507 ,  0.14204127,  0.17516151,
        0.2584842 ,  0.36942893,  0.46859708,  0.63442393,  0.73018295,  1.        ])
low = np.array([ 0.        , -0.02625276,  0.02090501,  0.10185435,  0.13461994,
        0.22115107,  0.33420499,  0.43402753,  0.59919975,  0.68943627,  1.        ])  
high = np.array([ 0.        ,  0.06218987,  0.10517804,  0.18024196,  0.21305171,
        0.29399166,  0.40292939,  0.50239248,  0.67209074,  0.77305887,  1.        ])
ci95 = np.vstack((low, high))

sigma = 0.159663808552398
sigmamns = 0.156704779330266
sigmaci95 = np.array([0.133561621842693, 0.181506665443706])

d = 2

aic = 570.3924
daf = 0.5187585
prob = 0.927

#

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
        
    
    @unittest.skip("skipping bootstrap")
    def test_bootstrap(self):
        obs = mlds.MLDSObject('test.csv', boot=True, save=False)
        obs.parallel=True
        obs.master= '"localhost"'
        obs.workers= ['"localhost"']*4
        
        obs.run()
        
        self.compare(obs)

        
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

        


if __name__ == '__main__':
    unittest.main()
    
