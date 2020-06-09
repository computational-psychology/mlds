# -*- coding: utf-8 -*-
"""
Unittests for sensoryfunctions module

@author: G. Aguilar, Sep 2015
"""

from mlds.sensoryfunctions import PowerSensoryFunc, QuadraticSensoryFunc
import unittest
import numpy as np


d = 2  # decimal positions to compare

class TestPowerfunctionClass(unittest.TestCase):
    
    def setUp(self):
        global fn
        fn = PowerSensoryFunc()
        
        
        fn.sigmamin = 0.0
        fn.sigmamax = 0.0
        fn.exponent = 3
        fn.a = 1
        fn.c = 0
    
        
    def test_call_single(self):
        x = 0.5
        self.assertAlmostEqual(fn(x), x**3, places=d)  
  

    def test_call_vector(self):
        x = np.linspace(0, 1, 10)
        np.testing.assert_almost_equal(fn(x), x**3)  
        
    def test_call_fail(self):
        self.assertRaises(ValueError, fn, 'a')
      


class TestQuadraticfunctionClass(unittest.TestCase):
    
    def setUp(self):
        global fn
        fn = QuadraticSensoryFunc()
        
        fn.sigmamin = 0.0
        fn.sigmamax = 0.0
        fn.a = 2
        fn.b = 2
        fn.c = 0
    
        
    def test_call_single(self):
        x = 0.5
        self.assertAlmostEqual(fn(x), 2*x**2+2*x, places=d)  
  

    def test_call_vector(self):
        x = np.linspace(0, 1, 10)
        np.testing.assert_almost_equal(fn(x), 2*x**2+2*x)  
        
    def test_call_fail(self):
        self.assertRaises(ValueError, fn, 'a')
        
        
        
    
if __name__ == '__main__':
    unittest.main()
    
