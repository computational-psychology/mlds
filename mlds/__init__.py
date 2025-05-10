# -*- coding: utf-8 -*-
"""
Python MLDS - based on R MLDS by Ken Knoblauch.
"""

__version__ = "0.8.0"

from .mlds import MLDSObject, simulateobserver, generate_quadruples, generate_triads, plotscale
from .threshold_prediction import predict_thresholds

#def test():
#    import pytest, os
#    cwd = os.path.abspath(os.path.dirname(__file__))
#    pytest.main(['-x', cwd])
