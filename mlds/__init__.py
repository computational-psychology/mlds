# -*- coding: utf-8 -*-
"""
Python MLDS - based on R MLDS by Ken Knoblauch.
"""

from .mlds import MLDSObject, simulateobserver, generate_quadruples, generate_triads, plotscale

def test():
    import pytest, os
    cwd = os.path.abspath(os.path.dirname(__file__))
    pytest.main(['-x', cwd])