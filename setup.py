#
# -*- coding: utf-8 -*-
#

import os
from setuptools import setup, find_packages

# Get the long description from the README file
cwd = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(cwd, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Programming Language :: Python :: 3',
    ],
    packages=find_packages(),
    package_data={'mlds.test': ['*.csv', '*.MLDS']},
    install_requires=[
        'scipy',
        'numpy',
        'matplotlib',
        'rpy2>=2.3.10',
        'pytest',
        'joblib',
        'pandas',
        'seaborn',
    ],

    name='mlds',
    description='Python wrapper for MLDS R package',
    long_description=long_description,
    author='Guillermo Aguilar',
    license='GPL2',
    version='0.2',
    url='https://github.com/computational-psychology/mlds',
)
