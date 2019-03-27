#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 01:02:20 2019

@author: cldebuys
"""

from sympy import Symbol
a=Symbol('a')
c=Symbol('c')

y=a+3*c

import pickle
with open("y.py", "w") as outf:
    pickle.dump(y, outf)
    
with open("y.py") as inf:
    x = pickle.load(inf)