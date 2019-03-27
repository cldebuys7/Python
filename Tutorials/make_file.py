#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 19:03:46 2019

@author: cldebuys
"""
fname = 'generated.py'
data = 6

# =============================================================================
with open(fname, 'w') as f:
    f.write('data = [{}]'.format(data))
# =============================================================================
# =============================================================================
# with open(fname, 'w') as f:
#     f.writelines(['data = [{}]'.format(data), "line2", "line3"])
# =============================================================================

import generated
print(generated.data)
