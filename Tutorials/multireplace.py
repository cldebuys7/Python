#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 01:51:38 2019

@author: cldebuys
This function is from bgusach (not from the author).
source => https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
"""

import re

def multireplace(string, replacements):
    """
    Given a string and a replacement map, it returns the replaced string.
    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary {value to find: value to replace}
    :rtype: str
    """
    # Place longer ones first to keep shorter substrings from matching where the longer ones should take place
    # For instance given the replacements {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
    # 'hey ABC' and not 'hey ABc'
    substrs = sorted(replacements, key=len, reverse=True)

    # Create a big OR regex that matches any of the substrings to replace
    regexp = re.compile('|'.join(map(re.escape, substrs)))

    # For each match, look up the new string in the replacements
    return regexp.sub(lambda match: replacements[match.group(0)], string)

################################## UNIT TEST ##################################
test = "random string q0 q1 dq1 dq0 dq2"
rep = {"q0": "q[0]","q1":"q[1]","dq0":"q[5]","dq1": "q[6]"}
testrep = multireplace(test,rep)

test2 = "Matrix([[-R*dq1**2*l*m1*(I1 + l**2*m1)*sin(q1)/(-R**2*l**2*m1**2*cos(q1)**2 + (I0 + R**2*(m0 + m1))*(I1 + l**2*m1)) - R*g*l**2*m1**2*sin(q1)*cos(q1)/(-R**2*l**2*m1**2*cos(q1)**2 + (I0 + R**2*(m0 + m1))*(I1 + l**2*m1))],[-R**2*dq1**2*l**2*m1**2*sin(q1)*cos(q1)/(-R**2*l**2*m1**2*cos(q1)**2 + (I0 + R**2*(m0 + m1))*(I1 + l**2*m1)) - g*l*m1*(I0 + R**2*(m0 + m1))*sin(q1)/(-R**2*l**2*m1**2*cos(q1)**2 + (I0 + R**2*(m0 + m1))*(I1 + l**2*m1))]])"
rep2 = {"q0": "x[0]","q1":"x[1]","dq0":"x[2]","dq1": "x[3]"}
testrep2 = multireplace(test2,rep2)
###############################################################################