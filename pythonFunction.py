#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 02:35:57 2019

@author: cldebuys
"""
import re

# input: string = string to be altered, 
#        replacemements = dict of replacements {'old0':'new0','old1':'new1'}
# output: new string with alterations

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



# input: funcname = name of new function [string]
#        filename = name of new function file [string]
#           - if left empty [], funcname.py is used
#        eqn = equation to evaluate in new function [string]
#        varsToReplace = symbols to replace [list of strings]
#        importString = any necessary imports [list of strings]
#        ode_opt = 1 if func is for ode (t,x), otherwise 0 (x)
# output: creates a new file containing the function
    
def pythonFunction(funcname,filename,eqn,varsToReplace,importString='',ode_opt=0):
    # use varsToReplace to make dictionary with the replacement map
    # example: # rep = {"q0": "x[0]","q1":"x[1]","dq0":"x[2]","dq1": "x[3]"}
    rep = {}
    for i in range(len(varsToReplace)):
        rep[varsToReplace[i]] = 'y[' + str(i) + ']'
        
    # replace old symbols q_i with new x[i]
    new_eqn = multireplace(eqn,rep) # this function is defined above
    
    # if no filename is given, use funcname
    if not filename:
        filename = funcname + ".py"
        
    # define function arguments (add 't' if this func will be called by an ode)
    arguments = "(y):"
    if ode_opt == 1:
        arguments = "(t,y,par):" # array of parameters to pass

    # create new .py file, where new_eqn is the function to be evaluated
    with open(str(filename), 'w') as filehandle:
        filebuffer = importString + ["",
                     "def " + funcname + arguments,
                     "",
                     "\treturn " + new_eqn]
        filehandle.writelines("%s\n" % line for line in filebuffer)
    return()



# This is an example of the file writing portion alone
# =============================================================================
# with open('M_matrix.py', 'w') as filehandle:  
#     filebuffer = ["from sympy import Matrix, cos, sin",
#                   "",
#                   "def M_matrix(t,x):",
#                   "",
#                   "\treturn " + "new_eqn"]
#     filehandle.writelines("%s\n" % line for line in filebuffer)
# =============================================================================
    
# add preceding code which replaces state variables for the user