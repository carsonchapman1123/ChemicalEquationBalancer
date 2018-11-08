#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg as LA
import sympy as sp
import re
from sympy.parsing.sympy_parser import parse_expr

# takes a string of a chemical equation and returns a string of the balanced equation.
# the equation must have no leading coefficients for it to work properly.
def balance(eq):
    # n is number of molecules
    n = len(re.findall(r'\+|\=', eq)) + 1
    
    # sides is a list containing strings of the two sides of the equations
    sides = re.findall(r'[^\=]+', eq)

    # separate the left and right hand sides into a list of compound names
    lcompounds = re.findall(r'[^\s\+]+', sides[0])
    rcompounds = re.findall(r'[^\s\+]+', sides[1])
    
    # separate the compound names into a list of elements and their counts
    elements = [] # a list of all the unique elements in the list
    lelements = [] # a list of the elements and their counts on the lhs of the equation
    for i in range(len(lcompounds)):
        compound = [] # a list for this compound to be added to the list of separated compounds
        separated = re.findall(r'[A-Z][a-z]*\d*', lcompounds[i])
        for j in range(len(separated)):
            elementName = re.search(r'[A-Z][a-z]*', separated[j]).group()
            # add the element name to the set of elements
            elements.append(elementName)
            elementCount = re.search(r'\d+', separated[j])
            # if the element has a subscript, add it to the list. otherwise set it to 1
            if elementCount:
                compound.append([elementName, int(elementCount.group())])
            else:
                compound.append([elementName, 1])
        # add the separated compound to the list
        lelements.append(compound)
    # do the same type of separation for the right hand side
    relements = []
    for i in range(len(rcompounds)):
        compound = []
        separated = re.findall(r'[A-Z][a-z]*\d*', rcompounds[i])
        for j in range(len(separated)):
            elementName = re.search(r'[A-Z][a-z]*', separated[j]).group()
            elements.append(elementName)
            elementCount = re.search(r'\d+', separated[j])
            if elementCount:
                compound.append([elementName, int(elementCount.group())])
            else:
                compound.append([elementName, 1])
        relements.append(compound)
    
    # remove duplicates from the list of elements
    elements = list(set(elements))
    
    # llength is the length of the left hand side, which will be used to add
    # numbers from the right hand side of the equation into the proper position
    # in the matrix
    llength = len(lelements)
    
    # now create the empty matrix
    matrix = [[0 for _ in range(n+1)] for _ in range(len(elements))]
    # go element by element adding it to the matrix
    # i is index of current element
    for i in range(len(elements)):
        # find the current element name to search for
        elementName = elements[i]
        # check for that element on left hand side and add its count to matrix
        # j is current molecule number
        for j in range(len(lelements)):
            for k in range(len(lelements[j])):
                # if the element we are searching for is in this molecule
                if lelements[j][k][0] == elementName:
                    # then add it to the corresponding row in the matrix
                    matrix[i][j] += lelements[j][k][1]
        # now check for element on right hand side and subtract it from matrix
        for j in range(len(relements)):
            for k in range(len(relements[j])):
                if relements[j][k][0] == elementName:
                    matrix[i][j+llength] -= relements[j][k][1]
    
    # now solve matrix for unknowns
    M = sp.Matrix(matrix)
    x=[parse_expr('x%d'%i) for i in range(n)]
    x=sp.symbols('x0:%d'%n)
    sols=sp.solve_linear_system(M,*x)
    
    # substitute in 1 for any free variable    
    for item in sols:
        for j in range(len(x)):
            sols[item] = sols[item].subs(x[j],1)
    
    # create a list of solutions and replace any free variable with 1
    solList = [0 for _ in range(n)]
    for i in range(n):
        varname=sp.symbols("x"+str(i))
        # if the variable is in the solution matrix, add it to the solution list
        if varname in sols.keys():
            solList[i] = sols[varname]
        else:
            # otherwise add in 1 as its value
            solList[i] = sp.Rational(1,1)
    # create a list of denominators to find the lcm
    denominators = []
    for i in range(len(solList)):
        denominators.append(solList[i].q)
    lcm = sp.lcm(denominators)
    # multiply each solution by the lcm
    for i in range(len(solList)):
        solList[i] *= lcm

    # create string to return
    returnString = ""
    for i in range(len(lcompounds)):
        if solList[i] == 1:
            returnString = returnString + lcompounds[i] + "+"
        else:
            returnString = returnString + str(solList[i]) + lcompounds[i] + "+"
    # remove the one extra "+" at the end of the return string
    returnString = returnString[:-1]
    # add the "="
    returnString = returnString + "="
    for i in range(len(rcompounds)):
        if solList[i + llength] == 1:
            returnString = returnString + rcompounds[i] + "+"
        else:
            returnString = returnString + str(solList[i+llength]) + rcompounds[i] + "+"
    # remove the one extra "+" at the end of the return string
    returnString = returnString[:-1]
    return returnString
