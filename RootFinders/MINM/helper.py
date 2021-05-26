import numpy as np
import sympy as symp
from sympy import *

x = symp.Symbol('x')

def interval_division(x,y):
    '''x/y for intervals x and y
    :param x: interval dividend or numerator
    :type x: list
    :param y: interval divisor or denominator
    :type y: list
    '''
    z = []
    boundaries = []
    for num in x:
        for den in y:
            if den != 0 :
                boundaries.append(num/den)
    z.append(min(boundaries))
    z.append(max(boundaries))
    return z

def interval_subtraction(x,y):
    '''x - y for intervals x and y
    :param x: interval minuend
    :type x: list
    :param y: interval subtrahend
    :type y: list
    '''
    z = []
    z.append(x[0]-y[1])
    z.append(x[1]-y[0])
    return z

def interval_addition(x,y):
    '''x + y for intervals x and y
    :param x: interval addend
    :type x: list
    :param y: interval addend
    :type y: list
    '''
    z = []
    z.append(x[0] + y[0])
    z.append(x[1] + y[1])
    return z

def interval_intersection(x,y):
    ''' intersection of interval x and interval y
    :param x: interval 
    :type x: list
    :param y: interval
    :type y: list
    ''' 
    z = []
    z.append(max([x[0],y[0]]))
    z.append(min([x[1],y[1]]))
    return z

def get_critical_points(f,a,b):
    '''Returns list of critical points of f over the interval [a,b]
    :param f: function (over symbol 'x')
    :type f: symp function
    :param a: left bound of interval
    :type a: float
    :param b: right bound of interval
    :type b: float
    '''
    f_prime = symp.diff(f)
    critical_points = solve(f_prime,x)
    for cp in critical_points:
        if cp < a or cp > b:
            critical_points.remove(cp)
    critical_points.extend([a,b])
    return critical_points

def function_image(f,a,b):
    '''Returns image of f under the interval [a,b]
    :param f: function (over symbol 'x')
    :type f: symp function
    :param a: left bound of interval
    :type a: float
    :param b: right bound of interval
    :type b: float
    '''
    f_cps = get_critical_points(f,a,b)
    cp_vals = []
    for cp in f_cps:
        cp_vals.append(f.subs(x,cp))
    f_image = [min(cp_vals),max(cp_vals)]
    return f_image