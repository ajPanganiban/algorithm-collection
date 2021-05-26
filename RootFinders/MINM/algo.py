import numpy as np
import sympy as symp
from sympy import *
from helper import *

x = symp.Symbol('x')

def minm(f,interval,epsilon):
    '''Implements the Modified Interval Newton Method by Ansary and Panda.
        This algorithm seeks to find a root of f inside a given interval.
    :param f: function (over the symbol 'x')
    :type f: symp function
    :param interval: close interval cosidered
    :type interval: list (of endpoints)
    :param epsilon: tolerance value
    :type epsilon: float
    :returns: interval of length smaller than epsilon containing a root of f
    '''
    f_prime = symp.diff(f)

    xk = interval
    dif = 1

    while abs(dif) > epsilon:
        f_prime_xk = function_image(f_prime,xk[0],xk[1])
        xk_mid = (xk[0] + xk[1])/2
        N_tilde = interval_subtraction([xk_mid,xk_mid],interval_division([f.subs(x,xk_mid),f.subs(x,xk_mid)],f_prime_xk))
        xk_tilde = interval_intersection(xk,N_tilde)
        f_prime_xk_tilde = function_image(f_prime,xk_tilde[0],xk_tilde[1])
        xk_tilde_mid = (xk_tilde[0] + xk_tilde[1])/2
        f_prime_xk_tilde = function_image(f_prime,xk_tilde[0],xk_tilde[1])
        N_tilde_hat = interval_subtraction([xk_tilde_mid,xk_tilde_mid],interval_division(
                                                                        [2*f.subs(x,xk_tilde_mid),2*f.subs(x,xk_tilde_mid)], \
                                                                        interval_addition(f_prime_xk,f_prime_xk_tilde)))
        xk = interval_intersection(xk,N_tilde_hat)
        if xk[1] < xk[0]:
            xk = [xk[1], xk[0]]
        dif = xk[1] - xk[0]
    return xk