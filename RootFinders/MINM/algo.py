import numpy as np
import sympy as symp
from sympy import *
from helper import *

x = symp.Symbol('x')
f = ((x-1)**5)*(x-3)
f_prime = symp.diff(f)

### Step 1 ###
xk = [0,2]
epsilon = 0.001
k = 0
dif = 1

### Step 2 ###
while abs(dif) > epsilon/100:
    f_prime_xk = function_image(f_prime,xk[0],xk[1])

### Step 3 ###
    xk_mid = (xk[0] + xk[1])/2
    N_tilde = interval_subtraction([xk_mid,xk_mid],interval_division([f.subs(x,xk_mid),f.subs(x,xk_mid)],f_prime_xk))
    xk_tilde = interval_intersection(xk,N_tilde)

### Step 4 ###
    f_prime_xk_tilde = function_image(f_prime,xk_tilde[0],xk_tilde[1])

### Step 5 ###
    xk_tilde_mid = (xk_tilde[0] + xk_tilde[1])/2
    f_xk_tilde = function_image(f,xk_tilde[0],xk_tilde[1])
    f_prime_xk_tilde = function_image(f_prime,xk_tilde[0],xk_tilde[1])
    N_tilde_hat = interval_subtraction([xk_tilde_mid,xk_tilde_mid],interval_division([2*f.subs(x,xk_tilde_mid),2*f.subs(x,xk_tilde_mid)], \
                                                                                    interval_addition(f_prime_xk,f_prime_xk_tilde)))
    xk = interval_intersection(xk,N_tilde_hat)
    if xk[1] < xk[0]:
        xk = [xk[1], xk[0]]
    dif = xk[1] - xk[0]
    print(xk)