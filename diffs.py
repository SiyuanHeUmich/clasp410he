#!/usr/bin/env python3
'''
Let's exploring numerical differenceing!

'''
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')


def fwd_diff(y,dx):
    '''
    Return forward diff approx of 1st derivative
    '''
    dydx = np.zeros(y.size)

    # Forward diff:
    dydx[:-1] = (y[1:] - y[:-1]) / dx

    # last point:
    dydx[-1] = (y[-1] - y[-2]) / dx
    return dydx

def bwd_diff(y,dx):
    '''
    Return backward diff approx of 1st derivative
    '''
    dydx = np.zeros(y.size)

    # backward diff:
    dydx[1:] = (y[1:] - y[:-1]) / dx

    # first point:
    dydx[0] = (y[1] - y[0]) / dx
    return dydx

deltax=0.1
x= np.arange(0, 4*np.pi, deltax)

fx = np.sin(x)
fxd1=np.cos(x)

fig,ax= plt.subplots(1,1)

#Plotting actual sin(x) and actual derivative (cos(x))
ax.plot(x, fx, '.', alpha=.6, label=r'$f(x) = \sin(x)$')
ax.plot(x, fxd1, label=r'$f(x) = \frac{d\sin(x)}{dx}$')

# Plotting two different metods of difference
ax.plot(x, fwd_diff(fx, deltax), label='Fwd Diff')
ax.plot(x, bwd_diff(fx, deltax), label='Bwd Diff')

ax.set_title('Forward and Backward Difference Approximations of $\sin(x)$ Derivative')
ax.legend(loc='upper right')