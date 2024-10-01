#!/usr/bin/env python3

'''
Draft for lab 02, competition equation only
'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

#define Lotka-Volterra eqns
def dNdt_comp(t, N, a=1, b=2, c=1, d=3): 
    '''
    This function calculates the Lotka-Volterra competition equations for
    two species. Given normalized populations, `N1` and `N2`, as well as the
    four coefficients representing population growth and decline,
    calculate the time derivatives dN_1/dt and dN_2/dt and return to the
    caller.
    This function accepts `t`, or time, as an input parameter to be
    compliant with Scipy's ODE solver. However, it is not used in this
    function.
    Parameters
    ----------
    t : float
        The current time (not used here).
    N : two-element list
        The current value of N1 and N2 as a list (e.g., [N1, N2]).
    a, b, c, d : float, defaults=1, 2, 1, 3
        The value of the Lotka-Volterra coefficients.
    Returns
    -------
    dN1dt, dN2dt : floats
        The time derivatives of `N1` and `N2`.
    '''
    # Here, N is a two-element list such that N1=N[0] and N2=N[1]
    dN1dt = a*N[0]*(1-N[0]) - b*N[0]*N[1]
    dN2dt = c*N[1]*(1-N[1]) - d*N[1]*N[0]

    return dN1dt, dN2dt

# define ODE solvers
def euler_solve(func, N1_init=.3, N2_init=.6, dt=1, t_final=100.0):
    '''
    <Your good docstring here>
    Parameters
    ----------
    func : function
        A python function that takes `time`, [`N1`, `N2`] as inputs and
        returns the time derivative of N1 and N2.
    N1_init : float
        <more good docstring here>
    '''
    # Create time array. We won't use that here, but will return it
    # to the caller for convenience.
    t = np.arange(0, t_final, dt)

    # Create container for the solution, set initial condition.
    N1 = np.zeros(t.size)
    N1[0] = N1_init
    N2 = np.zeros(t.size)
    N2[0] = N2_init

    for i in range(1, t.size):
        dN1, dN2 = func(i, [N1[i-1], N2[i-1]] ) # get time derivatives
        # get solutions with Euler formula
        N1[i] = N1[i-1] + dt * dN1
        N2[i] = N2[i-1] + dt * dN2

    return t, N1, N2


def solve_rk8(func, N1_init=.3, N2_init=.6, dT=10, t_final=100.0,
              a=1, b=2, c=1, d=3):
    '''
    Solve the Lotka-Volterra competition and predator/prey equations using
    Scipy's ODE class and the adaptive step 8th order solver.
    Parameters
    ----------
    func : function
        A python function that takes `time`, [`N1`, `N2`] as inputs and
        returns the time derivative of N1 and N2.
    N1_init, N2_init : float
        Initial conditions for `N1` and `N2`, ranging from (0,1]
    dT : float, default=10
        Largest timestep allowed in years.
    t_final : float, default=100
        Integrate until this value is reached, in years.
    a, b, c, d : float, default=1, 2, 1, 3
        Lotka-Volterra coefficient values
    Returns
    -------
    time : Numpy array
        Time elapsed in years.
    N1, N2 : Numpy arrays
        Normalized population density solutions.
    '''
    from scipy.integrate import solve_ivp
    # Configure the initial value problem solver
    result = solve_ivp(func, [0, t_final], [N1_init, N2_init],
                       args=[a, b, c, d], method='DOP853', max_step=dT)
    
    # Perform the integration
    time, N1, N2 = result.t, result.y[0, :], result.y[1, :] # Return values to caller.
    return time, N1, N2

t_comp_eu, N1_comp_eu, N2_comp_eu = euler_solve(dNdt_comp)
t_comp_rk, N1_comp_rk, N2_comp_rk = solve_rk8(dNdt_comp)

# Plot
n1_color = 'red'
n2_color = 'blue'
fig, ax = plt.subplots(1,1)
ax.plot(t_comp_eu,N1_comp_eu, label='N1 - euler', color=n1_color)
ax.plot(t_comp_eu,N2_comp_eu, label='N2 - euler', color=n2_color)
ax.plot(t_comp_rk, N1_comp_rk, ':', label='N1 - rk8', color='yellow')
ax.plot(t_comp_rk, N2_comp_rk, ':',label='N2 - rk8', color='green')
ax.set_ylabel('Normalized population density')
ax.set_xlabel('Time [years]')
ax.set_title('Competition')
fig.legend()

plt.show()