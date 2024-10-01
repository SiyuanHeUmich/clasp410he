#!/usr/bin/env python3

'''
This file performs low-order and high-order Ordinary Differential Equations solvers
in both the Lotka-Volterra Competition Scenarios and the Lotka-Volterra 
Predator-Prey Scemarios.

To get solution for lab 2: Run these commands:

>>> run lab02.py
>>> plot.ion()
>>> plot_q1_q2_results()
>>> plot_q3_results()

'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

def dNdt_comp(t, N, a=1, b=2, c=1, d=3):
    '''
    This function calculates the Lotka-Volterra competition equations for
    two species. Given normalized populations, `N1` and `N2`, as well as the
    four coefficients representing population growth and decline,
    calculate the time derivatives dN_1/dt and dN_2/dt and return to the caller.
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


def dNdt_prey(t, N, a=1, b=2, c=1, d=3):
    '''
    This function calculates the Lotka-Volterra predator-prey equations for
    two species. Given normalized populations, `N1` and `N2`, as well as the
    four coefficients representing population growth and decline,
    calculate the time derivatives dN_1/dt and dN_2/dt and return to the caller.
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
    dN1dt = a*N[0] - b*N[0]*N[1]
    dN2dt = -c*N[1] + d*N[1]*N[0]
    return dN1dt, dN2dt

def euler_solve(func, N1_init=0.3, N2_init=0.6, dt=1.0, t_final=100.0):
    """
    Euler method for solving first-order ODEs (N1, N2).
    
    Parameters
    ----------
    func : function
        A python function that takes `time` and `[N1, N2]` as inputs, and
        returns the time derivatives dN1/dt and dN2/dt as a tuple or list.
    N1_init : float, optional
        Initial value of N1 (default is 0.5).
    N2_init : float, optional
        Initial value of N2 (default is 0.5).
    dt : float, optional
        Time step for the Euler method (default is 0.1).
    t_final : float, optional
        The final time that equations will be solved (default is 100.0).
        
    Returns
    -------
    time : numpy.ndarray
        Array of time values.
    N1 : numpy.ndarray
        Array of N1 values over time.
    N2 : numpy.ndarray
        Array of N2 values over time.
    """
    
    # Initialize time array
    time = np.arange(0, t_final, dt)
    
    # Initialize arrays for N1 and N2
    N1 = np.zeros(time.size)
    N2 = np.zeros(time.size)
    
    # Set initial conditions
    N1[0] = N1_init
    N2[0] = N2_init
    
    # Loop over each time step and update N1 and N2
    for i in range(1, time.size):
        dN1, dN2 = func(time[i-1], [N1[i-1], N2[i-1]])
        N1[i] = N1[i-1] + dt * dN1
        N2[i] = N2[i-1] + dt * dN2
    
    # Return the time and the solutions for N1 and N2
    return time, N1, N2

def solve_rk8(func, N1_init=0.3, N2_init=0.6, dT=10.0, t_final=100.0,
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
    time, N1, N2 = result.t, result.y[0, :], result.y[1, :]
   
    # Return values to caller.
    return time, N1, N2

def plot_q1_q2_results(a=1, b=1, c=1, d=1):
    '''
    This function gathers all the data from the Euler method and RK8 calculation.
    After that, it can create plots that shows the answers.


    Parameters
    ----------
    a, b, c, d : float, defaults=1, 2, 1, 3
    Represent the value of the Lotka-Volterra coefficients. In this function, 
    they can also help to print the coefficients at teh bottom of the graph. 

    '''

    # Create a big plot to contain two subplots 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
    # Euler and RK8 solutions for Lotka-Volterra Competition model
    
    # Getting results for Competition Model
    t_comp_euler, N1_comp_euler, N2_comp_euler = euler_solve (dNdt_comp,
                                                              N1_init=0.3, 
                                                              N2_init=0.6, 
                                                              dt=0.1)
    t_comp_rk8, N1_comp_rk8, N2_comp_rk8 = solve_rk8(dNdt_comp)

    # Getting results for Predator-Prey Model
    t_prey_euler, N1_prey_euler, N2_prey_euler = euler_solve(dNdt_prey, dt=0.05)
    t_prey_rk8, N1_prey_rk8, N2_prey_rk8 = solve_rk8(dNdt_prey)

    # Plotting for Competition Model
    ax1.plot(t_comp_euler, N1_comp_euler, label='N1 Euler', 
             color='red', linewidth=3.2)
    ax1.plot(t_comp_euler, N2_comp_euler, label='N2 Euler', 
             color='skyblue', linewidth=3)
    ax1.plot(t_comp_rk8, N1_comp_rk8, ':', label='N1 RK8', 
             color='orange', linewidth=3)
    ax1.plot(t_comp_rk8, N2_comp_rk8, ':',label='N2 RK8', 
             color='green', linewidth=3)
    
    ax1.set_ylabel('Normalized Population/Carrying Capacity')
    ax1.set_xlabel('Time (years)')
    ax1.set_title('Lotka-Volterra Competition Model')
    ax1.legend()

    # Plotting for Predator-Prey Model
    ax2.plot(t_prey_euler, N1_prey_euler, label='N1 Euler', 
             color='red', linewidth=3)
    ax2.plot(t_prey_euler, N2_prey_euler, label='N2 Euler', 
             color='skyblue', linewidth=3)
    ax2.plot(t_prey_rk8, N1_prey_rk8, ':', label='N1 RK8', 
             color='orange', linewidth=3)
    ax2.plot(t_prey_rk8, N2_prey_rk8, ':',label='N2 RK8', 
             color='green', linewidth=3.5)
    
    ax2.set_ylabel('Normalized Population/Carrying Capacity')
    ax2.set_xlabel('Time (years)')
    ax2.set_title('Lotka-Volterra Predator-Prey Model')
    ax2.legend()

    fig.text(0.5, 0.01, f'Coefficients: a={a}, b={b}, c={c}, d={d}', 
             ha='center', fontsize=12)


    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_q3_results(a=1, b=2, c=1, d=3):
    '''
    This function gathers all the data from the Euler method and RK8 calculation
    for q3.
    After that, it can create plots (verus time diagram and pahse diagram) 
    that show the answers.


    Parameters
    ----------
    a, b, c, d : float, defaults=1, 2, 1, 3
    Represent the value of the Lotka-Volterra coefficients. In this function, 
    they can also help to print the coefficients at teh bottom of the graph. 
    Changing the value here won't change the parameter in the dNdt_comp or
    dNdt_prey functions! 

    '''

    # Create a big plot to contain two subplots 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))


    # Getting results for Predator-Prey Model
    t_prey_euler, N1_prey_euler, N2_prey_euler = euler_solve(dNdt_prey, 
                                                              N1_init=0.3, 
                                                              N2_init=0.6,
                                                              dt=0.05)
    t_prey_rk8, N1_prey_rk8, N2_prey_rk8 = solve_rk8(dNdt_prey)

    # Plotting phase diagram for Predator-Prey Model
    ax1.plot(N1_prey_euler, N2_prey_euler, label='Euler Method', 
             color='blue', linewidth=2)
    ax1.plot(N1_prey_rk8, N2_prey_rk8, ':', label='RK8 Method', 
             color='orange', linewidth=2)
    ax1.set_xlabel('Prey Population (N1)')
    ax1.set_ylabel('Predator Population (N2)')
    ax1.set_title('Phase Diagram: Prey vs. Predator Population')
    ax1.legend()

    # Plotting density vs time diagram for Predator-Prey Model
    ax2.plot(t_prey_euler, N1_prey_euler, label='N1 Euler', 
             color='red', linewidth=3)
    ax2.plot(t_prey_euler, N2_prey_euler, label='N2 Euler', 
             color='skyblue', linewidth=3)
    ax2.plot(t_prey_rk8, N1_prey_rk8, ':', label='N1 RK8', 
             color='orange', linewidth=3)
    ax2.plot(t_prey_rk8, N2_prey_rk8, ':',label='N2 RK8', 
             color='green', linewidth=3.5)
    
    ax2.set_ylabel('Normalized Population/Carrying Capacity')
    ax2.set_xlabel('Time (years)')
    ax2.set_title('Lotka-Volterra Predator-Prey Model')
    ax2.legend()



    fig.text(0.5, 0.01, f'Coefficients: a={a}, b={b}, c={c}, d={d}', 
             ha='center', fontsize=12)


    # Show the plot
    plt.tight_layout()
    plt.show()