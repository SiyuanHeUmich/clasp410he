#!/usr/bin/env python3
'''
A set of tools and routines for solving the N-layer atmosphere energy
balance problem and perform some useful analysis.

To get solution for lab 3, run these commands:

>>> run lab3.py
>>> plot.ion()
>>> q3_experiment1()
>>> q3_experiment2()
>>> q4()
>>> q5()
'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

#  Define some useful constants here.
sigma = 5.67E-8  # Steffan-Boltzman constant.

def n_layer_atmos(N, epsilon=1, S0=1350, albedo=0.33, debug=False):
    '''
    Solve the n-layer atmosphere problem and return temperature at each layer.

    Parameters
    ----------
    N : int
        Set the number of layers.
    epsilon : float, default=1.0
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0.33
        Set the planetary albedo from 0 to 1.
    S0 : float, default=1350
        Set the incoming solar shortwave flux in Watts/m^2.
    debug : boolean, default=False
        Turn on debug output.

    Returns
    -------
    temps : Numpy array of size N+1
        Array of temperatures from the Earth's surface (element 0) through
        each of the N atmospheric layers, from lowest to highest.
    '''

    # Create matrices:
    A = np.zeros([N+1, N+1])
    b = np.zeros(N+1)
    b[0] = -S0/4 * (1-albedo)

    if debug:
        print(f"Populating N+1 x N+1 matrix (N = {N})")

    # Populate our A matrix piece-by-piece.
    for i in range(N+1):
        for j in range(N+1):
            if debug:
                print(f"Calculating point i={i}, j={j}")
            # Diagonal elements are always -2 ('cept at Earth's surface.)
            if i == j:
                A[i, j] = -1*(i > 0) - 1
                # print(f"Result: A[i,j]={A[i,j]}")
            else:
                # Calculate the incoming radiation to the N layer
                m = np.abs(j-i) - 1
                A[i, j] = epsilon * (1-epsilon)**m
    # At Earth's surface, epsilon =1, breaking our pattern.
    # Divide by epsilon along surface to get correct results.
    A[0, 1:] /= epsilon

    # Verify our A matrix.
    if debug:
        print(A)

    # Get the inverse of our A matrix.
    Ainv = np.linalg.inv(A)

    # Multiply Ainv by b to get Fluxes.
    fluxes = np.matmul(Ainv, b)

    # Convert fluxes to temperatures.
    # Fluxes for all atmospheric layers
    temps = (fluxes/epsilon/sigma)**0.25
    temps[0] = (fluxes[0]/sigma)**0.25  # Flux at ground: epsilon=1.
    
    return temps


def n_layer_atmos_q5(N=5, epsilon=0.5, S0=1350, albedo=0, debug=False):
    '''
    Solve the nuclear winter problem and return temperature at each layer.

    Parameters
    ----------
    N : int, default=5
        Set the number of layers.
    epsilon : float, default=0.5
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0
        Set the planetary albedo from 0 to 1. Here, the top layer absorbs all 
        radiation so albedo = 0
    S0 : float, default=1350
        Set the incoming solar shortwave flux in Watts/m^2.
    debug : boolean, default=False
        Turn on debug output.

    Returns
    -------
    temps : Numpy array of size N+1
        Array of temperatures from the Earth's surface (element 0) through
        each of the N atmospheric layers, from lowest to highest.
    '''

    # Create matrices:
    A = np.zeros([N+1, N+1])
    b = np.zeros(N+1)

    # Since the top layer absord all the radiation, 
    # the top layer would be the first layer
    b[N] = -S0/4 * (1-albedo)

    if debug:
        print(f"Populating N+1 x N+1 matrix (N = {N})")

    # Populate our A matrix piece-by-piece.
    for i in range(N+1):
        for j in range(N+1):
            if debug:
                print(f"Calculating point i={i}, j={j}")
            # Diagonal elements are always -2 ('cept at Earth's surface.)
            if i == j:
                A[i, j] = -1*(i > 0) - 1
                # print(f"Result: A[i,j]={A[i,j]}")
            else:
                # Calculate the incoming radiation to the N layer
                m = np.abs(j-i) - 1
                A[i, j] = epsilon * (1-epsilon)**m
    # At Earth's surface, epsilon =1, breaking our pattern.
    # Divide by epsilon along surface to get correct results.
    A[0, 1:] /= epsilon

    # Verify our A matrix.
    if debug:
        print(A)

    # Get the inverse of our A matrix.
    Ainv = np.linalg.inv(A)

    # Multiply Ainv by b to get Fluxes.
    fluxes = np.matmul(Ainv, b)

    # Convert fluxes to temperatures.
    # Fluxes for all atmospheric layers
    temps = (fluxes/epsilon/sigma)**0.25
    temps[0] = (fluxes[0]/sigma)**0.25  # Flux at ground: epsilon=1.
    
    return temps


def q3_experiment1():
    '''
    Solve the first experiment in question 3, and plot surface temperature under
    the different emissivities.

    Parameters
    ----------
    N : int
        Set the number of layers.
    epsilon : float, default=1.0
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0.33
        Set the planetary albedo from 0 to 1.
    S0 : float, default=1350
        Set the incoming solar shortwave flux in Watts/m^2.
    '''

    # Set up steps for emissivity
    emissivities = np.arange(0.1, 1.0, 0.01)
    
    # Create an empty list to store the value from n_layer_atmos
    surface_temps=[]

    for epsilon in emissivities:
        temps = n_layer_atmos(N=1, epsilon=epsilon, S0=1350, albedo=0.33)
        
        # Only need to plot the themperature of the ground layer 
        surface_temps.append(temps[0])
        
        print(f"Emissivity: {epsilon:.2f}, Surface Temperatures: {temps[0]:.3f}K")
        
    
    # Plot Surface Temperature vs Emissivity
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(emissivities, surface_temps, color='r', marker='o', linewidth=3, 
            label='Surface Temperature')
    ax.set_xlabel('Emissivity')
    ax.set_ylabel('Surface Temperature (K)')
    ax.set_title('Surface Temperature vs Emissivity (Single Layer Atmosphere)')
    ax.grid(True)
    ax.legend()
    plt.show()


def q3_experiment2():
    '''
    Solve the second experiment in question 3, and plot surface temperature under
    the different number of atmosphere layers.

    Parameters
    ----------
    N : int
        Set the number of layers.
    epsilon : float, default=1.0
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0.33
        Set the planetary albedo from 0 to 1.
    S0 : float, default=1350
        Set the incoming solar shortwave flux in Watts/m^2.
    '''
    # Set up steps for the number of layers
    N_values = np.arange(1, 10, 1)
    surface_temps = []

    for N in N_values:
        temps = n_layer_atmos(N=N, epsilon=0.255, S0=1350, albedo=0.33)
        surface_temps.append(temps[0])
        print(f"Number of layers: {N:.0f}, Surface Temperatures: {temps[0]:.3f}K")

    # Plot Surface Temperature vs Number of Layers
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(N_values, surface_temps, color='r', marker='o', linewidth=3, 
            label='Surface Temperature')
    ax.set_xlabel('Number of Layers')
    ax.set_ylabel('Surface Temperature (K)')
    ax.set_title('Surface Temperature vs Number of Layers (Emissivity = 0.255)')
    ax.grid(True)
    ax.legend()
    plt.show()


def q4():
    '''
    Solve question 4 about the number of layers on Venus, 
    and print surface temperature under the different number of atmosphere layers.
    Therefore, I can find the number of layer that has the surface temperture of 
    700K
    
    Parameters
    ----------
    N : int
        Set the number of layers.
    epsilon : float, default=1.0
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0.76
        Set the planetary albedo from 0 to 1.
    S0 : float, default=2600
        Set the incoming solar shortwave flux in Watts/m^2.
    '''
    # Perfectly absorbing layers
    epsilon = 1.0  

    # Solar flux for Venus
    S0 = 2600 

    # Venus Bond albedo, from this paper: 
    # https://www.sciencedirect.com/science/article/pii/S0019103516001263?via%3Dihub  
    albedo = 0.76  
    
    # Try N from 75 to 95 and print the result
    N_values=np.arange(75, 95, 1)
    for N in N_values:  
        temps = n_layer_atmos(N, epsilon=epsilon, S0=S0, albedo=albedo)
        surface_temp = temps[0]
        print(f" When there are {N} layers, Surface Temperature = {surface_temp:f} K")


def q5():
    '''
    Solve question 5 about nuclear winter, and plot temperature of each layer
    under the nuclear winter scenario. 

    Parameters
    ----------
    N : int, default=5
        Set the number of layers.
    epsilon : float, default=0.5
        Set the emissivity of the atmospheric layers.
    albedo : float, default=0
        Set the planetary albedo from 0 to 1. Here, the top layer absorbs all 
        radiation so albedo = 0
    S0 : float, default=1350
        Set the incoming solar shortwave flux in Watts/m^2.
    '''
    N = 5
    S0 = 1350
    epsilon = 0.5
    albedo = 0.33

    # Call the nuclear winter model with the nuclear winter parameters
    temperatures = n_layer_atmos_q5(N=5, epsilon=0.5, S0=1350, albedo=0, 
                                    debug=False)

    # Print surface temperature
    print(f"Surface Temperature in Nuclear Winter: {temperatures[0]} K")

    # Plot the temperature profile
    altitude = np.arange(N+1) 
    
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(temperatures, altitude, marker='o', color='r', linewidth=3, 
            label='Temperature at Each Layer')
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Altitude (layers)')
    ax.set_title('Temperature vs Altitude in Nuclear Winter')
    ax.grid(True)
    ax.legend()
    plt.show()