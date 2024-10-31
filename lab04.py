#1/usr/bin/env python3
'''
This file contain tools and methods for solving our heat equation and heat 
diffusion

To get solution for lab 4: Run these commands:
>>> run diffusion1.py
>>> plt.ion()
>>> validate_solver()
>>> plot_heat_map()
>>> plot_heat_map_with_temp_shift(0)
>>> plot_heat_map_with_temp_shift(0.5)
>>> plot_heat_map_with_temp_shift(1)
>>> plot_heat_map_with_temp_shift(3)
>>> plot_ground_profile_map_all_scenarios()

'''

import numpy as np
import matplotlib.pyplot as plt

# Kangerlussuaq average temperature:
t_kanger = np.array([-19.7, -21.0, -17., -8.4, 2.3, 8.4,10.7, 8.5, 
                     3.1, -6.0, -12.0, -16.9])

def heatdiff(xmax=1, tmax=0.2, dx=0.1, dt=0.002, c2=1, debug=True):
    '''
    Solve the heat equation by using a forward-difference method

    Parameters
    ----------
    xmax : float, default=1
        Maximum spatial extent, unit: meters
    tmax : float, default=0.2
        Length of time for simulation, unit: seconds
    dx : float, default=0.1
        Spatial step size, unit: meters
    dt : float, default=0.002
        Time step size, unit: seconds
    c2 : float, default=1
        Thermal diffusivity constant
    debug : boolean, default=True
        Print additional information about the grid and stability.

    Returns
    -------
    tgrid : array
        Array of time used in the simulation.
    xgrid : array
        Array of depth used in the simulation.
    U : array
        2D array of temperatures at each depth and time step.

    '''
    if dt>(dx**2)/(2*c2):
        raise ValueError('dt is too large! Must be less than dx^2/ (2*c2) for stability')
    
    # Start by calculating size of array: MxN
    M = int(np.round(xmax / dx + 1))
    N = int(np.round(tmax / dt + 1))

    xgrid, tgrid = np.arange(0, xmax+dx, dx), np.arange(0, tmax+dt, dt)

    if debug:
        print(f'Our grid goes from 0 to {xmax}m and 0 to {tmax}s')
        print(f'Our spatial step is {dx} and time step is {dt}')
        print(f'There are {M} points in space and {N} points in time.')
        print('Here is our spatial grid:')
        print(xgrid)
        print('Here is our time grid:')
        print(tgrid)

    # Initialize our data array:
    U = np.zeros([M, N])

    # Set initial conditions:
    U[:, 0] = 4*xgrid - 4*xgrid**2

    # Set boundary conditions:
    U[0, :] = 0
    U[-1, :] = 0

    # Set our "r" constant.
    r = c2 * dt / dx**2

    # Solve forward differnce
    for j in range(N-1):
        U[1:-1, j+1] = (1-2*r) * U[1:-1, j] + \
            r*(U[2:, j] + U[:-2, j])

    # Return grid and result
    return tgrid, xgrid, U

def validate_solver():
    '''
    Validate the heat equation solver by comparing a test case provided in the handout. 

    Returns
    -------
    time : array
        Array of time step for the validation test case.
    x : array
        Array of depth for the validation test case.
    U : array
        2D array of temperatures for each depth and time step.
    '''
    #The first array is the time. The second array is the depth. The 2D array is
    # the temperture at certain depth.  
    time, x, U = heatdiff(xmax=1, tmax=0.2, dx=0.2, dt=0.02, c2=1, debug=False)
    
    return time, x, U

def plot_heat_map():
    '''
    Plot a heat map of temperature distribution over depth and 
    time for a heat diffusion problem. 
    Displays a heat map showing temperature changes over depth and time.
    '''
    # Get solution using your solver:
    time, x, heat = heatdiff()
    
    # Create a figure/axes object
    fig, axes = plt.subplots(1, 1)
    # Create a color map and add a color bar.
    map = axes.pcolor(time, x, heat, cmap='seismic', vmin=0, vmax=1, shading='nearest')
    plt.colorbar(map, ax=axes, label='Temperature (C)')
    plt.xlabel("Time (days)")
    plt.ylabel("Depth (m)")
    plt.title("Temperature Distribution Over Depth and Time")
    plt.show()

def temp_kanger(t, temp_shift=0):
    '''
    For an array of times in days, return timeseries of temperature for
    Kangerlussuaq, Greenland.

    Parameters
    ----------
    t : float
        Time in days for  the temperature.
    temp_shift : float, default=0
        Temperature shift to simulate warming effects in degrees Celsius.
    Returns
    -------
    float
        Surface temperature for Kangerlussuaq at time t.
    '''
    t_amp = (t_kanger - t_kanger.mean()).max()
    return t_amp*np.sin(np.pi/180 * t - np.pi/2) + t_kanger.mean()+ temp_shift

def permafrost_heatdiff_model(dx=1, dt=60*60*24, c2=0.25e-6, years=100, warming_shift=0):
    '''
    Solve the heat equation by using a forward-difference method for permafrost problem

    Parameters
    ----------
    xmax : float, default=100
        Maximum spatial extent, unit: meters
    tmax : float
        Length of time for simulation, unit: seconds
    dx : float, default=1
        Spatial step size, unit: meters
    dt : float, default=60*60*24
        Time step size, unit: seconds
    c2 : float, default=0.25e-6
        Thermal diffusivity constant

    Returns
    -------
    tgrid : array
        Array of time used in the simulation.
    xgrid : array
        Array of depth used in the simulation.
    U : array
        2D array of temperatures at each depth and time step.

    '''
    
    # Set 100m depth, 100-year simulation
    xmax = 100 
    tmax= years*365*24*60*60
    
    time = np.arange(0, tmax + dt, dt)
    xgrid = np.arange(0, xmax + dx, dx)
    
    # Start by calculating size of array: MxN
    M = int(np.round(xmax / dx + 1))
    N = int(np.round(tmax / dt + 1))
    
    #Initialize the values in the array
    U = np.zeros([M, N])
    
    # Initial ground temperature set to 0
    U[:, 0] = 0  

    # Set Upper boundary
    U[0, :] = [temp_kanger(t/86400, temp_shift=warming_shift) for t in time] 
    # Set Lower boundary at 5 degrees
    U[-1, :] = 5  

    # Set our "r" constant.
    r = c2 * dt / dx**2

    # Solve! Forward differnce ahoy.
    for j in range(N-1):
        U[1:-1, j+1] = (1-2*r) * U[1:-1, j] + \
            r*(U[2:, j] + U[:-2, j])

    return time, xgrid, U, dt

def plot_heat_map_with_temp_shift(temp_shift=0):
    '''
    Plot a heat map for the Kangerlussuaq permafrost temperature distribution 
    with time with a temperature shift. Displays a heat map of temperature 
    distribution over depth and time with the specified temperature shift.

    Parameters
    ----------
    temp_shift : float, default=0
        Temperature shift in degrees Celsius to simulate warming effects

    Return
    ----------
    The plots of heat map
    '''
    # Run permafrost model with the specified temperature shift
    # Collect data from the results
    time, xgrid, U, dt = permafrost_heatdiff_model(warming_shift=temp_shift)

    # Convert time, seconds, to years for x-axis labeling
    tgrid_years = time / (365 * 86400)

    # Plotting the heat map
    fig, ax = plt.subplots(figsize=(10, 6))
    heat_map = ax.pcolor(tgrid_years, xgrid, U, cmap='seismic', vmin=-25, vmax=25, shading='nearest')
    ax.invert_yaxis() 

    plt.colorbar(heat_map, ax=ax, label='Temperature (C)')
    plt.xlabel("Time (years)")
    plt.ylabel("Depth (m)")
    plt.title(f"Temperature Distribution with {temp_shift}C Shift in Kangerlussuaq, Greenland")    

def plot_ground_profile_map(U, xgrid, dt, temp_shift):
    '''
    Plot seasonal temperature for winter and summer at various depths

    Parameters
    ----------
    U : array
        2D array of temperatures at each depth and time step.
    xgrid : array
        Array of depths in meters.
    dt : float
        Time step size in seconds.
    temp_shift : float
        Temperature shift in degrees Celsius applied to show warming effects.

    '''

    # Set indexing for the final year of results:
    loc = int(-365 * 86400 / dt) 
    # Extract the min values for winter over the final year:
    winter = U[:, loc:].min(axis=1)
    # Extract the max values for summer over the final year:
    summer = U[:, loc:].max(axis=1)

    # Create a temperature data plot:
    fig, ax2 = plt.subplots(1, 1, figsize=(10, 8))
    ax2.plot(winter, xgrid, label='Winter', color='blue')
    ax2.plot(summer, xgrid, linestyle='--', label='Summer', color='red')
    ax2.invert_yaxis()

    plt.xlim(-15, 8)
    plt.xlabel('Temperature (C)')
    plt.ylabel('Depth (m)')
    plt.legend()
    plt.title(f"Seasonal Temperature Profile with Warming Shift {temp_shift}C")
    plt.grid(True)
    plt.show()

def plot_ground_profile_map_all_scenarios():
    '''
    Run the permafrost model under multiple scenarios and plot ground profile maps.

    This function uses a range of warming shifts to simulate increased surface 
    temperatures and plots the resulting seasonal temperature profiles 

    '''
    # Set a list of temperature for different scenarios
    for temp_shift in [0, 0.5, 1, 3]: 
        # Recording data from permafrost_heatdiff_model and use them in plot_ground_profile_map
        time, xgrid, U, dt = permafrost_heatdiff_model(warming_shift=temp_shift)
        plot_ground_profile_map(U, xgrid, dt=86400, temp_shift=temp_shift)