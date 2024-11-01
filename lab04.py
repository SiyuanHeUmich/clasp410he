#1/usr/bin/env python3
'''
This file contain tools and methods for solving our heat equation and heat 
diffusion

To get solution for lab 4: Run these commands:
>>> run lab04.py
>>> import matplotlib.pyplot as plt
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

solution = np.array(([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0.64, 0.48, 0.40, 0.32, 0.26, 0.21, 0.17, 0.1375, 0.11125, 0.090000, 0.072812],
                    [0.96, 0.80, 0.64, 0.52, 0.42, 0.34, 0.28, 0.2225, 0.18000, 0.145625, 0.117813],
                    [0.96, 0.80, 0.64, 0.52, 0.42, 0.34, 0.28, 0.2225, 0.18000, 0.145625, 0.117813],
                    [0.64, 0.48, 0.40, 0.32, 0.26, 0.21, 0.17, 0.1375, 0.11125, 0.090000, 0.072812],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))

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


# create heat diff function that works for both problems:
def heatdiff(xmax, tmax, dx, dt, c2, debug=False, conditions='permafrost', temp_shift=0):
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

    if conditions == 'wire':
        # Set initial conditions:
        U[:, 0] = 4*xgrid - 4*xgrid**2

        # Set boundary conditions:
        U[0, :] = 0
        U[-1, :] = 0

    elif conditions == 'permafrost':
        # Set boundary conditions:
        U[0, :] = temp_kanger(tgrid, temp_shift=temp_shift)
        U[-1, :] = 5

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
    time, x, U = heatdiff(xmax=1, tmax=0.2, dx=0.2, dt=0.02, c2=1, conditions='wire')
    diff = solution - U
    print(f'your solution: {U}')
    print(f'difference between expected and calculated solution: {diff}')

def plot_wire_heat_map():
    '''
    Plot a heat map of temperature distribution over length and 
    time for a heat diffusion problem in a wire. 
    Displays a heat map showing temperature changes over depth and time.
    '''
    # Get solution using your solver:
    time, x, heat = heatdiff(xmax=1, tmax=0.2, dx=0.2, dt=0.02, c2=1, conditions='wire')
    
    # Create a figure/axes object
    fig, axes = plt.subplots(1, 1)
    # Create a color map and add a color bar.
    map = axes.pcolor(time, x, heat, cmap='seismic', vmin=0, vmax=1, shading='nearest')
    plt.colorbar(map, ax=axes, label='Temperature (C)')
    plt.xlabel("Time (days)")
    plt.ylabel("Depth (m)")
    plt.title("Temperature Distribution Over Depth and Time")
    plt.show()

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
    time, xgrid, U, dt = heatdiff(warming_shift=temp_shift)

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
    plt.show() 

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

