#!/usr/bin/env python3

'''
This file contain tools and methods for solving code for Lab 5: SNOWBALL EARTH!!!

To get solution for lab 5: Run these commands:
>>> run lab05.py
>>> test_snowball()
>>> q2_plot()
>>> q3_plot()
>>> q4_plot()
'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

# Some constants:
radearth = 6357000.  # Earth radius in meters.
mxdlyr = 50.         # depth of mixed layer (m)
sigma = 5.67e-8      # Steffan-Boltzman constant
C = 4.2e6            # Heat capacity of water
rho = 1020           # Density of sea-water (kg/m^3)


def gen_grid(nbins=18):
    '''
    Generate a grid from 0 to 180 lat (where 0 is south pole, 180 is north)
    where each returned point represents the cell center.

    Parameters
    ----------
    nbins : int, defaults to 18
        Set the number of latitude bins.

    Returns
    -------
    dlat : float
        Grid spacing in degrees.
    lats : Numpy array
        Array of cell center latitudes.
    '''

    dlat = 180 / nbins  # Latitude spacing.
    lats = np.arange(0, 180, dlat) + dlat/2.

    # Alternative way to obtain grid:
    # lats = np.linspace(dlat/2., 180-dlat/2, nbins)

    return dlat, lats

def temp_warm(lats_in):
    '''
    Create a temperature profile for modern day "warm" earth.

    Parameters
    ----------
    lats_in : Numpy array
        Array of latitudes in degrees where temperature is required

    Returns
    -------
    temp : Numpy array
        Temperature in Celcius.
    '''

    # Get base grid:
    dlat, lats = gen_grid()

    # Set initial temperature curve:
    T_warm = np.array([-47, -19, -11, 1, 9, 14, 19, 23, 25, 25,
                       23, 19, 14, 9, 1, -11, -19, -47])
    coeffs = np.polyfit(lats, T_warm, 2)

    # Now, return fitting:
    temp = coeffs[2] + coeffs[1]*lats_in + coeffs[0] * lats_in**2

    return temp

def insolation(S0, lats):
    '''
    Given a solar constant (`S0`), calculate average annual, longitude-averaged
    insolation values as a function of latitude.
    Insolation is returned at position `lats` in units of W/m^2.

    Parameters
    ----------
    S0 : float
        Solar constant (1370 for typical Earth conditions.)
    lats : Numpy array
        Latitudes to output insolation. Following the grid standards set in
        the diffusion program, polar angle is defined from the south pole.
        In other words, 0 is the south pole, 180 the north.

    Returns
    -------
    insolation : numpy array
        Insolation returned over the input latitudes.
    '''

    # Constants:
    max_tilt = 23.5   # tilt of earth in degrees

    # Create an array to hold insolation:
    insolation = np.zeros(lats.size)

    #  Daily rotation of earth reduces solar constant by distributing the sun
    #  energy all along a zonal band
    dlong = 0.01  # Use 1/100 of a degree in summing over latitudes
    angle = np.cos(np.pi/180. * np.arange(0, 360, dlong))
    angle[angle < 0] = 0
    total_solar = S0 * angle.sum()
    S0_avg = total_solar / (360/dlong)

    # Accumulate normalized insolation through a year.
    # Start with the spin axis tilt for every day in 1 year:
    tilt = [max_tilt * np.cos(2.0*np.pi*day/365) for day in range(365)]

    # Apply to each latitude zone:
    for i, lat in enumerate(lats):
        # Get solar zenith; do not let it go past 180. Convert to latitude.
        zen = lat - 90. + tilt
        zen[zen > 90] = 90
        # Use zenith angle to calculate insolation as function of latitude.
        insolation[i] = S0_avg * np.sum(np.cos(np.pi/180. * zen)) / 365.

    # Average over entire year; multiply by S0 amplitude:
    insolation = S0_avg * insolation / 365

    return insolation

def snowball_earth(nbins=18, dt=1., tstop=10000, lam=100., spherecorr=True, debug=False, 
                   albedo=0.3, albedo_gnd=0.3, albedo_ice=0.6, emiss=1, S0=1370, gamma=1, 
                   hot_earth=False, cold_earth=False, da_earth=False, ff_earth=False):
    '''
    Perform snowball earth simulation.

    Parameters
    ----------
    nbins : int, defaults to 18
        Number of latitude bins.
    dt : float, defaults to 1
        Timestep in units of years
    tstop : float, defaults to 10,000
        Stop time in years
    lam : float, defaults to 100
        Diffusion coefficient of ocean in m^2/s
    spherecorr : bool, defaults to True
        Use the spherical coordinate correction term. This should always be
        true except for testing purposes.
    debug : bool, defaults to False
        Turn  on or off debug print statements.
    albedo : float, defaults to 0.3
        Set the Earth's albedo.
    emiss : float, defaults to 1.0
        Set ground emissivity. Set to zero to turn off radiative cooling.
    S0 : float, defaults to 1370
        Set incoming solar forcing constant. Change to zero to turn off
        insolation.

    Returns
    -------
    lats : Numpy array
        Latitude grid in degrees where 0 is the south pole.
    Temp : Numpy array
        Final temperature as a function of latitude.
    '''
    # Get time step in seconds:
    dt_sec = 365 * 24 * 3600 * dt  # Years to seconds.

    # Generate grid:
    dlat, lats = gen_grid(nbins)

    # Get grid spacing in meters.
    dy = radearth * np.pi * dlat / 180.

    # Generate insolation:
    insol = gamma*insolation(S0, lats)

    # Create initial condition:
    if hot_earth:
        Temp=60*np.ones(nbins)
    elif cold_earth:
        Temp=-60*np.ones(nbins)
    else:
        Temp = temp_warm(lats)
        if debug:
            print('Initial temp = ', Temp)

    # Get number of timesteps:
    nstep = int(tstop / dt)

    # Debug for problem initialization
    if debug:
        print("DEBUG MODE!")
        print(f"Function called for nbins={nbins}, dt={dt}, tstop={tstop}")
        print(f"This results in nstep={nstep} time step")
        print(f"dlat={dlat} (deg); dy = {dy} (m)")
        print("Resulting Lat Grid:")
        print(lats)

    # Build A matrix:
    if debug:
        print('Building A matrix...')
    A = np.identity(nbins) * -2  # Set diagonal elements to -2
    A[np.arange(nbins-1), np.arange(nbins-1)+1] = 1  # Set off-diag elements
    A[np.arange(nbins-1)+1, np.arange(nbins-1)] = 1  # Set off-diag elements
    # Set boundary conditions:
    A[0, 1], A[-1, -2] = 2, 2

    # Build "B" matrix for applying spherical correction:
    B = np.zeros((nbins, nbins))
    B[np.arange(nbins-1), np.arange(nbins-1)+1] = 1  # Set off-diag elements
    B[np.arange(nbins-1)+1, np.arange(nbins-1)] = -1  # Set off-diag elements
    # Set boundary conditions:
    B[0, :], B[-1, :] = 0, 0

    # Set the surface area of the "side" of each latitude ring at bin center.
    Axz = np.pi * ((radearth+50.0)**2 - radearth**2) * np.sin(np.pi/180.*lats)
    dAxz = np.matmul(B, Axz) / (Axz * 4 * dy**2)

    if debug:
        print('A = ', A)
    # Set units of A derp
    A /= dy**2

    # Get our "L" matrix:
    L = np.identity(nbins) - dt_sec * lam * A
    L_inv = np.linalg.inv(L)


    if debug:
        print('Time integrating...')
    for i in range(nstep):

        if hot_earth or cold_earth or da_earth:
            # Add spherical correction term:
            if spherecorr:
                Temp += dt_sec * lam * dAxz * np.matmul(B, Temp)


            # Update albedo dynamically based on temperature at each latitude
            loc_ice = Temp <= -10  # Mask where temperature is <= -10°C
            albedo = np.zeros_like(Temp)  # Create an empty array with the same shape as Temp

            # Set albedo values based on the temperature condition
            albedo[loc_ice] = albedo_ice  # Set albedo to ice for temperatures <= -10°C
            albedo[~loc_ice] = albedo_gnd  # Set albedo to ground/water for temperatures > -10°C


            # Apply insolation and radiative losses:
            # print('T before insolation:', Temp)
            radiative = (1-albedo) * insol - emiss*sigma*(Temp+273.15)**4
            # print('\t Rad term = ', dt_sec * radiative / (rho*C*mxdlyr))
            Temp += dt_sec * radiative / (rho*C*mxdlyr)
            # print('\t T after rad:', Temp)

            Temp = np.matmul(L_inv, Temp)            
        elif ff_earth:
            if spherecorr:
                Temp += dt_sec * lam * dAxz * np.matmul(B, Temp)

            # Apply insolation and radiative losses:
            # print('T before insolation:', Temp)
            radiative = (1-albedo_ice) * insol - emiss*sigma*(Temp+273.15)**4
            # print('\t Rad term = ', dt_sec * radiative / (rho*C*mxdlyr))
            Temp += dt_sec * radiative / (rho*C*mxdlyr)
            # print('\t T after rad:', Temp)

            Temp = np.matmul(L_inv, Temp)
        else:
            # Add spherical correction term:
            if spherecorr:
                Temp += dt_sec * lam * dAxz * np.matmul(B, Temp)

        # Apply insolation and radiative losses:
        # print('T before insolation:', Temp)
        radiative = (1-albedo) * insol - emiss*sigma*(Temp+273.15)**4
        # print('\t Rad term = ', dt_sec * radiative / (rho*C*mxdlyr))
        Temp += dt_sec * radiative / (rho*C*mxdlyr)
        # print('\t T after rad:', Temp)

        Temp = np.matmul(L_inv, Temp)
            

    return lats, Temp

def test_snowball(tstop=10000):
    '''
    Reproduce example plot in handout.

    Using our DEFAULT values (grid size, diffusion, etc.) and a warm-Earth
    initial condition, plot:
        - Initial condition
        - Plot simple diffusion only
        - Plot simple diffusion + spherical correction
        - Plot simple diff + sphere corr + insolation
    '''
    nbins = 18

    # Generate very simple grid
    # Generate grid:
    dlat, lats = gen_grid(nbins)

    # Create initial condition:
    initial = temp_warm(lats)

    # Get simple diffusion solution:
    lats, t_diff = snowball_earth(tstop=tstop, spherecorr=False, S0=0, emiss=0)

    # Get diffusion + spherical correction:
    lats, t_sphe = snowball_earth(tstop=tstop, S0=0, emiss=0)

    # Get diffusion + sphercorr + radiative terms:
    lats, t_rad = snowball_earth(tstop=tstop)

    # Create figure and plot!
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(lats, initial, label='Warm Earth Init. Cond.')
    ax.plot(lats, t_diff, label='Simple Diffusion')
    ax.plot(lats, t_sphe, label='Diffusion + Sphere. Corr.')
    ax.plot(lats, t_rad, label='Diffusion + Sphere. Corr. + Radiative')

    ax.set_xlabel('Latitude (0=South Pole)')
    ax.set_ylabel('Temperature ($^{\circ} C$)')
    ax.legend(loc='best')
    fig.tight_layout()

def q2_plot(tstop=10000):
    '''
    Reproduce example plot in  the handout.
    Using Diffusivity=100, Emissivity=0.7to genarate the temperature profile.
    With these parameter, I can get the temperature profile that fits the best
    with temp_warm. 
    '''
    nbins = 18

    # Generate very simple grid
    # Generate grid:
    dlat, lats = gen_grid(nbins)

    # Create initial condition:
    initial = temp_warm(lats)
    

    # Create figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Plot warm-Earth profile for reference
    ax.plot(lats, initial, label='Warm Earth (Reference)')
    

    lats, t_lamdiff = snowball_earth(nbins=18, dt=1., tstop=10000, lam=100, spherecorr=True,
                        debug=False, albedo=0.3, emiss=0.7, S0=1370)
    ax.plot(lats, t_lamdiff, label='Diffusivity=100, Emissivity=0.7', linestyle='--', color='black')

    ax.set_xlabel('Latitude (0=South Pole)')
    ax.set_ylabel('Temperature (degree C)')
    ax.legend(loc='best', fontsize='small', ncol=2)
    ax.set_title('Comparison of Warm Earth to Simulated Equilibria')
    fig.tight_layout()
    plt.show()

def q3_plot(tstop=10000, lam=100, emiss=0.7, albedo_ice=0.6, albedo_gnd=0.3, albedo=0.3):
    '''
    Plot temperature profiles for various scenarios in Q3.

    Parameters
    ----------
    tstop : float
        Total simulation time in years.
    lam : float
        Diffusion coefficient.
    emiss : float
        Emissivity of the Earth's surface.
    albedo_ice : float
        Albedo of ice-covered regions.
    albedo_gnd : float
        Albedo of ground/water.
    albedo : float
        Initial albedo for the static scenario.

    Returns
    '''
    nbins = 18
    dlat, lats = gen_grid(nbins)

    # Simulate Hot Earth
    lats, t_hot = snowball_earth(nbins=nbins, dt=1., tstop=tstop, lam=lam,spherecorr=True,
                                       emiss=emiss, albedo=albedo_gnd, S0=1370, hot_earth=True)

    # Simulate Cold Earth
    lats, t_cold = snowball_earth(nbins=nbins, dt=1., tstop=tstop, lam=lam,spherecorr=True,
                                        emiss=emiss, albedo=albedo_gnd, S0=1370, cold_earth=True)

    # Simulate Flash Freeze
    lats, t_flash = snowball_earth(nbins=nbins, dt=1., tstop=tstop, lam=lam,spherecorr=True,
                                        emiss=emiss, albedo=albedo_ice, S0=1370, ff_earth=True)

    # Simulate Warm Earth with static albedo
    lats, t_warm_origin = snowball_earth(nbins=nbins, dt=1., tstop=tstop, lam=lam,spherecorr=True,
                                        emiss=emiss, albedo=albedo, S0=1370)

    # Simulate Warm Earth with dynamic albedo
    lats, t_warm_dynamic = snowball_earth(nbins=nbins, dt=1., tstop=tstop, lam=lam,
                                        emiss=emiss, albedo=albedo_gnd, S0=1370, da_earth=True)

    # Plot all scenarios
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.plot(lats, t_hot, label="Hot Earth (60°C)", color='red', linestyle='-')
    ax.plot(lats, t_cold, label="Cold Earth (-60°C)", color='blue', 
            linestyle='--', linewidth=4.5)
    ax.plot(lats, t_flash, label="Flash Freeze (High Albedo=0.6)", color='pink')
    ax.plot(lats, t_warm_origin, label="Warm Earth (Static Albedo=0.3)", color='green')
    ax.plot(lats, t_warm_dynamic, label="Warm Earth (Dynamic Albedo)", 
            color='orange', linestyle='--')


    ax.set_xlabel('Latitude (0° = South Pole)', fontsize=12)
    ax.set_ylabel('Temperature (°C)', fontsize=12)
    ax.set_title('Equilibrium Temperature Profiles for Different Scenarios', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True)
    fig.tight_layout()
    plt.show()

def q4_plot(tstop=10000, dt=1., nbins=18, lam=100., emiss=0.7, 
                      albedo_ice=0.6, albedo_gnd=0.3,
                      gamma_start=0.4, gamma_stop=1.4, gamma_step=0.05):
    '''
    Explore the impact of solar forcing on snowball Earth by varying the solar multiplier (gamma),
    and plot the average global temperature versus the solar multiplier.

    Parameters
    ----------
    tstop : float, default=10000
        Total simulation time for each step in years.
    dt : float, default=1
        Time step size in years.
    nbins : int, default=18
        Number of latitude bins.
    lam : float, default=100
        Diffusion coefficient.
    emiss : float, default=0.7
        Emissivity.
    albedo_ice : float, default=0.6
        Albedo for ice-covered regions.
    albedo_gnd : float, default=0.3
        Albedo for ground/water regions.
    gamma_start : float, default=0.4
        Starting value for the solar multiplier gamma.
    gamma_stop : float, default=1.4
        Ending value for the solar multiplier.
    gamma_step : float, default=0.05
        Increment/Decrement step size for gamma.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Generate grid
    dlat, lats = gen_grid(nbins)

    # Initialize storage for results
    gamma_values_inc = []
    avg_temps_inc = []

    gamma_values_dec = []
    avg_temps_dec = []

    # Increasing gamma     
    for gamma in np.arange(gamma_start, gamma_stop + gamma_step, gamma_step):        
        lats, Temp = snowball_earth(nbins=nbins, dt=dt, tstop=tstop, lam=lam, emiss=emiss, 
                                    gamma=gamma,albedo_gnd=albedo_gnd, albedo_ice=albedo_ice,
                                    cold_earth=True)
        avg_temp = np.mean(Temp)
        gamma_values_inc.append(gamma)
        avg_temps_inc.append(avg_temp)
        #print(Temp)

    # Decreasing gamma: Start from the final state of increasing gamma, no longer cold_earth
    for gamma in np.arange(gamma_stop, gamma_start - gamma_step, -gamma_step):
        lats, Temp = snowball_earth(nbins=nbins, dt=dt, tstop=tstop, lam=lam, emiss=emiss, 
                                    gamma=gamma,albedo_gnd=albedo_gnd, albedo_ice=albedo_ice,
                                    da_earth=True)
        avg_temp = np.mean(Temp)
        gamma_values_dec.append(gamma)
        avg_temps_dec.append(avg_temp)
        #print(Temp)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot results for increasing gamma
    ax.plot(gamma_values_inc, avg_temps_inc, marker='o', 
             label='Increasing gamma', color='blue')

    # Plot results for decreasing gamma
    ax.plot(gamma_values_dec, avg_temps_dec, marker='o', 
             label='Decreasing gamma', color='red')

    ax.set_xlabel('Solar Multiplier gamma', fontsize=12)
    ax.set_ylabel('Average Global Temperature (degree C))', fontsize=12)
    ax.set_title('Impact of Solar Forcing on Snowball Earth Stability', fontsize=14)
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    plt.show()
