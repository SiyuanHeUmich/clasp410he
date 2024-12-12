#!/usr/bin/env python3

'''
This file contain tools and methods for solving code for Final Project

To get solution for final project Run these commands:
>>> pip install burnman
>>> run final_project.py
>>> q1()
>>> q2()
>>> q2_plot_error_vs_kprime()
>>> q2_plot_best_fit()
>>> q3()

'''
import numpy as np
import matplotlib.pyplot as plt
from burnman.eos.birch_murnaghan import volume

# Solid phase parameters
solid_params = {
    'K_0': 45e9, # Bulk modulus in Pa (solid phase)
    'Kprime_0': 4, # Pressure derivative of bulk modulus (solid phase)
    'V_0': 6.634e-5, # Reference molar volume in m^3/mol (solid phase)
    'P_0': 1e5 # Reference pressure in Pa
}

# Liquid phase parameters
liquid_params = {
    'K_0': 3.8e9,  # Bulk modulus in Pa (liquid phase)
    'Kprime_0': 4, # Pressure derivative of bulk modulus (to be varied)
    'V_0': 7.15e-5,  # Reference molar volume in m^3/mol (liquid phase)
    'P_0': 1e5  # Reference pressure in Pa
}

# Constants
dH_fusion = 27.6e3  # Enthalpy of fusion (J/mol)
T_m0 = 1170  # Melting temperature at reference pressure (K)

# Pressure range (Pa)
pressures = np.linspace(1e5, 2e10, 100)

# Experimental estimates
estimate_p = np.array([3e9, 6e9, 7e9, 10e9]) #Pa
estimate_temp = np.array([1368+273, 1480+273, 1510+273, 1670+273]) #K

def compute_melting_curve(Kprime_0, pressures, solid_params, liquid_params, T_m0, dH_fusion):
    """
    Compute the melting curve as a function of pressure for given phase parameters.

    Parameters:
        Kprime_0 (float): 
            Pressure derivative of bulk modulus for the liquid phase.
        pressures (array): 
            Array of pressure values (Pa).
        solid_params: 
            Call solid phase parameters.
        liquid_params: 
            Call liquid phase parameters.
        T_m0 (float): 
            Melting temperature at reference pressure (K).
        dH_fusion (float): 
            Enthalpy of fusion (J/mol).

    Returns:
        Temperatures corresponding to the pressure values.
    """
    liquid_params['Kprime_0'] = Kprime_0
    temperatures = [T_m0]
    for i in range(1, len(pressures)):
        P_prev, P_curr = pressures[i - 1], pressures[i]
        V_solid = volume(P_curr, solid_params)
        V_liquid = volume(P_curr, liquid_params)
        delta_V = V_liquid - V_solid

        # avoid divide by 0 or a nagetive number
        if abs(delta_V) < 1e-10:
            dT = 0
        else:
            delta_S = dH_fusion / temperatures[-1]
            dT = (P_curr - P_prev) * delta_V / delta_S

        temperatures.append(temperatures[-1] + dT)
    return np.array(temperatures)

def q1():
    """
    Plot melting curves for different values of the pressure derivative (K0')
    of the liquid phase's bulk modulus.
    """
    plt.figure(figsize=(10, 6))
    kprime_values = np.arange(0, 200, 20)
    for kprime in kprime_values:
        temperatures = compute_melting_curve(kprime, pressures, solid_params, 
                                             liquid_params, T_m0, dH_fusion)
        plt.plot(pressures / 1e9, temperatures, label=f"K0' = {kprime}")

    plt.scatter(estimate_p / 1e9, estimate_temp, color='red', label="Estimate Data from Experiments")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Pressure (GPa)')
    plt.ylim(500, 2000)
    plt.title("Melting Curves for Different K' Values of Liquid Phase")
    plt.grid()
    plt.legend()
    plt.show()

def q2():
    """
    Find the best-fit value of K0' by minimizing the error between experimental
    and computed melting temperatures.

    Returns:
        float: Best-fit value of K0'.
    """
    best_kprime = None
    min_error = float('inf')
    kprime_values = np.arange(60, 100, 0.01)
    for kprime in kprime_values:
        model_temperatures = compute_melting_curve(kprime, pressures, solid_params, 
                                                   liquid_params, T_m0, dH_fusion)
        interpolated_temps = np.interp(estimate_p, pressures, model_temperatures)
        error = np.sum((interpolated_temps - estimate_temp) ** 2)
        if error < min_error:
            min_error = error
            best_kprime = kprime
    
    best_temperatures = compute_melting_curve(best_kprime, pressures, solid_params, 
                                              liquid_params, T_m0, dH_fusion)

    plt.figure(figsize=(10, 6))
    plt.plot(pressures / 1e9, best_temperatures, label=f"Best-fit K0' = {best_kprime:.2f}", color='blue')
    plt.scatter(estimate_p / 1e9, estimate_temp, color='red', label="Experimental Data")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Pressure (GPa)')
    plt.title("Best-Fit Melting Curve")
    plt.legend()
    plt.grid()
    plt.show()        
    return best_kprime

def q2_plot_error_vs_kprime():
    """
    Plot the sum of squared errors as a function of K0' values.
    """
    kprime_values = np.arange(60, 100, 0.01) 
    errors = []

    for kprime in kprime_values:
        model_temperatures = compute_melting_curve(kprime, pressures, solid_params, 
                                                   liquid_params, T_m0, dH_fusion)
       
        # Since the pressure value in pressure range different from the estimate 
        # pressure value, I decide to interpolated midpoint data
        interpolated_temps = np.interp(estimate_p, pressures, model_temperatures)
        error = np.sum((interpolated_temps - estimate_temp) ** 2)
        errors.append(error)

    plt.figure(figsize=(10, 6))
    plt.plot(kprime_values, errors, marker='.', color='blue')
    plt.xlabel("K0' Value")
    plt.ylabel("Sum of Squared Errors")
    plt.title("Error vs K0' for Grid Search Method")
    plt.grid()
    plt.show()

def q2_plot_best_fit():
    """
    Plot the best-fit melting curve and compare it with experimental data.
    """
    best_kprime = q2()
    print(f"Best-fit K0' value: {best_kprime:.2f}")
    best_temperatures = compute_melting_curve(best_kprime, pressures, solid_params, 
                                              liquid_params, T_m0, dH_fusion)

    plt.figure(figsize=(10, 6))
    plt.plot(pressures / 1e9, best_temperatures, label=f"Best-fit K0' = {best_kprime:.2f}", color='blue')
    plt.scatter(estimate_p / 1e9, estimate_temp, color='red', label="Experimental Data")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Pressure (GPa)')
    plt.title("Best-Fit Melting Curve")
    plt.legend()
    plt.grid()
    plt.show()

def objective_function(Kprime_0):
    """
    Setup for simulated annealing method
    Calculate the sum of squared errors between experimental and computed melting temperatures.

    Parameters:
        Kprime_0 (float): Pressure derivative of bulk modulus for the liquid phase.

    Returns:
        float: Sum of squared errors.
    """
    model_temperatures = compute_melting_curve(Kprime_0, pressures, solid_params, liquid_params, T_m0, dH_fusion)
    interpolated_temps = np.interp(estimate_p, pressures, model_temperatures)
    return np.sum((interpolated_temps - estimate_temp) ** 2)

def simulated_annealing(initial_guess, temp_initial=1.0, temp_final=1e-6, cooling_rate=0.99, max_iter=100000):
    """
    Perform simulated annealing to find the best-fit K0' value golobally. 

    Parameters:
        initial_guess (float): 
            Initial guess for K0'.
        temp_initial (float): 
            Initial temperature for the annealing process.
        temp_final (float): 
            Final temperature for stopping criterion.
        cooling_rate (float): 
            Cooling rate for temperature reduction.
        max_iter (int): 
            Maximum number of iterations.

    Returns:
        Best K0' value and history of errors during the annealing process.
    """
    current_solution = initial_guess
    current_error = objective_function(current_solution)
    best_solution = current_solution
    best_error = current_error
    temperature = temp_initial

    error_history = []

    for iteration in range(max_iter):
        new_solution = current_solution + np.random.uniform(-1, 1)
        # Keep new_solution always bigger than 0
        new_solution = max(new_solution, 1e-3)

        new_error = objective_function(new_solution)

        if new_error < current_error or np.random.rand() < np.exp((current_error - new_error) / temperature):
            current_solution = new_solution
            current_error = new_error

            if new_error < best_error:
                best_solution = new_solution
                best_error = new_error

        error_history.append(best_error)
        temperature *= cooling_rate

        if temperature < temp_final:
            break

    return best_solution, error_history

def q3():
    """
    Use simulated annealing to find the best-fit K0' value and plot error convergence.
    """
    initial_guess = 10
    best_kprime_0, error_history = simulated_annealing(initial_guess)

    print(f"Best-fit K0' value: {best_kprime_0:.2f}")

    plt.figure(figsize=(10, 6))
    plt.plot(error_history, marker='.', color='blue')
    plt.xlabel("Iteration")
    plt.ylabel("Best Sum of Squared Errors")
    plt.title("Error vs Simulated Annealing Method Iteration")
    plt.grid()
    plt.show()
