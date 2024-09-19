#!/usr/bin/env python3

'''
This file performs fire/disease spread simulations.

To get solution for lab 1: Run these commands:

>>> run lab01.py
>>> plt.ion()
>>> question_1_3x3()
>>> question_1_5x7()
>>> question_2_pspread_analysis()
>>> question_2_pbare_analysis()
>>> question_3_pfatal_high_analysis()
>>> question_3_pfatal_mid_analysis()
>>> question_3_pfatal_low_analysis()

'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')


def question_1_3x3(nNorth=3, nEast=3, maxiter=4, pspread=1):
    '''
    This function performs a fire spread simultion.

    Parameters
    ==========
    nNorth, nEast : integer, defaults to 3
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : int, defaults to 4
        Set the maximum number of iterations including initial condition
    pspread: float, value varies for different senarios
        For this function, fire will be 100% spread to the adjacent grids
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2

    # Set fire! To the center of the forest.
    forest[0, nNorth//2, nEast//2] = 3

    # Set pspread to a value
    pspread=1

    # Plot initial condition
    fig, ax = plt.subplots(1, 1)
    contour = ax.matshow(forest[0, :, :], vmin=1, vmax=3)
    ax.set_title(f'Iteration = {0:03d}')
    plt.colorbar(contour, ax=ax)

    # Propagate the solution.
    for k in range(maxiter-1):
       
        # Use current step to set next step:
        forest[k+1, :, :] = forest[k, :, :]

        # Burn in each cardinal direction.
        #From north to south:
        for i in range(nNorth - 1):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i+1, j] == 2) &\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i+1, j] = 3

        #From south to north:
        for i in range(1, nNorth):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i-1, j] == 2)&\
                    (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i-1, j] = 3

        #From east to west
        for i in range(nNorth):
            for j in range(nEast-1):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j+1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j+1] = 3

        #From west to east 
        for i in range(nNorth):
            for j in range(1,nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j-1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j-1] = 3


        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3  # Find cells that WERE burning
        forest[k+1, wasburn] = 1       # ...they are NOW bare.

        fig, ax = plt.subplots(1, 1)
        contour = ax.matshow(forest[k+1, :, :], vmin=1, vmax=3)
        ax.set_title(f'Iteration = {k+1:03d}')
        plt.colorbar(contour, ax=ax)

def question_1_5x7(nNorth=5, nEast=7, maxiter=7, pspread=1):
    '''
    This function performs a fire spread simultion.

    Parameters
    ==========
    nNorth, nEast : integer, defaults to 5
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : int, defaults to 7
        Set the maximum number of iterations including initial condition
    pspread: float, value varies for different senarios
        For this function, fire will be 100% spread to the adjacent grids
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2

    # Set fire! To the center of the forest.
    forest[0, nNorth//2, nEast//2] = 3

    # Set pspread to a value
    pspread=1

    # Plot initial condition
    fig, ax = plt.subplots(1, 1)
    contour = ax.matshow(forest[0, :, :], vmin=1, vmax=3)
    ax.set_title(f'Iteration = {0:03d}')
    plt.colorbar(contour, ax=ax)

    # Propagate the solution.
    for k in range(maxiter-1):
       
        # Use current step to set next step:
        forest[k+1, :, :] = forest[k, :, :]

        # Burn in each cardinal direction.
        #From north to south:
        for i in range(nNorth - 1):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i+1, j] == 2) &\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i+1, j] = 3

        #From south to north:
        for i in range(1, nNorth):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i-1, j] == 2)&\
                    (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i-1, j] = 3

        #From east to west
        for i in range(nNorth):
            for j in range(nEast-1):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j+1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j+1] = 3

        #From west to east 
        for i in range(nNorth):
            for j in range(1,nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j-1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j-1] = 3


        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3  # Find cells that WERE burning
        forest[k+1, wasburn] = 1       # ...they are NOW bare.

        fig, ax = plt.subplots(1, 1)
        contour = ax.matshow(forest[k+1, :, :], vmin=1, vmax=3)
        ax.set_title(f'Iteration = {k+1:03d}')
        plt.colorbar(contour, ax=ax)

def fire_simulation_pspread(nNorth, nEast, maxiter, pspread):
    '''
    This function performs a fire spread simultion. The simulation will show how 
    different value of the probablity of fire spreading (pspread) will impact the 
    wildfire evolution

    Parameters
    ==========
    nNorth, nEast : integer
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : integer 
        Set the maximum number of iterations including initial condition
    pspread: float, value varies for different senarios
        Set the probablity of fire spreading to the surrounding. 
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2

    # Set fire in the center of the forest.
    forest[0, nNorth // 2, nEast // 2] = 3

    # Propagate the solution.
    for k in range(maxiter-1):
       # Count the number of burning, bare, and forested grids
        burning_grid = np.sum(forest[k, :, :] == 3)
        bare_grid = np.sum(forest[k, :, :] == 1)
        forest_grid = np.sum(forest[k, :, :] == 2)

        # Use current step to set next step:
        forest[k+1, :, :] = forest[k, :, :]

        # Burn in each cardinal direction.
        #From north to south:
        for i in range(nNorth - 1):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i+1, j] == 2) &\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i+1, j] = 3

        #From south to north:
        for i in range(1, nNorth):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i-1, j] == 2)&\
                    (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i-1, j] = 3

        #From east to west
        for i in range(nNorth):
            for j in range(nEast-1):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j+1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j+1] = 3

        #From west to east 
        for i in range(nNorth):
            for j in range(1,nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j-1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j-1] = 3


        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3  # Find cells that WERE burning
        forest[k+1, wasburn] = 1       # ...they are NOW bare.
     
    return burning_grid, bare_grid, forest_grid

def question_2_pspread_analysis(nNorth=30, nEast=30, maxiter=100, num_runs=20):
    '''
    Analyze how varying the probability of fire spreading 
    (pspread) impacts wildfire evolution at the end of the max iteration

    Parameters
    ==========
    nNorth, nEast : integer, default to 30
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : integer, default to 100
        Set the maximum number of iterations including initial condition.
        A large number of interation can reduce the chance to have burning 
        grid at the end of the max iteration
    num_runs : integer, default to 20
        Set the times for running the code. The result will be the average of
        all the runs, which can reduce the potential outliers. 

    '''
    pspread_values = np.arange(0.1, 1.1, 0.1)
   
    # Initialize variables to store results
    avg_burning_counts = {}
    avg_bare_counts = {}
    avg_forest_counts = {}
    
    for pspread in pspread_values:
        total_burning_count = 0
        total_bare_count = 0
        total_forest_count = 0

        for i in range(num_runs):
            # Get the counts for each iteration
            burning_grid, bare_grid, forest_grid = \
                fire_simulation_pspread(nNorth, nEast, maxiter, pspread)
            total_burning_count += burning_grid
            total_bare_count += bare_grid
            total_forest_count += forest_grid
            i += 1

        # Calculate the average counts
        avg_burning_count = total_burning_count / num_runs
        avg_bare_count = total_bare_count / num_runs
        avg_forest_count = total_forest_count / num_runs
        
        # Save the results to variables for the later plots
        avg_burning_counts[pspread] = avg_burning_count
        avg_bare_counts[pspread] = avg_bare_count
        avg_forest_counts[pspread] = avg_forest_count

        # Print the numbers in the terminal
        print(f'pspread = {pspread:.1f}')
        print('Average Burning count:', avg_burning_count)
        print('Average Bare count:', avg_bare_count)
        print('Average Forest count:', avg_forest_count)
        
    
    # Plot the results
    plt.figure(figsize=(15, 5))
    plt.suptitle('Fire Spread Simulation pspread Analysis', fontsize=20)
    # Plot for Average Burning Count
    plt.subplot(1, 3, 1)
    # For the y-axis, the code iterates over each value in pbare_values 
    # and retrieves the corresponding average burning 
    # count from the avg_burning_counts
    plt.plot(pspread_values, [avg_burning_counts[pspread] \
                              for pspread in pspread_values], marker='o')
    plt.title('Average Burning Count')
    plt.xlabel('pspread')
    plt.ylabel('Number of Grids')

    # Plot for Average Bare Count
    plt.subplot(1, 3, 2)
    plt.plot(pspread_values, [avg_bare_counts[pspread] \
                              for pspread in pspread_values], marker='o')
    plt.title('Average Bare Count')
    plt.xlabel('pspread')
    plt.ylabel('Number of Grids')

    # Plot for Average Forest Count
    plt.subplot(1, 3, 3)
    plt.plot(pspread_values, [avg_forest_counts[pspread] \
                              for pspread in pspread_values], marker='o')
    plt.title('Average Forest Count')
    plt.xlabel('pspread')
    plt.ylabel('Number of Grids')

    plt.tight_layout()
    plt.show()
    
def fire_simulation_pbare(nNorth, nEast, maxiter, pspread, pbare):
    '''
    This function performs a fire spread simultion. The simulation will show how 
    different value of the probablity of bare land (pbare) will impact the 
    wildfire evolution

    Parameters
    ==========
    nNorth, nEast : integer
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : integer 
        Set the maximum number of iterations including initial condition
    pspread: float, default to 1
        Set the probablity of fire spreading to the surrounding. 
    pbare: float, varys for different senarios
        Set the probablity of bare gird for the initial condition.
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2
    
    # Set fire in the center of the forest.
    forest[0, nNorth // 2, nEast // 2] = 3

    # Create bare land and set initial condition
    for i in range(nNorth):
            for j in range(nEast):
                 if (np.random.rand()<=pbare) & (forest[0, i, j] != 3):
                    forest[0,i,j]= 1

    

    # Propagate the solution.
    for k in range(maxiter-1):
       # Count the number of burning, bare, and forested grids
        burning_grid = np.sum(forest[k, :, :] == 3)
        bare_grid = np.sum(forest[k, :, :] == 1)
        forest_grid = np.sum(forest[k, :, :] == 2)

        # Use current step to set next step:
        forest[k+1, :, :] = forest[k, :, :]

        # Burn in each cardinal direction.
        #From north to south:
        for i in range(nNorth - 1):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i+1, j] == 2) &\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i+1, j] = 3

        #From south to north:
        for i in range(1, nNorth):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i-1, j] == 2)&\
                    (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i-1, j] = 3

        #From east to west
        for i in range(nNorth):
            for j in range(nEast-1):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j+1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j+1] = 3

        #From west to east 
        for i in range(nNorth):
            for j in range(1,nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j-1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j-1] = 3


        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3  # Find cells that WERE burning
        forest[k+1, wasburn] = 1       # ...they are NOW bare.
     
    return burning_grid, bare_grid, forest_grid

def question_2_pbare_analysis(nNorth=30, nEast=30, maxiter=100, num_runs=20):
    '''
    Analyze how varying the probability of bare land 
    impacts wildfire evolution at the end of the max iteration.

    Parameters
    ==========
    nNorth, nEast : integer, default to 30
        Set the north-south (i) and east-west (j) size of grid.
    maxiter : integer, default to 100
        Set the maximum number of iterations including initial condition.
    num_runs : integer, default to 20
        Set the times for running the code. The result will be the average of
        all the runs, which can reduce the potential outliers.
    '''
    pbare_values = np.arange(0.1, 1.1, 0.1)

    # Initialize variables to store results
    avg_burning_counts = {}
    avg_bare_counts = {}
    avg_forest_counts = {}

    for pbare in pbare_values:
        total_burning_count = 0
        total_bare_count = 0
        total_forest_count = 0

        for i in range(num_runs):
            # Get the counts for each iteration
            burning_grid, bare_grid, forest_grid = \
                fire_simulation_pbare(nNorth, nEast, maxiter, pspread=1,\
                                       pbare=pbare)
            total_burning_count += burning_grid
            total_bare_count += bare_grid
            total_forest_count += forest_grid
            i += i

        # Calculate the average counts
        avg_burning_count = total_burning_count / num_runs
        avg_bare_count = total_bare_count / num_runs
        avg_forest_count = total_forest_count / num_runs

        # Save the results to variables for the later plots
        avg_burning_counts[pbare] = avg_burning_count
        avg_bare_counts[pbare] = avg_bare_count
        avg_forest_counts[pbare] = avg_forest_count

        # Print the numbers in the terminal
        print(f'pbare = {pbare:.1f}')
        print('Average Burning count:', avg_burning_count)
        print('Average Bare count:', avg_bare_count)
        print('Average Forest count:', avg_forest_count)
        

    # Plot the results
    plt.figure(figsize=(15, 5))
    plt.suptitle('Fire Spread Simulation pbare Analysis', fontsize=20)
    # Plot for Average Burning Count
    plt.subplot(1, 3, 1)
    # For the y-axis, the code iterates over each value in pspread_values 
    # and retrieves the corresponding average burning 
    # count from the avg_burning_counts
    plt.plot(pbare_values, [avg_burning_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Burning Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Bare Count
    plt.subplot(1, 3, 2)
    plt.plot(pbare_values, [avg_bare_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Bare Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Forest Count
    plt.subplot(1, 3, 3)
    plt.plot(pbare_values, [avg_forest_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Forest Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    plt.tight_layout()
    plt.show()

def disease_simulation_pfatal(nNorth, nEast, maxiter, pspread, pbare, pfatal):
    '''
    This function performs a fire spread simultion. The simulation will show how 
    different value of the probablity of fire spreading (pspread) will impact the 
    wildfire evolution

    Parameters
    ==========
    nNorth, nEast : integer
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : integer 
        Set the maximum number of iterations including initial condition
    pspread: float, default to 1
        Set the probablity of fire spreading to the surrounding. 
    pbare: float, varys for different senarios
        Set the probablity of bare gird for the initial condition.
    pfatal: float, varys for different senarios
        Set the probablity of death after being sick.
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2
    
    # Set fire in the center of the forest.
    forest[0, nNorth // 2, nEast // 2] = 3

    # Create bare land and set initial condition
    for i in range(nNorth):
            for j in range(nEast):
                 if (np.random.rand()<=pbare) & (forest[0, i, j] != 3):
                    forest[0,i,j]= 1

    

    # Propagate the solution.
    for k in range(maxiter-1):
       # Count the number of burning, bare, and forested grids
        burning_grid = np.sum(forest[k, :, :] == 3)
        bare_grid = np.sum(forest[k, :, :] == 1)
        forest_grid = np.sum(forest[k, :, :] == 2)
        fatal_grid = np.sum(forest[k, :, :] == 0)

        # Use current step to set next step:
        forest[k+1, :, :] = forest[k, :, :]

        # Burn in each cardinal direction.
        #From north to south:
        for i in range(nNorth - 1):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i+1, j] == 2) &\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i+1, j] = 3

        #From south to north:
        for i in range(1, nNorth):
            for j in range(nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i-1, j] == 2)&\
                    (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i-1, j] = 3

        #From east to west
        for i in range(nNorth):
            for j in range(nEast-1):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j+1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j+1] = 3

        #From west to east 
        for i in range(nNorth):
            for j in range(1,nEast):
                # Is current patch burning AND adjacent forested?
                if (forest[k, i, j] == 3) & (forest[k, i, j-1] == 2)&\
                      (np.random.rand()<=pspread):
                        # Spread fire to new square:
                        forest[k+1, i, j-1] = 3


        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3  # Find cells that WERE burning
        
        for i in range(nNorth):
            for j in range(nEast):
                if wasburn[i, j]:
                    if np.random.rand()<=pfatal:
                        forest[k+1, wasburn] = 0       # People are dead.
                    else:
                        forest[k+1, wasburn] = 1       # People are survived
     
    return burning_grid, bare_grid, forest_grid, fatal_grid

def question_3_pfatal_high_analysis(nNorth=30, nEast=30, maxiter=20, num_runs=20):
    '''
    Analyze how varying the probability of fire spreading 
    (pspread) impacts wildfire evolution at the end of the max iteration.
    pfatal=0.9


    Parameters
    ==========
    nNorth, nEast : integer, default to 30
        Set the north-south (i) and east-west (j) size of grid.
    maxiter : integer, default to 100
        Set the maximum number of iterations including initial condition.
    num_runs : integer, default to 20
        Set the times for running the code. The result will be the average of
        all the runs, which can reduce the potential outliers.
    '''
    pbare_values = np.arange(0.1, 1.1, 0.1)

    # Initialize variables to store results
    avg_burning_counts = {}
    avg_bare_counts = {}
    avg_forest_counts = {}
    avg_fatal_counts = {}

    for pbare in pbare_values:
        total_burning_count = 0
        total_bare_count = 0
        total_forest_count = 0
        total_fatal_count = 0

        for i in range(num_runs):
            # Get the counts for each iteration
            burning_grid, bare_grid, forest_grid, fatal_grid = \
                disease_simulation_pfatal(nNorth, nEast, maxiter, pspread=1,\
                                       pbare=pbare, pfatal=0.9)
            total_burning_count += burning_grid
            total_bare_count += bare_grid
            total_forest_count += forest_grid
            total_fatal_count += fatal_grid
            i += i

        # Calculate the average counts
        avg_burning_count = total_burning_count / num_runs
        avg_bare_count = total_bare_count / num_runs
        avg_forest_count = total_forest_count / num_runs
        avg_fatal_count = total_fatal_count / num_runs

        # Save the results to variables for the later plots
        avg_burning_counts[pbare] = avg_burning_count
        avg_bare_counts[pbare] = avg_bare_count
        avg_forest_counts[pbare] = avg_forest_count
        avg_fatal_counts[pbare] = avg_fatal_count

        # Print the numbers in the terminal
        print(f'pbare = {pbare:.1f}')
        print('Average Burning count:', avg_burning_count)
        print('Average Bare count:', avg_bare_count)
        print('Average Forest count:', avg_forest_count)
        print('Average Fatal count:', avg_fatal_count)
        print('Average Mortality Rate:', avg_fatal_count/(nNorth*nEast))
        

    # Plotting the results
    plt.figure(figsize=(25, 5))
    plt.suptitle('Disease Spread Simulation pbare Analysis, pfatal=0.9', fontsize=20)
    
    # Plot for Average Burning Count
    plt.subplot(1, 5, 1)
    
    # For the y-axis, the code iterates over each value in pspread_values 
    # and retrieves the corresponding average burning 
    # count from the avg_burning_counts
    plt.plot(pbare_values, [avg_burning_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Sick Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Bare Count
    plt.subplot(1, 5, 2)
    plt.plot(pbare_values, [avg_bare_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Immune Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Forest Count
    plt.subplot(1, 5, 3)
    plt.plot(pbare_values, [avg_forest_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Healthy Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Fatal Count
    plt.subplot(1, 5, 4)
    plt.plot(pbare_values, [avg_fatal_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Fatal Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

     # Plot for Average Mortality Rate
    plt.subplot(1, 5, 5)
    plt.plot(pbare_values, [avg_fatal_counts[pbare]/(nNorth*nEast) \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Mortality Rate')
    plt.xlabel('pbare')
    plt.ylabel('Ratio')

    plt.tight_layout()
    plt.show()

def question_3_pfatal_mid_analysis(nNorth=30, nEast=30, maxiter=20, num_runs=20):
    '''
    Analyze how varying the probability of disease spreading 
    (pspread) impacts wildfire evolution at the end of the max iteration.
    pfatal=0.5

    Parameters
    ==========
    nNorth, nEast : integer, default to 30
        Set the north-south (i) and east-west (j) size of grid.
    maxiter : integer, default to 100
        Set the maximum number of iterations including initial condition.
    num_runs : integer, default to 20
        Set the times for running the code. The result will be the average of
        all the runs, which can reduce the potential outliers.
    '''
    pbare_values = np.arange(0.1, 1.1, 0.1)

    # Initialize variables to store results
    avg_burning_counts = {}
    avg_bare_counts = {}
    avg_forest_counts = {}
    avg_fatal_counts = {}

    for pbare in pbare_values:
        total_burning_count = 0
        total_bare_count = 0
        total_forest_count = 0
        total_fatal_count = 0

        for i in range(num_runs):
            # Get the counts for each iteration
            burning_grid, bare_grid, forest_grid, fatal_grid = \
                disease_simulation_pfatal(nNorth, nEast, maxiter, pspread=1,\
                                       pbare=pbare, pfatal=0.5)
            total_burning_count += burning_grid
            total_bare_count += bare_grid
            total_forest_count += forest_grid
            total_fatal_count += fatal_grid
            i += i

        # Calculate the average counts
        avg_burning_count = total_burning_count / num_runs
        avg_bare_count = total_bare_count / num_runs
        avg_forest_count = total_forest_count / num_runs
        avg_fatal_count = total_fatal_count / num_runs

        # Save the results to variables for the later plots
        avg_burning_counts[pbare] = avg_burning_count
        avg_bare_counts[pbare] = avg_bare_count
        avg_forest_counts[pbare] = avg_forest_count
        avg_fatal_counts[pbare] = avg_fatal_count

        # Print the numbers in the terminal
        print(f'pbare = {pbare:.1f}')
        print('Average Burning count:', avg_burning_count)
        print('Average Bare count:', avg_bare_count)
        print('Average Forest count:', avg_forest_count)
        print('Average Fatal count:', avg_fatal_count)
        print('Average Mortality Rate:', avg_fatal_count/(nNorth*nEast))


    # Plotting the results
    plt.figure(figsize=(25, 5))
    plt.suptitle('Disease Spread Simulation pbare Analysis, pfatal=0.5', fontsize=20)
    
    # Plot for Average Burning Count
    plt.subplot(1, 5, 1)
    
    # For the y-axis, the code iterates over each value in pspread_values 
    # and retrieves the corresponding average burning 
    # count from the avg_burning_counts
    plt.plot(pbare_values, [avg_burning_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Sick Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Bare Count
    plt.subplot(1, 5, 2)
    plt.plot(pbare_values, [avg_bare_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Immune Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Forest Count
    plt.subplot(1, 5, 3)
    plt.plot(pbare_values, [avg_forest_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Healthy Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Fatal Count
    plt.subplot(1, 5, 4)
    plt.plot(pbare_values, [avg_fatal_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Fatal Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

     # Plot for Average Mortality Rate
    plt.subplot(1, 5, 5)
    plt.plot(pbare_values, [avg_fatal_counts[pbare]/(nNorth*nEast) \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Mortality Rate')
    plt.xlabel('pbare')
    plt.ylabel('Ratio')

    plt.tight_layout()
    plt.show()

def question_3_pfatal_low_analysis(nNorth=30, nEast=30, maxiter=20, num_runs=20):
    '''
    Analyze how varying the probability of disease spreading 
    (pspread) impacts wildfire evolution at the end of the max iteration.
    pfatal=0.1

    Parameters
    ==========
    nNorth, nEast : integer, default to 30
        Set the north-south (i) and east-west (j) size of grid.
    maxiter : integer, default to 100
        Set the maximum number of iterations including initial condition.
    num_runs : integer, default to 20
        Set the times for running the code. The result will be the average of
        all the runs, which can reduce the potential outliers.
    '''
    #create a varible for the later loop
    pbare_values = np.arange(0.1, 1.1, 0.1)

    # Initialize variables to store results
    avg_burning_counts = {}
    avg_bare_counts = {}
    avg_forest_counts = {}
    avg_fatal_counts = {}

    for pbare in pbare_values:
        total_burning_count = 0
        total_bare_count = 0
        total_forest_count = 0
        total_fatal_count = 0

        for i in range(num_runs):
            # Get the counts for each iteration
            burning_grid, bare_grid, forest_grid, fatal_grid = \
                disease_simulation_pfatal(nNorth, nEast, maxiter, pspread=1,\
                                       pbare=pbare, pfatal=0.1)
            total_burning_count += burning_grid
            total_bare_count += bare_grid
            total_forest_count += forest_grid
            total_fatal_count += fatal_grid
            i += i

        # Calculate the average counts
        avg_burning_count = total_burning_count / num_runs
        avg_bare_count = total_bare_count / num_runs
        avg_forest_count = total_forest_count / num_runs
        avg_fatal_count = total_fatal_count / num_runs

        # save the results to variables for the later plots
        avg_burning_counts[pbare] = avg_burning_count
        avg_bare_counts[pbare] = avg_bare_count
        avg_forest_counts[pbare] = avg_forest_count
        avg_fatal_counts[pbare] = avg_fatal_count

        # Print the numbers in the terminal
        print(f'pbare = {pbare:.1f}')
        print('Average Burning count:', avg_burning_count)
        print('Average Bare count:', avg_bare_count)
        print('Average Forest count:', avg_forest_count)
        print('Average Fatal count:', avg_fatal_count)
        print('Average Mortality Rate:', avg_fatal_count/(nNorth*nEast))
        print('---')

    # Plotting the results
    plt.figure(figsize=(25, 5))
    plt.suptitle('Disease Spread Simulation pbare Analysis pfatal=0.1', fontsize=20)
    
    # Plot for Average Burning Count
    plt.subplot(1, 5, 1)
    
    # For the y-axis, the code iterates over each value in pspread_values 
    # and retrieves the corresponding average burning 
    # count from the avg_burning_counts
    plt.plot(pbare_values, [avg_burning_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Sick Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Bare Count
    plt.subplot(1, 5, 2)
    plt.plot(pbare_values, [avg_bare_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Immune Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Forest Count
    plt.subplot(1, 5, 3)
    plt.plot(pbare_values, [avg_forest_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Healthy Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

    # Plot for Average Fatal Count
    plt.subplot(1, 5, 4)
    plt.plot(pbare_values, [avg_fatal_counts[pbare] \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Fatal Count')
    plt.xlabel('pbare')
    plt.ylabel('Number of Grids')

     # Plot for Average Mortality Rate
    plt.subplot(1, 5, 5)
    plt.plot(pbare_values, [avg_fatal_counts[pbare]/(nNorth*nEast) \
                            for pbare in pbare_values], marker='o')
    plt.title('Average Mortality Rate')
    plt.xlabel('pbare')
    plt.ylabel('Ratio')

    plt.tight_layout()
    plt.show()
