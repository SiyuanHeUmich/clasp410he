#!/usr/bin/env python3
'''
This file contains tools and scripts for completing Lab 1 for CLaSP 410.
To reproduce the plots shown in the lab report, do this...
'''
import numpy as np
import matplotlib.pyplot as plt 

plt.style.use('fivethirtyeight')
n_north, n_east = 3, 3 # Number of cells in X and Y direction.
prob_spread = 1.0 # Chance to spread to adjacent cells.
prob_bare = 0.0 # Chance of cell to start as bare patch.
prob_start = 0.0 # Chance of cell to start on fire.

# # Create an initial grid, set all values to "2". dtype sets the value
# # type in our array to integers only.
# forest = np.zeros([max iter, n_north, n_east], dtype=int) + 2
# # Set the center cell to "burning":
# forest[1, 1] = 3


# for k in range(max iter):
# #burn north to south
# # Loop over every cell in the x and y directions.
#     for i in range(n_north-1):
#         for j in range(n_east):
#             if forest[i,j]==3:
#                 if forest(i+1,j)==2:
#                     forest[i,j,k+1]=1
#                     forest[i,j,k+1]=3
#         # Roll our "dice" to see if we get a bare spot:
#             if np.random.rand() < prob_bare:
#                 forest[i, j] = 1 # 1 is a bare spot.

def fire_spread (n_north=3, n_east=3, maxiter=4):
    '''
    This function performs a fire/disease spread simulation.

    parameters
    =========
    n_north. n_east : integer, default to 3
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxier : int, default to 4
        Ser the maximum number of iterations

    '''
    # Create forest and set initial condition
    forest = np.zeros([maxiter, n_north, n_east]) +2

    # Set fire! to the center of the forest.
    forest[0,1,1]=3

    # Plot initial condition

    #propagate the solution
    for k in range(maxiter):
        
        # USe current step to set next step:
        forest[k+1, :, :]=forest[k, :, :]
        
        # burn in each cardinal direction.
        #from north to south:
        for i in range (n_north -1):
            for j in range (n_east -1):
                #is current patch burning and adjacent forested?
                if forest [k, i, j] ==3 & forest[k, i+1, j]==2:
                        #spread fire to new square
                        forest[k+1, i+1, j] = 3
        
        
        # Set currently burning to bare:
        wasburn = forest[k, :, :] == 3 
        #Find cells that were burning, now they are bare
        forest[k+1, wasburn]=1         

        fig, ax = plt.subplots(1, 1)
        contour = ax.matshow(forest[k+1,:,:], vmin=1, vmax=3)
        ax.set_title(f'Iteration = {k+1:03d}')
        plt.colorbar(contour, ax=ax)
fire_spread()
plt.show()