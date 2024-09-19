#!/usr/bin/env python3

'''
This file performs fire/disease spread simulations.

To get solution for lab 1: Run these commands:

>>> blah
>>> blah blah

'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')


def fire_spread(nNorth=3, nEast=5, maxiter=7, pspread=1):
    '''
    This function performs a fire/disease spread simultion.

    Parameters
    ==========
    nNorth, nEast : integer, defaults to 3
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : int, defaults to 4
        Set the maximum number of iterations including initial condition
    pspread
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

        # fig.savefig(f'fig{k:04d}.png')
        # plt.close('all')

# fig=fire_spread()
# plt.show()

def fire_spread_2(nNorth=3, nEast=5, maxiter=5, pspread=0.5, pbare=0.3):
    '''
    This function performs a fire/disease spread simultion.

    Parameters
    ==========
    nNorth, nEast : integer, defaults to 3
        Set the north-south (i) and east-west (j) size of grid.
        Default is 3 squares in each direction.
    maxiter : int, defaults to 4
        Set the maximum number of iterations including initial condition
    pspread
    '''

    # Create forest and set initial condition
    forest = np.zeros([maxiter, nNorth, nEast]) + 2

    # Create bare land and set initial condition
    for i in range(nNorth):
            for j in range(nEast):
                 if (np.random.rand()<=pbare):
                    forest = np.zeros([maxiter, nNorth, nEast]) + 1

    # Set fire! To the center of the forest.
    if (forest[k, i, j] == 3):
        forest[0, nNorth//2, nEast//2] = 3

    # Set pspread to a value
    #pspread=.8

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

        # fig.savefig(f'fig{k:04d}.png')
        # plt.close('all')