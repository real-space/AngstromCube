#!/usr/bin/env python3

import sys # argv
print(len(sys.argv), 'arguments:', str(sys.argv))

import json # json.load
import matplotlib.pyplot as plt
import numpy as np

chemSymbol = 'Cu'
if len(sys.argv) > 1:
    chemSymbol = sys.argv[1] 

with open(chemSymbol+'.json') as f:
    data = json.load(f)
    # print(data)
    r   = np.array(data['radial grid']['values'])
    rho = np.array(data['density'])
    rV  = np.array(data['r*potential'])
    states = list(data['eigenstates'].keys())
    print(chemSymbol, states)
    occ  = np.zeros((len(states)))
    wave = np.zeros((len(states),r.size))
    for i,state in enumerate(states):
        occ[i] = data['eigenstates'][state]['occupation']
        wave[i] = np.array(data['eigenstates'][state]['wave'])

    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    ax.plot(r, -rV, label='-r*potential')
    ax.plot(r, r**2*rho, label='r^2*rho')
    for i,state in enumerate(states):
        ax.plot(r, occ[i]*(r*wave[i])**2, label='rho_'+state)
    ax.set_xlabel('Radius (Bohr)')  # Add an x-label to the axes.
    # ax.set_ylabel('y label')  # Add a y-label to the axes.
    ax.set_title(chemSymbol) # Add a title to the axes.
    ax.legend() # Add a legend.
    plt.savefig(chemSymbol+'.svg') # Save in vector format
    plt.savefig(chemSymbol+'.png') # Save in pixel format
