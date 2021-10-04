#!/usr/bin/env python3

import sys # argv
print(len(sys.argv), 'arguments:', str(sys.argv))

import json # json.load
import matplotlib.pyplot as plt
import numpy as np

chemSymbol = 'Cu'
if len(sys.argv) > 1:
    chemSymbol = sys.argv[1] 

with open(chemSymbol+'.json') as jsonfile:
    data = json.load(jsonfile)
    # print(data)
    r   = np.array(data['radial grid']['values'])
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    ax.set_xlabel('Radius ('+str(data['length unit'])+')')  # Add an x-label to the axes.
    # ax.set_ylabel('y label')  # Add a y-label to the axes.
    ax.set_title(chemSymbol+', Z = '+str(data['atomic number'])) # Add a title to the axes.
    try:
        partial_waves = data['partial waves']
    except KeyError:
        # no key 'partial waves' in JSON file
        print('load data exported by a43/src/atom_core.cxx')
        rho = np.array(data['density'])
        rV  = np.array(data['r*potential'])
        states = list(data['eigenstates'].keys())
        print(chemSymbol, 'eigenstates', states)

        ax.plot(r, -rV, label='-r*potential')
        ax.plot(r, r**2*rho, label='r^2*rho')

        occ  = np.zeros((len(states)))
        wave = np.zeros((len(states),r.size))
        for i,state in enumerate(states):
            occ[i] = data['eigenstates'][state]['occupation']
            wave[i] = np.array(data['eigenstates'][state]['wave'])
            ax.plot(r, occ[i]*(r*wave[i])**2, label='rho_'+str(state))
    else:
        # key 'partial waves' found in JSON file
        print('load data exported by a43/src/single_atom.cxx')
        rV_smt = np.array(data['r*true potential'])
        rV_tru = np.array(data['r*smooth potential'])
        ax.plot(r, rV_tru, label='r*true potential',   color='black')
        ax.plot(r, rV_smt, label='r*smooth potential', color='black', linestyle='dashed')
        states = list(partial_waves.keys())
        print(chemSymbol, 'partial waves', states)
        for i,state in enumerate(states):
            psi_tru = np.array(partial_waves[state]['true wave'])
            psi_smt = np.array(partial_waves[state]['smooth wave'])
            color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(r, psi_tru, label=str(state), color=color)
            ax.plot(r, psi_smt, label=str(state), color=color, linestyle='dashed')

    ax.legend() # Add a legend.
    plt.savefig(chemSymbol+'.svg') # Save in vector format
    plt.savefig(chemSymbol+'.png') # Save in pixel format
