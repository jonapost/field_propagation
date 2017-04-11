#!/usr/local/bin/python3
import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def plot(fname):
    px = []
    py = []
    pz = []
    with open(fname) as f:
        data = f.read()
        row = data.split('\n')
        for irow in row:
            if len(irow) > 0:
                px.append(float(irow.split()[0]))
                py.append(float(irow.split()[1]))
                pz.append(float(irow.split()[2]))
	
        mpl.rcParams['legend.fontsize'] = 10
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        plt.plot(px, py, pz,label=fname)
        plt.legend()

files = input("file name: ")
for fname in files.split():
    plot(fname)

plt.show()
