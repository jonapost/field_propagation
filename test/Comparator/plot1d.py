#!/usr/local/bin/python3
import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def plot(fname, loglog):
    array = []
    with open(fname) as f:
        data = f.read()
        row = data.split('\n')
        for irow in row:
            if len(irow) > 0:
                array.append(float(irow))

        mpl.rcParams['legend.fontsize'] = 10
        if loglog:
            plt.loglog(array, label=fname)
        else:
            plt.plot(array, label=fname)
        plt.legend()

plotType = input("plot type: (loglog/standard)")
files = input("file name: ")
for fname in files.split():
    plot(fname, plotType == "loglog")

plt.show()
