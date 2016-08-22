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
        vecList = data.split('\n')
        for row in vecList:
            if len(row) > 0:
					#print "row: ",row
					begin = end = 0
					end = row.find(',')
					val = row[begin+1:end]
					#print "begin: ",begin, " end: ", end, "val: ",val
					px.append(float(row[begin+1:end]))
					begin = end
					end = row.find(',',end+1)
					val = row[begin+1:end]
					#print "begin: ",begin, " end: ", end, "val: ",val
					py.append(float(row[begin+1:end]))
					begin = end
					end = row.find(')')
					val = row[begin+1:end]
					#print val
					pz.append(float(row[begin+1:end]))
		
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.plot(px, py, pz,label=fname)
    plt.legend()


plot("pos_log.txt")
plt.show()
