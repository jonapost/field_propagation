import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


px = []
py = []
pz = []
with open("out.txt") as f:
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
plt.plot(px, py, pz)

plt.show()
