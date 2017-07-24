import numpy as np
import matplotlib.pyplot as plt

N = 3

ind = np.arange(N)  # the x locations for the groups
width = 0.25        # the width of the bars

fig, ax = plt.subplots()

eq1 = (217686114, 780247592, 2405482710)
eq1_errors = (0.000527, 0.00149, 0.00544)
rects1 = ax.bar(ind, eq1, width, color='r', yerr=eq1_errors)

eq2 = (222704371, 779841565, 2525802727)
eq2_errors = (0.000424, 0.00137, 0.00515)
rects2 = ax.bar(ind + width, eq2, width, color='g', yerr=eq2_errors)

eq3 = (217951569, 815614582, 2664384060)
eq3_errors = (0.000424, 0.00137, 0.00515)
rects3 = ax.bar(ind + 2 * width, eq3, width, color='b', yerr=eq3_errors)

# add some text for labels, title and axes ticks
ax.set_ylabel('number of field evaluations')
ax.set_title('test NTST')
ax.set_xticks(ind + width)
ax.set_xticklabels(('run2a.mac', 'run2b.mac', 'run2c.mac'))

#ax.set_yscale('log')

ax.legend((rects1[0], rects2[0], rects3[0]), ('RK547FEq1', 'RK547FEq2', 'RK547FEq3'))


plt.show()
