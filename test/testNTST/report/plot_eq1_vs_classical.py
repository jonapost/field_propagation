import numpy as np
import matplotlib.pyplot as plt

N = 3

ind = np.arange(N)  # the x locations for the groups
width = 0.35        # the width of the bars

fig, ax = plt.subplots()

eq1 = (217686114, 780247592, 2405482710)
classical = (259251794, 1216830112, 4298927706)

data = []
for i in range(3):
    data.append((classical[i] - eq1[i]) / classical[i] * 100)
    print(data[i])

rects1 = ax.bar(ind, data, width, color='b')


# add some text for labels, title and axes ticks
ax.set_ylabel('field evaluations decrease %')
ax.set_title('test NTST')
ax.set_xticks(ind)
ax.set_xticklabels(('run2a.mac', 'run2b.mac', 'run2c.mac'))

#ax.set_yscale('log')

plt.show()
