import numpy as np
import matplotlib.pyplot as plt

N = 3

ind = np.arange(N)  # the x locations for the groups
width = 0.35        # the width of the bars

fig, ax = plt.subplots()

nonfsal = (217686114, 780247592, 2405482710)
fsal =    (217178591, 767103732, 2342675019)

data = []
error = []
for i in range(3):
    data.append((nonfsal[i] - fsal[i])/ nonfsal[i] * 100)
    print (data[i])

rects2 = ax.bar(ind + width, data, width, color='b')

# add some text for labels, title and axes ticks
ax.set_ylabel('field evaluations decrease %')
ax.set_title('test NTST')
ax.set_xticks(ind + width) #+ width / 2)
ax.set_xticklabels(('run2a.mac', 'run2b.mac', 'run2c.mac'))

#ax.set_yscale('log')

#ax.legend((rects1[0], rects2[0]), ('G4MagInt_Driver', 'G4IntegrationDriver'))


plt.show()
