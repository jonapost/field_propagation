import numpy as np
import matplotlib.pyplot as plt

N = 3

ind = np.arange(N)  # the x locations for the groups
width = 0.35        # the width of the bars

fig, ax = plt.subplots()

G4MagInt_Driver = (0.0266, 0.0803, 0.36)
G4MagInt_Driver_errors = (0.000527, 0.00149, 0.00544)
#rects1 = ax.bar(ind, G4MagInt_Driver, width, color='g', yerr=G4MagInt_Driver_errors)

G4IntegrationDriver = (0.0237, 0.0779, 0.345)
G4IntegrationDriver_errors = (0.000424, 0.00137, 0.00515)

def calc_error(val1, error1, val2, error2):
    return val1/val2 * (error1/val1 + error2/val2)

data = []
error = []
for i in range(3):
    data.append((G4MagInt_Driver[i] - G4IntegrationDriver[i])/ G4MagInt_Driver[i] * 100)
    print (data[i])
    error.append(
        calc_error(
            G4IntegrationDriver[i],
            G4IntegrationDriver_errors[i],
            G4MagInt_Driver[i],
            G4MagInt_Driver_errors[i]))

rects2 = ax.bar(ind + width, data, width, color='b', yerr=error)

# add some text for labels, title and axes ticks
ax.set_ylabel('event user time decrease %')
ax.set_title('test NTST')
ax.set_xticks(ind + width) #+ width / 2)
ax.set_xticklabels(('run2a.mac', 'run2b.mac', 'run2c.mac'))

#ax.set_yscale('log')

#ax.legend((rects1[0], rects2[0]), ('G4MagInt_Driver', 'G4IntegrationDriver'))


plt.show()
