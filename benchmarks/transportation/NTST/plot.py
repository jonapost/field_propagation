import matplotlib.pyplot as plt
x = [2**i for i in xrange(4,14)]
y = [i**2 for i in x]
plt.loglog(x,y,'ro',basex=2,basey=2)
plt.xlim([0, 2**14]) # <--- this line does nothing
plt.show()
