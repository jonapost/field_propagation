from matplotlib import pyplot as plt

files = ["ClassicalRK.txt",
		 "CashKarp.txt",
		 "BogackiShampine45.txt",
		 "BulirschStoer.txt",
		 "BogackiShampine45Dense.txt",
		 "BulirschStoerDense.txt",
		 "BoostBulirschStoer.txt"
		]

for ifile in files:
	#print "read",ifile
	arcLength = []
	pos_error = []
	mom_error = []
	field_calls = []
	total_calls = 0
	with open(ifile) as f:
		data = f.read()
		line = data.split('\n')
		for iline in line:
			row = iline.split()
			if len(row) == 4:
				arcLength.append(float(row[0]))
				pos_error.append(float(row[1]))
				mom_error.append(float(row[2]))
				field_calls.append(float(row[3]))
				total_calls = total_calls + float(row[3])
	print ifile, "total calls: ", total_calls	
	plt.plot(arcLength, field_calls,linestyle='-',marker='o', label = ifile, linewidth = 5, markersize = 10)
plt.legend(loc='lower right',fontsize= 25)
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('arc length', fontsize = 30)
plt.ylabel('momentum error', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize = 30)
plt.show()	
