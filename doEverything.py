import os

massBins = ['b:33:34']
tBins    = ['t:0']
modes    = ['-f', '-bw']
#uses     = ['u:0','u:1','u:0:1']
uses     = ['u:0:1:2:3','u:0:1:2','u:0:1:3', 'u:0:2:3' , 'u:1:2:3', 'u:0:1', 'u:0:2', 'u:0:3', 'u:1:2', 'u:1:3', 'u:2:3', 'u:0', 'u:1', 'u:2', 'u:3']

#for mb in massBins:
#	for tb in tBins:
#		for mode in modes:
#			for use in uses:
#				if mode == '-bw' and ':0' in use:
#					continue
#				os.system("python consistencyStudies.py "+mb+" "+tb+" "+mode+" "+use)

for mb in massBins:
	for tb in tBins:
		os.system("python consistencyStudies.py "+mb+" "+tb+" u:0 -P r:0:0.0:0.94")



