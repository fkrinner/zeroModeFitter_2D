import os

massBins = ['b:31:32']
tBins    = ['t:0']
modes    = ['-f', '-bw']
uses     = ['u:0','u:1','u:0:1']

for mb in massBins:
	for tb in tBins:
		for mode in modes:
			for use in uses:
				if mode == '-bw' and not use == 'u:1':
					continue
				os.system("python consistencyStudies.py "+mb+" "+tb+" "+mode+" "+use)



