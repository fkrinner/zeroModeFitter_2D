import numpy as np

def writeComa(vals, fileName, normed = False):
	coma = np.cov(vals)
	dim  = len(coma)
	with open(fileName, 'w') as outFile:
		for i in range(dim):
			for j in range(dim):
				val = coma[i,j]
				if normed:
					if not val == 0.:
						val/=(coma[i,i]*coma[j,j])**.5
				outFile.write(str(val) + ' ' )
			outFile.write('\n')
