# pnGize.py
import os, sys

def main():
	folder = sys.argv[1]
	while folder.endswith(os.sep):
		folder = folder[:-1]
	fns    = []
	for fn in os.listdir(folder):
		if fn.endswith(".pdf"):
			fns.append(fn)
	newFolderName = folder+"_pngs"
	os.system("mkdir "+newFolderName)
	for fn in fns:
		print "Converting",fn
		os.system("convert " + folder + os.sep + fn + ' ' + newFolderName + os.sep + fn.replace(".pdf",".png"))

if __name__ == "__main__":
	main()
