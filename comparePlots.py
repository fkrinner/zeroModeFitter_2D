import os, sys


def main():
	fileName = sys.argv[1]
	if "comparisonResultsData" in fileName:
		replacer  = "comparisonResultsData"
		folderList =["comparisonResultsData_1mp","comparisonResultsData_3pp","comparisonResultsData_4mp","comparisonResultsData_4pp","comparisonResultsData_6mp","comparisonResultsData_bigger1pp","comparisonResultsData_bigger2pp","comparisonResultsData_bigger2mp"]
	else:
		replacer  = "comparisonResults_rest"
		folderList = ["comparisonResults_rest_1mp","comparisonResults_rest_3pp","comparisonResults_rest_4mp","comparisonResults_rest_4pp","comparisonResults_rest_6mp","comparisonResults_rest_bigger1pp","comparisonResults_rest_bigger2pp","comparisonResults_rest_bigger2mp"]

	stdFolder = "./"+replacer+"/"
	if not os.path.isfile(fileName):
		raise IOError(fileName+ " does not exist")
	for studyFolder in folderList:
		studyFileName = fileName.replace(replacer,studyFolder)
		if not os.path.isfile(studyFileName):
			print "Skipping", studyFileName
			continue
		else:
			print "Comparing", studyFileName
		os.system("evince "+studyFileName + " " + fileName)
		

if __name__ == "__main__":
	main()
