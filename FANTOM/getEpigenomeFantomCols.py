#!/nfs/users/nfs_j/js29/anaconda/bin/python
##
import argparse
import sys
import subprocess
import re
import os
import os.path

parser = argparse.ArgumentParser(description="Writes a table mapping Roadmap epigenomes to FANTOM columns and weightings.")
parser.add_argument("--roadmap", type=file, required=True, metavar='FILE', help="Roadmap metadata table")
parser.add_argument("--fantom", type=file, required=True, metavar='FILE', help="Fantom metadata table")
parser.add_argument("--enhancer", type=file, metavar='FILE', help="Output table for enhancers rather than TSS (TSS is default)")
parser.add_argument("--output", type=str, default="out", metavar='STR', help="Path and prefix for output files")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")
args = parser.parse_args()

# Stores the full line from hg19.cage_peak_tpm_ann.osc.colweights.txt keyed by FantomID
fantomTssDict = {}
# Stores just the column ID from human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz.cols keyed by FantomID
fantomEnhDict = {}

def main():
	global fantomTssDict
	global fantomEnhDict
	
	# Read entire FANTOM TSS metadata table into a dict keyed by FantomID
	for line in args.fantom:
		lineVals = line.strip().split("\t")
		fantomTssDict[lineVals[3]] = lineVals
	
	if args.enhancer:
		for line in args.enhancer:
			lineVals = line.strip().split("\t")
			fantomEnhDict[lineVals[0]] = lineVals[1]
	
	tssColsFname = args.output + ".tsscols.txt"
	tssColsFile = open(tssColsFname, 'w')
	
	if args.enhancer:
		enhColsFname = args.output + ".enhcols.txt"
		enhColsFile = open(enhColsFname, 'w')
	
	for line in args.roadmap:
		lineVals = line.strip("\n").split("\t")
		if lineVals[0] == "Epigenome ID":
			continue
		if lineVals[5].lower() == "exclude":
			continue
			
		# Col 5 has a comma-separated list of Fantom IDs to use to define TSSes,
		# e.g. "CNhs12824,CNhs12837". Col 6 has the same for enhancers.
		colsStr = None
		weightsStr = None
		fantomTssIDs = lineVals[5].split(",")
		for fantomID in fantomTssIDs:
			if not fantomID in fantomTssDict:
				die("Fantom ID '{0}' not found in Fantom TSS table".format(fantomID))
			tssLine = fantomTssDict[fantomID]
			if colsStr is None:
				colsStr = tssLine[0]
				weightsStr = tssLine[2]
			else:
				colsStr += "," + tssLine[0]
				weightsStr += "," + tssLine[2]

		tssColsFile.write("\t".join([lineVals[0], colsStr, weightsStr]) + "\n")
				
		if args.enhancer:
			colsStr = None
			weightsStr = None
			if len(lineVals[6]) == 0:
				if args.debug:
					sys.stderr.write("No fantom enhancers for epigenome {0}\n".format(lineVals[0]))
				continue
			fantomEnhIDs = lineVals[6].split(",")
			for fantomID in fantomEnhIDs:
				if not fantomID in fantomEnhDict:
					die("Fantom ID '{0}' not found in Fantom Enhancers table".format(fantomID))
				tssLine = fantomTssDict[fantomID]
				if colsStr is None:
					colsStr = fantomEnhDict[fantomID]
					weightsStr = tssLine[2]
				else:
					colsStr += "," + fantomEnhDict[fantomID]
					weightsStr += "," + tssLine[2]
			enhColsFile.write("\t".join([lineVals[0], colsStr, weightsStr]) + "\n")
		

def die(msg):
	sys.stderr.write(msg + "\n")
	exit(1)

def docall(cmd):
	sys.stderr.write(cmd + "\n")
	retval = subprocess.call(cmd, shell=True)
	if (retval != 0):
		die("Error running command.")

def docheckoutput(cmd):
	sys.stderr.write(cmd + "\n")
	retval = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	return retval


main()
