#!/usr/bin/python
import argparse
import sys
import os
import os.path
import re
import string

def die(msg):
	sys.stderr.write(msg)
	exit(1)

parser = argparse.ArgumentParser(description="Gets the peaks for each gene with a minimum TPM threshold.")

parser.add_argument("infile", type=file, nargs='?', metavar='FILE',
					help="input file name (or use stdin)")
parser.add_argument("--mintpm", type=float, required = True, metavar='NUM',
					help="do not include any peak with less than this TPM")
parser.add_argument("--genes", action='store_true', help="only include TSS from genes")
args = parser.parse_args()

if (args.mintpm < 0): exit(0)

infile = args.infile
if (infile == None):
	infile = sys.stdin

# Store all gene details in a hash until we have read the whole file and know that
# we have all TSS for each gene.
geneTssDict = {}

lineNum = 0
tssArray = []
header = infile.readline()
for line in infile:
	lineNum += 1
	if (re.search("\s*#", line) is not None):
		continue
	line = line.strip()
	#sys.stderr.write(line + "\n")

	fields = line.split()
	detailFields = fields[1].split(',')
	if (len(detailFields) < 1):
		die("Error: unable to parse peak/score field at line %(lineNum).\n")

	tpm = float(fields[-1])
	# Ignore all peaks with zero TPM in this cell line
	if (tpm <= 0 and args.mintpm > 0):
		continue

	peakFields = detailFields[0].split('@')
	if (len(peakFields) < 2):
		die("Error: unable to parse peak/gene field at line %(lineNum).\n")
	if (peakFields[0][:1] != 'p'):
		die("Error: expected peak name to begin with 'p': line %(lineNum).\n")
	# Ignore peak fields that don't correspond to a particular gene (i.e. have no peak number)
	peakNum = peakFields[0][1:]
	if (args.genes and (peakNum == "" or int(peakNum) <= 0)):
		continue

	geneName = peakFields[1]
	if geneName in geneTssDict:
		geneTssDict[geneName].append([line,tpm])
	else:
		# Save the whole line as being associated with this gene, along with the tpm count
		geneTssDict[geneName] = [[line,tpm]]


for geneName in geneTssDict:
	# Sort the lines for this gene by their TPM (col 1)
	sortedGeneLines = sorted(geneTssDict[geneName], key=lambda x:float(x[1]), reverse=True)

	# Go through lines and output them as long as the TPM sum to this point
	# is not above the threshold desired.
	tpmSumFraction = 0
	for geneLine in sortedGeneLines:
		tpm = geneLine[1]
		if (tpm < args.mintpm):
			break

		fields = geneLine[0].split()
		posFields = re.split('\:|\.\.|,', fields[0])
		#sys.stderr.write("posFields: " + posFields[0] + "," + posFields[1] + "," + posFields[2] + "," + posFields[3] + "\n")
		tss = (int(posFields[1]) + int(posFields[2])) / 2

		detailFields = fields[1].split(',')
		peakFields = detailFields[0].split('@')
		geneName = peakFields[1]
		enstID = ''
		hgncID = ''
		if (re.search('ENST', geneName) is not None):
			enstID = geneName
		else:
			hgncID = geneName

		peakNum = peakFields[0][1:]

		print '\t'.join(map(str, [fields[0], fields[1], posFields[0], tss, peakNum, enstID, hgncID, tpm]))
