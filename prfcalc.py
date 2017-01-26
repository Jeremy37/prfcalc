#!/software/python-2.7.10/bin/python
##
## prfcalc.py
## Created 2015-7-17 by Jeremy Schwartzentruber
##
import argparse
import sys
import os
import os.path
import gc
import re
import string
import gzip
import subprocess
import math
import pandas as pd
import numpy as np
import StringIO
import datetime
import bx.bbi.bigwig_file


parser = argparse.ArgumentParser(description="Given a set of SNPs and file listing annotations and their enrichments, calculates PRF scores for the SNPs.")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--snpbed", type=file, metavar='FILE',
					help="BED file of SNP positions to score")
group.add_argument("--range", type=str, help="Genomic range to get scores for, in the format 1:10000-11000")

group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument("--annotparams", type=file, metavar='FILE',
					help="file with list of annotations and their enrichments")
group2.add_argument("--annotparamfiles", type=file, metavar='FILE',
					help="file with list of annotparams files to use (i.e. names of files you would pass to --annotparams). ")
parser.add_argument("--tssfile", type=file, metavar='FILE',
					help="file with list of tss locations and associated genes")					
parser.add_argument("--bsub", action='store_true', help="use bsub to submit annotation jobs for --annotparamfiles separately")

parser.add_argument("--pergene", action='store_true',
					help="output PRF scores per gene in a separate output file. PRF 'singlescore' (max PRF score across genes) are still output per SNP.")
parser.add_argument("--precision", type=int, default=4, metavar='NUM',
					help="maximum number of significant digits to output")
parser.add_argument("--maxtssdist", type=int, default=1e6, metavar='NUM',
					help="maximum TSS dist for SNPs to be considered for a gene; default 1 million")
parser.add_argument("--annotate", type=str,
					help="output annotation values pre- and post-normalization. This argument specifies the base of the file path where annotations should be saved")
parser.add_argument("--qctest", type=str,
					help="run quality checks/testing and outputs in the file specified; assumes that the input snpbed has 5 columns, with the extra 5th col having the relevant gene ID")
parser.add_argument("--output", type=str, default="out", help="root name/path for output files")
parser.add_argument("--xoutput", action='store_true', help="output X values per annotation for the gene with highest PRF score")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")
parser.add_argument('--verbose', '-v', action='count')
parser.add_argument("--norun", action='store_true', help="don't run external commands (e.g. bedtools)")

#os.chdir('C:\\Users\\Jeremy\\Google Drive\\_Gaffney Lab\\ipython\\prfcalc')
#os.chdir('/Users/js29/Google Drive/_Gaffney Lab/ipython/prfcalc')
#os.chdir('/nfs/users/nfs_j/js29/src/prfcalc/test')
#argslist = ["--qctest", "snpQC.1k.x", "--snpbed", "geuvadis.pc.1e-6.full.snps.qc1k.txt", '--annotparams', 'annotparams.full.txt', '--debug']
#args = parser.parse_args(argslist)

args = parser.parse_args()


codingMaskBed = '/lustre/scratch117/cellgen/team170/js29/annotations/geneModels/gencode.v19.coding.bed'
transformsDict = {}
tempFiles = []
snpDF = None
snpLocChrDF = None
snpBedNoChr = None
snpBedChr = None
snpBedFile = None
paramsDF = None
paramVals = None
paramDistBins = None
snpGeneIDs = None
tssDict = {}
MAX_TSS_DIST = args.maxtssdist
#MAX_TSS_DIST = 5e5

def main():
	global tempFiles
	try:
		if args.range:
			tempSnpBedName = makeBedFromRange("tempSnpBed", args.range)
			tempFiles.append(tempSnpBedName)
			snpBedFile = openg(tempSnpBedName)
		else:
			snpBedFile = openg(args.snpbed)
		
		if args.annotparamfiles:
			handleAnnotParamFiles(snpBedFile, args.annotparamfiles)
		else:
			# args.annotparams should be set
			doPrfcalc(snpBedFile, args.annotparams, args.output)
			
		if args.verbose:
			if args.verbose > 1: sys.stderr.write(str(datetime.datetime.now())+"\n")
			sys.stderr.write("\nDone!\n")
	
	finally:
		for fname in tempFiles:
			os.remove(fname)


def handleAnnotParamFiles(snpBedFile, annotparamfiles):
	global tempFiles
	annotparamfileLines = [l.strip() for l in annotparamfiles.readlines()]
	if len(annotparamfileLines) < 1:
		die("No lines in --annotparamfiles file: {0}".format(annotparamfiles.name))
		
	if args.bsub:
		# First create the "chr" or "nochr" version of the Bed file in a single
		# process before multiple jobs are submitted, to avoid race condition
		readSnpFile(snpBedFile)
	
	for i,line in enumerate(annotparamfileLines):
		epigenomeID,annotparamFile,tssFile = line.split("\t")
		outputArg = args.output + (".{0}".format(epigenomeID))
		perGeneArg = "--pergene" if args.pergene else ""
		verboseArg = "-" + 'v' * args.verbose
		debugArg = "--debug" if args.debug else ""
		prfcalcCmd = "prfcalc.py --output {0} --snpbed {1} --tssfile {2} --annotparams {3} {4} {5} {6}".format(outputArg, snpBedFile.name, tssFile, annotparamFile, perGeneArg, verboseArg, debugArg)
		
		outputFname = args.output + ".{0}.singlescore.txt.gz".format(epigenomeID)
		if os.path.isfile(outputFname):
			sys.stderr.write("Output file {0} already exists. Skipping command:\n{1}\n".format(outputFname, prfcalcCmd))
			continue
			
		if args.bsub:
			jobname = "prfcalc.sub.{0}".format(epigenomeID)
			if not os.path.exists("farmOut"):
				os.makedirs("farmOut")
			farmOut = "farmOut/{0}.%J.txt".format(jobname)
			mem = 2000
			memoryString = '-R"span[hosts=1] select[mem>' + str(mem) +'] rusage[mem=' + str(mem) + ']" -M ' + str(mem)
			#bsubCmd = "bsub -G team170 -J {0} -o {1} -q normal {2} {3}".format(jobname, farmOut, memoryString, prfcalcCmd)
			bsubCmd = "bsub -J {0} -o {1} -q normal {2} {3}".format(jobname, farmOut, memoryString, prfcalcCmd)
			docall(bsubCmd)
		else:
			docall(prfcalcCmd)
			#with open(annotparamFile) as paramsFile:
			#	doPrfcalc(snpBedFile, paramsFile, args.output + (".%d" % i))

	if args.bsub:
		# Wait for all submitted jobs to finish, so that we can merge the results together
		jobname = "prfcalc.bsub.wait"
		waitCmd = "echo 'done'"
		memoryString = "-R\"span[hosts=1] select[mem>{0}] rusage[mem={0}]\" -M {0}".format(str(100))
		#bsubCmd = "bsub -G team170 -J {0} -o farmOut/{0}.%J.txt -q normal {1} -K -w \"ended('prfcalc.sub.*')\" {2}".format(jobname, memoryString, waitCmd)
		bsubCmd = "bsub -J {0} -o farmOut/{0}.%J.txt -q normal {1} -K -w \"ended('prfcalc.sub.*')\" {2}".format(jobname, memoryString, waitCmd)
		docall(bsubCmd)
		
	# All jobs should be done now. We can merge together the output files.
	header = "chr\tpos\tsnp"
	pasteCmd = "paste <(cut -f 1,3,4 {0})".format(snpBedFile.name)
	for i,line in enumerate(annotparamfileLines):
		epigenomeID,annotparamFile,tssFile = line.split("\t")
		outputFname = args.output + ".{0}.singlescore.txt.gz".format(epigenomeID)
		header += "\t{0}_PRF\t{0}_geneID".format(epigenomeID)
		pasteCmd += " <(zcat {0} | sed '1d' | cut -f 4,5)".format(outputFname)
		#getGenesCmd += " <(zcat {0} | cut -f 5)".format(outputFname)
		if not os.path.isfile(outputFname):
			die("Expected output file from prfcalc.py subcommand is missing: {0}".format(outputFname))
		if not args.debug:
			tempFiles.append(outputFname)
	
	mergedOutputFname = args.output + ".merged.singlescore.txt.gz"
	echoCmd = "echo -e '{0}' | gzip > {1}".format(header, mergedOutputFname)
	pasteCmd += " | gzip >> {0}".format(mergedOutputFname)
	docall(echoCmd)
	docall(pasteCmd)
	

def doPrfcalc(snpBedFile, annotparamsfile, outputRoot):
	if args.annotate:
		readSnpFile(snpBedFile)
		readParams(annotparamsfile, outputRoot)
		computeSnpAnnotations(outputRoot)
	
	elif args.qctest:
		readSnpFile(snpBedFile)
		readParams(annotparamsfile)
		computeSnpAnnotations(outputRoot)
		computeScoresQC()

	else:
		if not args.tssfile:
			die("Argument --tssfile is required.")
		readSnpFile(snpBedFile)
		readParams(annotparamsfile, outputRoot)
		readTssFile(args.tssfile, outputRoot)	
		computeSnpAnnotations(outputRoot)
		computeScores(outputRoot)


def isgzfile(fname):
	namelen = len(fname)
	return fname[namelen-3:namelen] == ".gz"

def openg(f, mode=None):
	if type(f) is file:
		if isgzfile(f.name):
			return gzip.GzipFile(fileobj=f)
		return f
	elif type(f) is str:
		if isgzfile(f):
			if mode is None:
				mode = 'rb'
			return gzip.open(f, mode)
		if mode is None:
			mode = 'r'
		return open(f, mode)
	else:
		die("openg: unrecognized type for parameter 'f': {0}".format(type(f)))
 

def makeBedFromRange(basename, rangeStr):
	# Create a temporary snpBed file to store positions for overlap with annotations
	snpRangeVals = re.split(":|-", rangeStr)
	# Check that the range values make sense
	chr = snpRangeVals[0]
	bedfilename = basename + "." + ".".join(snpRangeVals) + ".bed.gz"
	if os.path.isfile(bedfilename):
		sys.stderr.write("Creating bed file for range, but file already exists: {0}. New file will not be created.".format(bedfilename))
	else:
		with gzip.open(bedfilename, "w") as f:
			[startpos, endpos] = [int(x) for x in snpRangeVals[1:3]]
			for pos in range(startpos, endpos):
				f.write("{0}\t{1}\t{2}\t.\n".format(chr, pos, pos+1))
	return bedfilename
	

#@profile
def readSnpFile(snpBedFile):
	global snpDF
	global snpLocChrDF
	global snpBedNoChr
	global snpBedChr
	# Read in all SNPs, checking that they are in sorted order
	snpDF = pd.read_table(snpBedFile, header=None, dtype={0: str, 1: np.int32, 2: np.int32, 3: str} )
	if args.qctest:
		colnames = ['chr', 'start', 'end', 'name', 'gene', 'tssdist']
	else:
		colnames = ['chr', 'start', 'end', 'name']
		
	if len(snpDF.columns) == len(colnames):
		snpDF.columns = colnames
	else:
		die("SNP bed file has {0} columns, but {1} are expected".format(len(snpDF.columns), len(colnames)))
	
	# Check that the SNP bed file is in sorted order. If we are only annotating this isn't important.
	snpChrs = snpDF.chr.values
	snpPositions = snpDF.start.values
	for i in xrange(1, len(snpDF.index)):
		if snpChrs[i] == snpChrs[i-1] and int(snpPositions[i]) < int(snpPositions[i-1]):
			sys.stderr.write("Warning: SNP bed file coordinates are not in sorted order. Line number {0} found out of order:\n{1}\n".format(i, str(snpDF.iloc[i,:])))
			break
	if args.verbose: sys.stderr.write("Finished reading SNP bed file: {0}. Got {1} SNPs.\n".format(snpBedFile.name, len(snpDF.index)))
	
	# Make the alternate "chr" version of the SNP file, so that we have either format
	# to use later in overlapping with bed/bigwig annotations.
	hasChr = (re.search("chr", snpDF.chr[0]) is not None)
	if (hasChr):
		snpBedChr = snpBedFile.name
		if isgzfile(snpBedChr):
			namebase = os.path.splitext(snpBedChr)[0]
			snpBedNoChr = namebase + ".nochr.gz"
		else:
			snpBedNoChr = snpBedChr + ".nochr.gz"
		#tempFiles.append(snpBedNoChr)
		if not os.path.isfile(snpBedNoChr):
			if args.verbose: sys.stderr.write("Writing SNP bed file with 'chr' removed: " + snpBedNoChr + "\n")
			if isgzfile(snpBedChr):
				cmd = "zcat {0} | perl -ne 's/chr//g; print' | gzip > {1}".format(snpBedChr, snpBedNoChr)
			else:
				cmd = "cat {0} | perl -ne 's/chr//g; print' | gzip > {1}".format(snpBedChr, snpBedNoChr)
			docall(cmd)
		snpLocChrDF = snpDF.loc[:,['chr','start']]
		# Read in the version without 'chr' and save in snpDF
		snpDF = pd.read_table(snpBedNoChr, header=None, dtype={0: str, 1: np.int32, 2: np.int32, 3: str} )
		snpDF.columns = colnames
	else:
		snpBedNoChr = snpBedFile.name
		if isgzfile(snpBedNoChr):
			namebase = os.path.splitext(snpBedNoChr)[0]
			snpBedChr = namebase + ".chr.gz"
		else:
			snpBedChr = snpBedNoChr + ".chr.gz"
		#tempFiles.append(snpBedChr)
		if not os.path.isfile(snpBedChr):
 			if args.verbose: sys.stderr.write("Writing SNP bed file with 'chr' added: " + snpBedChr + "\n")
			if isgzfile(snpBedNoChr):
				cmd = "zcat {0} | perl -ne 'print \"chr\".$_' | gzip > {1}".format(snpBedNoChr, snpBedChr)
			else:
				cmd = "cat {0} | perl -ne 'print \"chr\".$_' | gzip > {1}".format(snpBedNoChr, snpBedChr)
			docall(cmd)
			# The original code below is incredibly slow (e.g. 40 min for 40 M SNPs)
			# Calling perl, in contrast, takes seconds
 			#snpDF_chr = snpDF.copy()
 			#for i in range(0, len(snpDF_chr.index)):
 			#	snpDF_chr.iat[i,0] = 'chr' + str(snpDF_chr.iat[i,0])
 			#if args.verbose: sys.stderr.write("Writing SNP bed file with 'chr' added: " + snpBedChr + "\n")
 			#snpDF_chr.to_csv(snpBedChr, sep='\t', header=False, index=False)
		# We save just the first two columns with 'chr' to use in overlapping bigwig files
		snpLocChrDF = pd.read_table(snpBedChr, header=None, dtype={0: str, 1: np.int32, 2: np.int32, 3: str} )
		snpLocChrDF.columns = colnames
		snpLocChrDF = snpLocChrDF.loc[:,['chr','start']]
		
	# snpDF itself now has the version without 'chr'


#@profile
def readParams(paramsFile, outputRoot):
	global paramsDF
	global paramVals
	global paramDistBins
	# Read all parameter values into an array. Also store an array that for each annotation
	# keeps the name and type of the file associated.
	if args.verbose: sys.stderr.write("Reading annotation params file: " + paramsFile.name + "\n")
	paramsDF = pd.read_table(paramsFile, header=None, dtype={0: str, 1: str, 2: str, 3: str, 4: str, 5: str}, keep_default_na=False)
	paramsDF.columns = ['paramName', 'params', 'file', 'type', 'transform', 'haschr']
	paramVals = list()
	paramDistBins = list()
	
	# Do sanity checks on the params file
	for index, row in paramsDF.iterrows():
		annotType = row.type
		if annotType == "tssdist":
			annotParams = [float(v) for v in row.params.split(',')]
			paramVals.append(annotParams)
			distBins = [ [float(v) for v in bin.split("-")] for bin in row.transform[8:].split(',') ]
			# Check that the distance bins are ordered with the farthest (and largest) bin first
			lastDist = 1e9
			for bin in distBins:
				if bin[1] > lastDist:
					die("tssdist bins in params file should be ordered with farthest bin first. Bin {0}-{1} is out of order.".format(bin[0], bin[1]))
			paramDistBins.append(distBins)
			
		else:
			# Check that each annotation file listed exists
			if not os.path.isfile(row.file):
				die("Annotation file \"{0}\" doesn't exist".format(row.file))
		
			# Save and check annotation enrichment params
			annotParams = None
			if re.search("bigwig|bedcol", annotType) is not None:
				annotParams = [float(v) for v in row.params.split(',')]
			else:
				annotParams = [float(row.params)]
			paramVals.append(annotParams)
			distBins = []
			if isinstance(row.transform, basestring) and row.transform[0:7] == "tssdist":
				distBins = [float(v) for v in row.transform[8:].split("-")]
			paramDistBins.append(distBins)
		
			if row.haschr:
				# Check that haschr is true/false
				if row.haschr.lower() == "t" or row.haschr.lower() == "true":
					paramsDF.loc[index,'haschr'] = True
				else:
					die("Unrecognized value for column 6 ('haschr') of params file:\n{0}\n".format(str(row)))
			else:
				# Haschr is not yet defined
				if annotType == 'bigwig':
					# First check whether the annotation file has "chr" in its chr field
					bigWigHead = docheckoutput("head -n 1 <(bigWigToWig {0} /dev/stdout)".format(row.file))
					## The first way I tried to do this is below, but it didn't work -- I guess because
					## of some interaction of subprocess with passing /dev/stdout to bigWigToWig?
					##bigWigHead = docheckoutput("bigWigToWig {0} /dev/stdout | head -n 1".format(row.file))
					hasChr = (re.search("chr", bigWigHead) is not None)
					paramsDF.loc[index,'haschr'] = hasChr
			
				elif annotType == 'bed' or annotType == 'diffgene' or annotType == 'samegene' or re.search('bedcol', annotType) is not None:
					# First check whether the annotation file has "chr" in its chr field
					if re.search('.gz$', row.file) is not None:
						cmd = "gzip -cd {0} | head -n 1".format(row.file)
					else:
						cmd = "head -n 1 {0}".format(row.file)
					bedHead = docheckoutput(cmd)
					bedHeadChr = bedHead.split(' ')
					hasChr = (re.search("chr", bedHeadChr[0]) is not None)
					paramsDF.loc[index,'haschr'] = hasChr
			
				else:
					die("Unrecognized annotation type: {0}".format(row.type))

		
	if args.verbose:
		sys.stderr.write("Param dist bins: " + str(paramDistBins) + "\n")
		sys.stderr.write("Param vals: " + str(paramVals) + "\n")
	if args.verbose > 2:
		paramsDF.to_csv(outputRoot + ".params.used.txt", sep='\t')
	

#@profile
def readTssFile(tssFile, outputRoot):
	global tempFiles
	global tssDict
	global snpGeneIDs
	tssChrs = {}
	
	if args.verbose:
		sys.stderr.write("Reading TSS file: " + tssFile.name + "\n")
	
	# First read in all TSSes. Check that they are sorted by position
	for linestr in tssFile:
		lineVals = linestr.strip().split('\t')
		#if (len(lineVals) < 3):
		#	continue
		#chr = lineVals[0]
		#pos = int(lineVals[1])
		geneID = lineVals[2]
		if (geneID in tssDict):
			tssDict[geneID].append(int(lineVals[1]))
		else:
			tssDict[geneID] = list([int(lineVals[1])])
			tssChrs[geneID] = lineVals[0]
	
	# Write a BED file that has, for each gene, a window extending 1 Mb from the lowest
	# TSS and 1 Mb from the highest TSS.
	tssTable = []
	rowindex = 0
	for geneID in tssDict:
		minTSS = min(tssDict[geneID])
		minPos = int(max(0, minTSS - MAX_TSS_DIST))
		maxPos = int(max(tssDict[geneID]) + MAX_TSS_DIST)
		tssTable.append([tssChrs[geneID], minPos, maxPos, geneID, minTSS])
	
	tssDF = pd.DataFrame(tssTable, columns=['chr','start','end','geneID','minTSS'])
	tssDF[['chr','geneID']] = tssDF[['chr','geneID']].astype(str)
	tssDF[['start','end','minTSS']] = tssDF[['start','end','minTSS']].astype(np.int32)
	
	generangeFname = outputRoot + '.tss.generange.bed'
	tssDF_sorted = tssDF.sort(['chr', 'minTSS'], ascending=[1, 1])
	tssDF_sorted[['chr','start','end','geneID']].to_csv(generangeFname, sep='\t', header=False, index=False)
	tempFiles.append(generangeFname)
	
	# Use the generated file with bedtools to get the set of genes overlapping each SNP
	hasChr = (re.search("chr", tssDF.chr[0]) is not None)
	snpBed = (snpBedChr if hasChr else snpBedNoChr)
	
	global snpDF
	cmd = "bedtools intersect -a {0} -b {1} -loj".format(snpBed, generangeFname)
	output = docheckoutput(cmd)
	snpGeneIDs = accumulateGeneOverlaps(output, len(snpDF.columns))
	gc.collect()
	numSnpRows = len(snpDF.index)
	if (len(snpGeneIDs) != numSnpRows):
		die("Unexpected number of gene ID values returned from accumulateGeneOverlaps. Got {0}, expected {1}. This can happen if there are duplicate SNP positions in the input. All SNP positions should be unique.".format(len(snpGeneIDs), numSnpRows))
	
	snpGenesDF = snpDF.copy()
	#snpGenesDF.loc[:,numSnpCols] = snpGeneIDs
	snpGenesDF.loc[:,'geneIDs'] = snpGeneIDs
	if args.verbose > 2:
		snpGenesDFFname = outputRoot + '.genesDF.txt'
		snpGenesDF.to_csv(snpGenesDFFname, sep='\t', index=False)
	

#@profile
def computeSnpAnnotations(outputRoot):
	# First initialize the array of SNPs x annotations
	global snpDF
	numSnps = len(snpDF.index)
	numSnpCols = len(snpDF.columns)
	numAnnot = len(paramsDF.index)
	columnLabels = paramsDF.iloc[:,0].values
	annotDF = pd.DataFrame([[0]*numAnnot]*numSnps, columns=columnLabels, dtype=np.float32)
	snpDF = pd.concat([snpDF, annotDF], axis=1)
	
	###snpDF.to_csv(outputRoot + ".snpDF.txt", sep='\t', header=False, index=False)
	###if args.debug: print snpDF.dtypes
	# First get a vector indicating for each SNP whether it is in the near, middle, or far
	# TSS distance bin
	
	if args.verbose:
		sys.stderr.write("\nGetting SNP annotations\n")
	
	# For each annotation, fill the appropriate column of values
	snpDFColIndex = numSnpCols-1
	for index, row in paramsDF.iterrows():
		if args.debug: sys.stderr.write("Param: " + str(row) + "\n")
		if args.verbose > 1: sys.stderr.write(str(datetime.datetime.now())+"\n")
		snpDFColIndex += 1
		snpBed = (snpBedChr if row.haschr else snpBedNoChr)			
	
		if row.type == 'bed':
			# For a binary (bed) annotation, use bedtools to overlap with all SNPs and store
			# the annot enrichment for each overlap
			enrichment = row.params
			countCol = numSnpCols + 1
			cmd = "bedtools intersect -c -a {0} -b {1} | awk -v enrichment=\"{2}\" '{{ if(${3}>=1) {{ ${3}=enrichment }}; print ${3} }}'".format(snpBed, row.file, enrichment, countCol)
			output = docheckoutput(cmd).strip()
			annotValues = np.float32(output.strip().split('\n'))
			if len(annotValues) != numSnps:
				die("Did not get correct number of SNP values for annotation {0}".format(row.paramName))
				
			snpDF.iloc[:,snpDFColIndex] = annotValues
			#f = open('output.{0}.txt'.format(row.paramName),'w')
			#f.write('\n'.join([("%f" % float(x)) for x in snpDF.iloc[:,snpDFColIndex].values]) + '\n')
			#f.close()
			
		elif row.type == 'bigwig':
			# For a quantitative (bigwig) annotation, use bigWigAverageOverBed to get annot values
			# tmpBedIn = row.paramName + ".bedin.tmp"
			# tmpBedOut = row.paramName + ".bedOut.tmp"
			# cmd = "cat {0} | awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,NR}}' > {1}".format(snpBed, tmpBedIn)
			# docall(cmd)
			# cmd = "(bigWigAverageOverBed {0} {1} output.junk.tmp -bedOut={1}.bedOut.tmp; cat {2}) | awk -v col={3} '{{print $col}}'".format(row.file, tmpBedIn, tmpBedOut, numSnpCols+1)
			# output = docheckoutput(cmd).strip()
			# os.remove(tmpBedIn)
			# os.remove(tmpBedOut)
			# os.remove("output.junk.tmp")
			# annotValues = np.float32(output.split('\n'))
			# if len(annotValues) != numSnps:
				# die("Did not get correct number of SNP values for annotation {0}".format(row.paramName))
			# snpDF.iloc[:,snpDFColIndex] = annotValues
			# f = open('output.{0}.txt'.format(row.paramName),'w')
			# f.write(output)
			# f.close()
			
			if row.haschr:
				annotValues = getBigWigValues(snpLocChrDF, row.file)
			else:
				annotValues = getBigWigValues(snpDF.loc[:,['chr','start']], row.file)
			
			# Gerp is special in that we set all coding bases to NA
			if row.paramName == 'GerpRS':
				if args.verbose:
					sys.stderr.write("Handling GerpRS annotation - masking coding values with NaN, using file: {0}\n".format(codingMaskBed))
				countCol = numSnpCols + 1
				cmd = "bedtools intersect -c -a {0} -b {1} | awk '{{ print ${2} }}'".format(snpBedNoChr, codingMaskBed, countCol)
				mask = np.int32(docheckoutput(cmd).strip().split('\n'))
				annotValuesTmp = annotValues
				annotValues = [np.nan if mask[k] else annotValuesTmp[k] for k in xrange(numSnps)]
				if len(annotValues) != len(annotValuesTmp):
					die("Did not get correct number of annotation values for GerpRS after masking")
				
			snpDF.iloc[:,snpDFColIndex] = annotValues
			#output = '\n'.join([("%f" % float(x)) for x in annotValues]) + '\n'
			#f = open('output.{0}.txt'.format(row.paramName),'w')
			#f.write(output)
			#f.close()
			
		elif row.type == 'diffgene' or row.type == 'samegene':
			# For a gencode annotation, use bedtools as above, but rather than storing the
			# the annot enrichment, store the gene ID of the overlapping gene. This lets us
			# later determine if the annotation corresponds to the "same gene" or not for
			# a given SNP-gene combination.
			cmd = "bedtools intersect -a {0} -b {1} -loj".format(snpBed, row.file)
			output = docheckoutput(cmd)
			geneIDs = accumulateGeneOverlaps(output, numSnpCols)
			gc.collect()
			if (len(geneIDs) != numSnps):
				die("Unexpected number of gene ID values returned from accumulateGeneOverlaps")
			snpDF.iloc[:,snpDFColIndex] = geneIDs
			
			#f = open('output.{0}.txt'.format(row.paramName),'w')
			#f.write('\n'.join(geneIDs) + '\n')
			#f.close()
			
		elif re.search('bedcol', row.type) is not None:
			# The index of the column to retrieve should appear right after 'bedcol'
			colIndex = int(row.type[6:]) + numSnpCols 
			if args.verbose: sys.stderr.write("Getting column index {0}\n".format(colIndex))
			
			cmd = "bedtools intersect -a {0} -b {1} -loj | awk -v col={2} '{{print $col}}'".format(snpBed, row.file, colIndex)
			output = docheckoutput(cmd).strip()
			outputArray = [np.nan if (x == '.' or x == "-1") else x for x in output.split('\n')]
			annotValues = np.float32(outputArray)
			if len(annotValues) != numSnps:
				die("Did not get correct number of SNP values for annotation {0}".format(row.paramName))
			
			snpDF.iloc[:,snpDFColIndex] = annotValues
			#f = open('output.{0}.txt'.format(row.paramName),'w')
			#f.write('\n'.join([("%7f" % x) for x in annotValues]) + '\n')
			#f.close()
		
		elif row.type == 'tssdist':
			# Do nothing; tssdist enrichments are determined later with respect to each gene
			# We store a small flag value so that the TSSdist annotation is nonzero... an
			# optimisation for when we calculate the scores
			snpDF.iloc[:,snpDFColIndex] = 1e-9
			
		else:
			die("Unrecognized annotation type: {0}".format(row.type))
			
		sys.stderr.flush()
	
	if args.annotate:
		snpDF.to_csv(outputRoot + ".snpDF.raw.txt", sep='\t', index=False, mode="w", float_format='%f')
	
	if args.debug:
		sys.stderr.write(str(snpDF.dtypes) + "\n")
	
	if args.verbose:
		if args.verbose > 1: sys.stderr.write(str(datetime.datetime.now())+"\n")
		sys.stderr.write("\nTransforming quantitative SNP annotations\n")
	
	# Transform quantitative annotation values. Transforms based on distance or same/diff
	# gene can only be done later with respect to a specific gene.
	snpDFColIndex = numSnpCols-1
	for index, row in paramsDF.iterrows():
		snpDFColIndex += 1
		if row.type == 'bigwig' or re.search('bedcol', row.type) is not None:
			if not os.path.isfile(row.transform):
				die("File does not exist: {0}\n".format(row.transform))
			else:
				if args.verbose > 1:
					sys.stderr.write("Transforming {0}\n".format(row.paramName))
					sys.stderr.flush()
				snpDF.iloc[:,snpDFColIndex] = quantileTransform(snpDF.iloc[:,snpDFColIndex].values, row.transform)
		
	if args.annotate:
		snpDF.to_csv(outputRoot + ".snpDF.transformed.txt", sep='\t', index=False, mode="w", float_format='%f')
	

#@profile
def getBigWigValues(snpDF, fname):
	if args.verbose:
		sys.stderr.write("Opening bigwig file: {0}\n".format(fname))
	with open(fname, "rb") as f:
		bwh = bx.bbi.bigwig_file.BigWigFile(f)
		outputValues = []
		
		# We want to get blocks of up to 10 M SNPs at once, on an individual chromosome.
		# To do this, we step through the list of SNPs to determine a set of SNPs within
		# 10 Mb of each other, and then get an array of all values within the block.
		blockChr = snpDF.chr[0]
		blockStartPos = snpDF.start[0]
		maxpos = blockStartPos
		blockStartIndex = 0
		i = 1
		if args.debug:
			sys.stderr.write("len(snpDF.index): {0}\n".format(len(snpDF.index)))
		
		snpPositions = snpDF.start.values
		snpChrs = snpDF.chr.values
		numSnps = len(snpDF.index)
		while i < numSnps:
			#snppos = snpDF.iat[i,1]
			snppos = snpPositions[i]
			# We have a starting genomic position saved in blockStartPos. We keep going
			# through the SNP array as long as SNPs are not more than 10 Mb past this, and
			# also not upstream.
			if snpChrs[i] != blockChr or snppos > (blockStartPos + 1e7) or snppos < blockStartPos:
				blockEndPos = maxpos + 1
				if args.verbose > 1:
					sys.stderr.write("Block indices: {0} - {1}\tBlock pos: {2}:{3} - {4}\n".format(blockStartIndex, i-1, blockChr, blockStartPos, blockEndPos))
				data = bwh.get_as_array(blockChr, blockStartPos, blockEndPos)
				if data is None:
					die("Failed to get data from bigWig file. One reason this can happen is if the chromosome requested is not in the bigWig.")
					
				# data now has an array of the bigwig value for every genomic position
				# from blockStartPos to blockEndPos. So to get the values we use the position
				# offsets from blockStartPos as the index into the data array
				indices = [(pos - blockStartPos) for pos in snpPositions[blockStartIndex:i]]
				outputValues.extend([(0.0 if np.isnan(data[j]) else data[j]) for j in indices])
				# Note that bigWigAverageOverBed (which we used to prep annots in training) returns
				# 0 at bigWig file locations that are undefined, so this is what I also do here.
				
				# Reset variables for the start of a new block
				blockChr = snpChrs[i]
				blockStartPos = snppos
				maxpos = snppos
				blockStartIndex = i
			else:
				if snppos > maxpos:
					maxpos = snppos

			i += 1
		# Handle last block

		blockEndPos = maxpos+1
		if args.verbose > 1:
			sys.stderr.write("Block indices: {0} - {1}\tBlock pos: {2}:{3} - {4}\n".format(blockStartIndex, i-1, blockChr, blockStartPos, blockEndPos))
		data = bwh.get_as_array(blockChr, blockStartPos, blockEndPos)
		if data is None:
			die("Failed to get data from bigWig file. One reason this can happen is if the chromosome requested is not in the bigWig.")
		indices = [(pos - blockStartPos) for pos in snpDF.start[blockStartIndex:i]]
		outputValues.extend([(0.0 if np.isnan(data[j]) else data[j]) for j in indices])
	return outputValues	
	

#@profile
def accumulateGeneOverlaps(bedtoolsOutput, numSnpCols):
	# We assume that the gene ID was the 4th column of the bed file used to intersect with
	# SNPs, as is the case for the Gencode bed files I use
	geneIDcol = numSnpCols + 3
	geneIDs = []
	curGeneIDs = None
	lastSNP = None
	for line in StringIO.StringIO(bedtoolsOutput):
		lineVals = line.strip().split('\t')
		snp = '\t'.join(lineVals[0:4])
		if (lastSNP is None):
			lastSNP = snp
			curGeneIDs = lineVals[geneIDcol]
		elif (lastSNP == snp):
			curGeneIDs += "," + lineVals[geneIDcol]
		else:
			# Store the previous SNP's geneIDs, as this is a different SNP
			geneIDs.append(curGeneIDs)
			curGeneIDs = lineVals[geneIDcol]
			lastSNP = snp
	geneIDs.append(curGeneIDs)
	return geneIDs


# To efficiently do a quantile transform for many values (that is based on externally
# determined quantiles), we take as input a file specifying intervals that map to a
# normal quantile. An example of such a file is:
# 0		-1.96
# 5		-0.7
# 10	0.9
# NA	1.8
# This would indicate that values <= 0 map to -1.96; values from >0-5 map to -0.7;
# values from >5-10 map to 0.9, and values above 10 map to 1.8. It doesn't matter what
# the value is (first column) for the last entry in the table, as any values larger than
# the previous value (in this case 10) map to the largest quantile anyway.
#
#@profile
def quantileTransformOld(annotValues, transformFName):
	global transformsDict
	# Get the dataframe of quantile intervals and corresponding quantile values for this annot
	transformTable = pd.read_table(transformFName, header=None, dtype=np.float32).values
	searchVals = transformTable[:,0]
	transformVals = transformTable[:,1]
	transformedValues = [np.nan] * len(annotValues)
	
	# Function to get transformed annotation value, while leaving NaNs alone
	getTransform = lambda x: np.nan if np.isnan(x) else transformVals[binarySearch(searchVals, x)] 
	
	# Transform the first value in the array
	transformedValues[0] = getTransform(annotValues[0])
	for i in xrange(1, len(annotValues)):
		if annotValues[i] == annotValues[i-1]:
			transformedValues[i] = transformedValues[i-1]
		else:
			transformedValues[i] = getTransform(annotValues[i])
	return transformedValues
	
#@profile
def quantileTransform(annotVals, transformFName):
	global transformsDict
	# We first sort the values, because this speeds up the transform. Many of the same
	# values will occur one after the other, which means we don't have to redo the lookup.
	annotValDF = pd.DataFrame(annotVals, columns=['val'])
	annotValDFSorted = annotValDF.sort(columns=['val'], inplace=False)
	annotValues = annotValDFSorted.values
	numValues = len(annotValues)
	annotValDFSorted['transformedIndex'] = 0
	annotValDFSorted['transformedIndex'] = range(numValues)
	#sys.stderr.write("annotValDFSorted: {0}\n".format(str(annotValDFSorted)))
	#annotValDFSorted.to_csv(args.output + ".annotValDFSorted.1.txt", sep='\t', na_rep="NA", float_format='%.3f')
	
	# Get the dataframe of quantile intervals and corresponding quantile values for this annot
	transformTable = pd.read_table(transformFName, header=None, dtype=np.float32).values
	searchVals = transformTable[:,0]
	transformVals = transformTable[:,1]
	transformedValues = np.array([np.nan] * numValues)
	# Function to get transformed annotation value, while leaving NaNs alone
	getTransform = lambda x: np.nan if np.isnan(x) else transformVals[binarySearch(searchVals, x)] 
	
	# Transform the first value in the array
	transformedValues[0] = getTransform(annotValues[0])
	for i in xrange(1, numValues):
		if annotValues[i] == annotValues[i-1]:
			transformedValues[i] = transformedValues[i-1]
		else:
			transformedValues[i] = getTransform(annotValues[i])
	
	annotValDFSorted = annotValDFSorted.sort_index()
	#annotValDFSorted.to_csv(args.output + ".annotValDFSorted.2.txt", sep='\t', na_rep="NA", float_format='%.3f')
	#sys.stderr.write("annotValDF: {0}\n".format(str(annotValDF.iloc[0:10])))
	#sys.stderr.write("annotValDFSorted: {0}\n".format(str(annotValDFSorted.iloc[0:10])))
	#sys.stderr.write("transformedValuesSorted: {0}\n".format(str(transformedValues[0:10])))
	#sys.stderr.write("transformedValues: {0}\n".format(str(transformedValues[ [int(x) for x in annotValDFSorted.transformedIndex[0:10]] ])))
	
	# Use the sorted index of the dataframe to get back the (transformed)
	# values as per the original annotation value ordering.
	return transformedValues[ [int(x) for x in annotValDFSorted.transformedIndex] ]

	
#@profile
def binarySearch(vec, value):
	#return linearSearch(vec, value)
	imin = 0
	imax = len(vec)-1
	while imax > imin:
		imid = (imin + imax) / 2
		if value <= vec[imid]:
			imax = imid
		else:
			imin = imid+1
	return imax # Same as imin now


def linearSearch(vec, value):
	for i in xrange(len(vec)-1):
		if value < vec[i]:
			return i
	# The value was larger than all entries in the table, so return the maximum index
	return len(vec)-1


#@profile
def calcQuantAnnot(lambdaE, b0, b1, annotval):
	exponent = -b1 * (annotval - b0)
	# Threshold -- if 1/(1+exp(-x)) is <0.01 or >0.99, return 0 or 1.
	# This corresponds to exponent of +/- 4.59512
	# This greatly reduces the number of calls to exp for most annotations
	# without significantly changing the results.
	x = 0
	if (exponent < -4.59512):
		x = lambdaE
	elif (exponent < 4.59512):
		x = lambdaE / (1 + math.exp(exponent))
	return x


def getSnpMinDist(geneID, snpPos):
	global tssDict
	return min([abs(x-snpPos) for x in tssDict[geneID]])


#@profile
def computeScores(outputRoot):
	global snpDF
	global snpGeneIDs
	numSnps = len(snpDF.index)
	numAnnot = len(paramsDF.index)
	numSnpCols = len(snpDF.columns) - numAnnot
	
	if args.verbose:
		if args.verbose > 1: sys.stderr.write(str(datetime.datetime.now())+"\n")
		sys.stderr.write("\nComputing scores\n")
	if args.debug:
		snpDF.to_csv(outputRoot + ".snpDF.computeScores.1.txt", sep='\t', index=False, na_rep="NA", float_format='%.3f')
	
	# First compute the individual contributions (x value) of quantitative annotations
	# and store in snpDF
	snpDFColIndex = numSnpCols-1
	for index, row in paramsDF.iterrows():
		snpDFColIndex += 1
		###sys.stderr.write("index: {0}, snpDFColIndex: {1}, row: {2}\n".format(index, snpDFColIndex, str(row)))
		if row.type == 'bigwig' or re.search('bedcol', row.type) is not None:
			params = paramVals[index]
			snpDF.iloc[:,snpDFColIndex] = [calcQuantAnnot(params[0], params[1], params[2], annotval) for annotval in snpDF.iloc[:,snpDFColIndex]]
			
	if args.debug:
		snpDF.to_csv(outputRoot + ".snpDF.computeScores.2.txt", sep='\t', index=False, na_rep="NA", float_format='%.3f')
	
	# Make a list of the annotations that need updating for each SNP. These are either
	# tssdist-dependent annotations or geneID-dependent annotations.
	geneDepAnnots = []
	snpDFColIndex = numSnpCols-1
	for annotindex, row in paramsDF.iterrows():
		if row.type == 'tssdist':
			geneDepAnnots.append([annotindex, row.type])
			tssDistBins = paramDistBins[annotindex]
		elif row.type == 'bed':
			if row.transform:
				geneDepAnnots.append([annotindex, row.type])
		elif row.type == 'diffgene' or row.type == 'samegene':
			geneDepAnnots.append([annotindex, row.type])
	if args.verbose:
		sys.stderr.write("Gene-dependent annotations:\n" + ','.join(paramsDF.iloc[annot[0],:].paramName for annot in geneDepAnnots) + "\n")
	
	giantSnpGeneDF = None
	if args.debug:
		giantSnpGeneDF = pd.DataFrame(None, columns=snpDF.columns, dtype='float')
		giantSnpGeneDF['X'] = 0
		giantSnpGeneDF['GeneIDs'] = ''
		#np.append(giantSnpGeneArray, [np.round(np.random.random(8), 3)], axis=0)
	
	# Use numpy array snpDF.values as it is MUCH MUCH faster than using a pandas DataFrame
	snpVals = snpDF.values
	snpXList = []
	
	# For each SNP
	geneIDList = None
	genesStr = None
	###sys.stderr.write("snpDF:\n" + str(snpDF) + "\n")
	for i in xrange(len(snpDF.index)):
		snpPos = snpVals[i, 1]
		# snpGeneIDs[i] has a list of genes, e.g. "ENSG0001,ENSG0002". Split this to get
		# a list of individual gene IDs. In many cases the set of genes for one SNP will
		# be the same as for the next, so we don't need to do the split.
		if snpGeneIDs[i] != genesStr:
			genesStr = snpGeneIDs[i]
			if genesStr:
				geneIDList = genesStr.split(',')
			else:
				geneIDList = ""
		
		# For each gene, get the minimum distance of the SNP to any of the gene's TSSes
		# We use the 1-based SNP position as that is what was used while training the model
		# FANTOM TSSes are not given as individual base positions anyway, but a range usually
		# covering 10 - 20 bp.
		#sys.stderr.write("genesStr:'{0}'".format(genesStr))
		if not geneIDList:
			# There are no genes within range of this SNP! So we don't output anything
			# for the SNP
			#if args.debug:
			#	sys.stderr.write("Warning: no gene within range of snp {0}:{1}\n".format(snpVals[i, 0], snpVals[i,1]))
			numGenes = 0
			tssdistList = [1e6]
			#snpX = np.array([0] * numAnnot)
			#sys.stderr.write("HERE1:" + str(geneIDList) + "\n")
			#continue
		else:
			numGenes = len(geneIDList)
			tssdistList = [getSnpMinDist(geneID, snpVals[i, 2]) for geneID in geneIDList]
		
		# Duplicate the annotation values once for each gene. This is actually faster as
		# it lets us modify the array of values and use sum() for each line. Now each row
		# has all the annotations for a given gene (for this SNP)
		snpAnnot = snpVals[i, numSnpCols:].copy()
		snpXRow = snpAnnot.copy()
		for annot in geneDepAnnots:
			# Set all gene-dependent annot values initially to zero, so that we only actually SET
			# the smaller number of non-zero values (for efficiency)
			annotIndex = annot[0]
			snpXRow[annotIndex] = 0
		
		snpX = np.array([snpXRow] * max(numGenes, 1))
		
		# Update all gene-dependent annotations
		# annot[0] is the index of the annotation (in snpX), and annot[1] is the annot row from paramsDF
		for annot in geneDepAnnots:
			annotIndex = annot[0]
			###sys.stderr.write("annotIndex: " + str(annotIndex) + "\n")
			if not snpAnnot[annotIndex]:
				continue
			
			if annot[1] == 'tssdist':
				for j in range(numGenes):
					tssdist = tssdistList[j]
					# Assume that the tssdist bins are ordered with the farthest bin first
					# This is more efficient since distal bins are much larger
					if tssdist > tssDistBins[0][1]:
						# Farther than greatest tssdist bin, so enrichment is zero
						snpX[j, annotIndex] = 0
					else:
						for k,distbin in enumerate(tssDistBins):
							if distbin[0] <= tssdist and tssdist <= distbin[1]:
								snpX[j, annotIndex] = paramVals[annotIndex][k]
								break
				
			elif annot[1] == 'bed':
				[startdist, enddist] = paramDistBins[annotIndex]
				# The snpX array already has the enrichment for every overlap of a SNP
				# with the annotation. 
				for j in range(numGenes):
					#tssdist = tssdistList[j]
					if startdist <= tssdistList[j] and tssdistList[j] <= enddist:
						snpX[j, annotIndex] = snpAnnot[annotIndex]
				
			elif annot[1] == 'samegene':
				for j in range(numGenes):
					#if geneIDList[j] == snpAnnot[annotIndex]:
					if re.search(geneIDList[j], snpAnnot[annotIndex]) is not None:
						snpX[j, annotIndex] = paramVals[annotIndex][0]
				
			elif annot[1] == 'diffgene':
				for j in range(numGenes):
					#if geneIDList[j] != snpAnnot[annotIndex]:
					if re.search(geneIDList[j], snpAnnot[annotIndex]) is None:
						snpX[j, annotIndex] = paramVals[annotIndex][0]
				
			else:
				die("The list of gene-dependent annotations has an unrecognized annotation type")
	
		# Now that all annotation enrichment values have been calculated, sum them to get
		# the final pi value for this SNP
		#snpDF.loc[i,'X'] = ','.join(["%.3f" % x for x in sumX])
		#snpDF.loc[i,'GeneIDs'] = snpGeneIDs[i]
		
		sumX = np.sum(snpX, axis=1)
		snpXList.append(','.join(["%.3f" % x for x in sumX]))
		
		if args.xoutput:
			snpDF.iloc[i, numSnpCols:] = snpX[np.argmax(sumX)]
		
		if giantSnpGeneDF is not None:
			numGenes = max(numGenes, 1)
			commonSnpCols = np.array([snpVals[i,0:numSnpCols]]*numGenes)
			sumXCol = sumX.reshape((numGenes,1))
			geneIDCol = np.reshape(geneIDList, (numGenes,1))
			# Change the TSSDist column to store both the distance and the enrichment
			# Currently assuming that the TSSDist column is the first one!
			snpX[:,0] = ["{},{}".format(a, b) for a, b in zip(tssdistList, snpX[:,0])]
			curSnpGenesDF = pd.DataFrame(np.hstack((commonSnpCols, snpX, sumXCol, geneIDCol)), columns=giantSnpGeneDF.columns)
			giantSnpGeneDF = giantSnpGeneDF.append(curSnpGenesDF, ignore_index=True)
	
	if giantSnpGeneDF is not None:
		fname = outputRoot + ".giantSnpGeneDF.txt"
		sys.stderr.write("Saving table of annotation enrichments for all SNP-Gene combinations in {0}\n".format(fname))
		paramTypes = snpDF.dtypes
		for i,t in enumerate(paramTypes):
			if giantSnpGeneDF.columns[i] != "TSSDist":
				giantSnpGeneDF.iloc[:,i] = giantSnpGeneDF.iloc[:,i].astype(t)
		#giantSnpGeneDF['X'] = giantSnpGeneDF['X'].astype(float)
		giantSnpGeneDF.to_csv(fname, sep='\t', index=False, na_rep="NA", float_format='%f')
	
	snpDF.loc[:,'X'] = snpXList
	snpDF.loc[:,'GeneIDs'] = snpGeneIDs
	
	# For each SNP, get the max PRF score across genes
	# Be sure to properly handle/ignore cases where there is no X value for a SNP, which
	# happens when there is no gene in range of the SNP
	floatXValues = [[(float(x) if len(x) > 0 else np.nan) for x in xStr.split(',')] for xStr in snpXList]
	maxXValues = [max(floatXs) for floatXs in floatXValues]
	maxXIndexes = [floatXValues[i].index(maxXValues[i]) for i in range(0,len(maxXValues))]
	snpDF.loc[:,'maxX'] = maxXValues
	snpGeneIDsSeparated = [idStr.split(',') for idStr in snpGeneIDs]
	snpDF.loc[:,'maxGeneID'] = [snpGeneIDsSeparated[i][maxXIndexes[i]] for i in range(0,len(maxXValues))]
	
	snpDF.rename(columns={ 'end' : 'pos'}, inplace=True)
	with gzip.open(outputRoot + ".singlescore.txt.gz", 'wb') as f:
		colsToOutput = ['chr','pos','name','maxX','maxGeneID']
		if args.xoutput:
			x_cols_old = list(snpDF.columns.values[numSnpCols:len(snpDF.columns.values)])
			for col in ['maxX','maxGeneID','X','GeneIDs']: x_cols_old.remove(col)
			col_rename_hash = {colname:'x.' + colname for colname in x_cols_old}
			snpDF.rename(columns=col_rename_hash, inplace=True)
			colsToOutput.extend(['x.' + colname for colname in x_cols_old])
		snpDF.to_csv(f, sep='\t', columns=colsToOutput, index=False, header=True, na_rep="NA", float_format='%.3f')
	
	if args.pergene:
		with gzip.open(outputRoot + ".genescores.txt.gz", 'wb') as f:
			snpDF.to_csv(f, sep='\t', columns=['chr','pos','name','maxX','maxGeneID','X','GeneIDs'], index=False, header=True, na_rep="NA", float_format='%.3f')
	

# computeScoresQC is used with a different snpDF format than normal. Rather than having
# unique SNPs in a region, it is similar to the input file format used when training fgwas,
# i.e. each SNP may appear multiple times, in the cis-windows around different genes. The
# gene ID is the last fifth column of snpDF, and tssDist is the sixth.
# Thus, we go through the snpDF in order and do not iterate over genes explicitly.
#@profile
def computeScoresQC():
	global snpDF
	numSnps = len(snpDF.index)
	numAnnot = len(paramsDF.index)
	numSnpCols = len(snpDF.columns) - numAnnot
	
	sys.stderr.write("computeScoresQC()\n")
	
	# First just get the annotation values for all annotations, taking gene and tssdist
	# into consideration
	
	snpDF_annot = snpDF.values
	paramTypes = paramsDF.type.values
	transforms = paramsDF.transform.values
		
	# For each SNP
	for i in xrange(len(snpDF.index)):
		for j in range(len(paramsDF.index)):
			annotIndex = j + numSnpCols
			#paramVals[j]
			paramType = paramTypes[j]
			#transform = transforms[j]
			if paramType == 'bed':
				# Convert annotation back to 0/1, which is what we have stored in the annotation
				# input file during training. This facilitates comparison.
				annotval = snpDF_annot[i, annotIndex]
				if annotval:
					if len(paramDistBins[j]) == 2:
						[startDist, endDist] = paramDistBins[j]
						tssDist = snpDF_annot[i, 5]
						if startDist <= tssDist and tssDist <= endDist:
							snpDF_annot[i, annotIndex] = 1
						else:
							snpDF_annot[i, annotIndex] = 0
					else:
						snpDF_annot[i, annotIndex] = 1
				else:
					snpDF_annot[i, annotIndex] = 0
					
			elif paramType == 'bigwig' or re.search('bedcol', paramType) is not None:
				# Do nothing - transformed annot value already stored in snpDF
				continue
				#enrichment = calcQuantAnnot(paramVals[j][0], paramVals[j][1], paramVals[j][2], snpDF_annot[i, annotIndex]):
				#snpDF_annot[i, annotIndex] = enrichment
				
			elif paramType == 'diffgene':
				if snpDF_annot[i, annotIndex]:
					#if re.search(snpDF_annot[i,4], snpDF_annot[i, annotIndex]) is not None:
					if snpDF_annot[i,4] == snpDF_annot[i, annotIndex]:
						# This gene matches
						snpDF_annot[i, annotIndex] = 0
					else:
						# It must be a different gene matching
						snpDF_annot[i, annotIndex] = 1
				else:
					snpDF_annot[i, annotIndex] = 0
			elif paramType == 'samegene':
				if snpDF_annot[i, annotIndex] and snpDF_annot[i,4] == snpDF_annot[i, annotIndex]:
					snpDF_annot[i, annotIndex] = 1
				else:
					snpDF_annot[i, annotIndex] = 0

	# Make sure 0/1 annotation columns are represented and output as int
	sys.stderr.write("Writing annotated DF for QC\n")
	for j in range(len(paramsDF.index)):
		annotIndex = j + numSnpCols
		paramType = paramsDF.type[j]
		if paramType == 'bed' or paramType == 'diffgene' or paramType == 'samegene':
			snpDF_annot[:,annotIndex] = snpDF_annot[:,annotIndex].astype(int)
	
	snpDF.loc[:,:] = snpDF_annot
	snpDF.to_csv(args.qctest + ".annotated.txt", sep='\t', index=False, na_rep="NA", float_format='%.3f')


def die(msg):
	sys.stderr.write(msg + "\n")
	exit(1)

def docall(cmd):
	if args.verbose:
		sys.stderr.write(cmd + "\n")
	retval = 0
	if (not args.norun):
		retval = subprocess.call(cmd, shell=True, executable='/bin/bash')
	if (retval != 0):
		die("Error running command.")

def docheckoutput(cmd):
	if args.verbose:
		sys.stderr.write(cmd + "\n")
	retval = ""
	if (not args.norun):
		retval = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	return retval


main()

#import cProfile
#cProfile.run('main()', 'profile.stats')
