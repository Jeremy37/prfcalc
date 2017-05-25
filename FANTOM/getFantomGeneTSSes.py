#!/nfs/users/nfs_j/js29/anaconda/bin/python
##
## getFantomGeneTSSes.py
## Created 2015-8-17 by Jeremy Schwartzentruber
##
## Given the FANTOM file of potential TSSes and the TPM (transcripts per million)
## of those for a set of samples, outputs a file with the TSSes above a given
## threshold for a given tissue type.
import argparse
import sys
import subprocess
import re
import os
import os.path

parser = argparse.ArgumentParser(description="Processes the FANTOM tissue TPM file to output a file with TSSes above a certain TPM threshold for a given tissue type.")
parser.add_argument("--fantomfile", type=str, required=True, metavar='FILE', help="hg19.cage_peak_tpm_ann.osc.txt.gz")
parser.add_argument("--threshold", type=float, default=0.0, metavar='FLOAT', help="Threshold TPM above which a TSS should be included")
parser.add_argument("--cols", type=str, required=True, metavar='STR', help="Comma-separated list of integer column indexes indicating which columns to combine (average) for TSS TPM levels, e.g. '34,35,36'")
parser.add_argument("--weights", type=str, required=True, metavar='STR', help="Comma-separated list of floating-point weights for columns, in the same order as --cols, e.g. '1.234,0.45,1.78'")
parser.add_argument("--genes", action='store_true', help="Only include genes (i.e. having ENSG ID)")
parser.add_argument("--proteincoding", action='store_true', help="Only include protein-coding genes")
parser.add_argument("--output", type=str, default="out", metavar='STR', help="Path and prefix for output files")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")
parser.add_argument("--norun", action='store_true', help="don't run external commands (e.g. bedtools)")
args = parser.parse_args()

sys.stderr.write("getFantomGeneTSSes.py --fantomfile {0} --threshold {1} --cols {2}{3} --output {4}{5}{6}{7}{8}\n".format(
					args.fantomfile,
					args.threshold,
					args.cols,
					" --weights " + args.weights if args.weights else "",
					args.output,
					" --genes" if args.genes else "",
					" --proteincoding" if args.proteincoding else "",
					" --debug" if args.debug else "",
					" --norun" if args.norun else ""))

ENST_ENSG_HGNC_FILE="ENST.ENSG81.HGNC.txt"
GENCODE_PC_IDS_FILE="gencode.v19.pc.ids.txt"
scriptDir = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(scriptDir)
os.environ["PATH"] += os.pathsep + scriptDir

tempFiles = []

def main():
	global tempFiles
	
	try:
		# First use awk to average the specified columns from the Fantom file for each TSS
		cols = args.cols.strip().split(',')
		awkColStr = None
		if args.weights:
			# Check that there are the same number of weights as cols
			weightStrs = args.weights.strip().split(',')
			if len(weightStrs) != len(cols):
				die("Expected the same number of values for --weights ({0}) as --cols ({1}).".format(args.weights, args.cols))
			weights = [float(x) for x in weightStrs]
			totalWeight = sum(weights)
			awkColSum = '+'.join(["${0}*{1}".format(cols[i], weightStrs[i]) for i in range(len(cols))])
			awkColStr = "(({0})/{1})".format(awkColSum, totalWeight)
		else:
			awkColStr = '+'.join(["$%s" % x for x in cols])
			awkColStr += "/{0}".format(len(cols))
			
		cmd = ""
		if re.search(".gz$", args.fantomfile) is not None:
			cmd = "gzip -cd {0} | ".format(args.fantomfile)
		else:
			cmd = "cat {0} | "
		baseTpmFname = args.output + ".tpm.txt"
		tempFiles.append(baseTpmFname)
		cmd += "awk 'BEGIN {{OFS=\"\t\"; FS=\"\t\"}} {{ if ($1 !~ /^#|^0/) {{ print $1,$2,({0}) }} }}' > {1}".format(awkColStr, baseTpmFname)
		docall(cmd)
		
		# Use secondary script get_tissue_tss_tpm.py to order peaks by TPM and split
		# the Fantom fields based on gene ID (or ENST, or HGNC ID)
		cmd = "cat " + baseTpmFname + " | get_tissue_tss_tpm.py --mintpm 0.0"
		if args.genes:
			cmd += " --genes"
		splitTpmFname = args.output + ".tpm.split.txt"
		tempFiles.append(splitTpmFname)
		cmd += " > " + splitTpmFname
		docall(cmd)
		if not os.path.exists(splitTpmFname) or os.path.getsize(splitTpmFname) <= 0:
			die("Failed!")
		
		# Add a new column that has the latest HGNC symbols for the genes. We keep the old
		# column of symbols because other tables that help us match IDs might have the old
		# symbols!
		splitTpmHgnc = args.output + ".tpm.split.hgncupdated.txt"
		tempFiles.append(splitTpmHgnc)
		getHgncCmd = "cut -f 7 " + splitTpmFname + " | get.updated.hgnc.symbols.pl --file stdin"
		cmd = "paste " + splitTpmFname + " <(" + getHgncCmd + ") | awk 'BEGIN {OFS=\"\t\";FS=\"\t\"} {print $3,$4,$5,$6,$7,$8,$9}' > " + splitTpmHgnc
		docall(cmd)
		if os.path.getsize(splitTpmHgnc) <= 0:
			die("Failed!")

		# Add in the HGNC symbols in the cases where it was missing
		# Match on col 4 of hash file, col 2 of scan file -- and
		# annotate col 5 of the hash file using col 3 of the scan file.
		tmpFile1 = args.output + ".tpm.split.fixed.1.txt"
		tempFiles.append(tmpFile1)
		cmd = ("hash.annotate.pl --hashFile " + splitTpmHgnc +
								" --scanFile " + ENST_ENSG_HGNC_FILE +
								" --matchColHash 4" +
								" --matchColScan 2" +
								" --annotateColHash 5" +
								" --annotateColScan 3" +
								" > " + tmpFile1)
		docall(cmd)
		if os.path.getsize(tmpFile1) <= 0:
			die("Failed!")
		
		# Add in ENST transcript IDs in col 4 of the hash file by matching on col 5 (HGNC ID)
		# in cases where the ENST were missing
		tmpFile2 = args.output + ".tpm.split.fixed.2.txt"
		tempFiles.append(tmpFile2)
		cmd = ("hash.annotate.pl --hashFile " + tmpFile1 +
								" --scanFile " + ENST_ENSG_HGNC_FILE +
								" --matchColHash 5 " +
								" --matchColScan 3 " +
								" --annotateColHash 4 " +
								" --annotateColScan 2 " +
								" > " + tmpFile2)
		docall(cmd)
		if os.path.getsize(tmpFile2) <= 0:
			die("Failed!")

		# Also update the ENST transcript IDs in col 4 using the updated HGNC symbols to catch
		# a few remaining missing cases. 
		tmpFile3 = args.output + ".tpm.split.fixed.3.txt"
		tempFiles.append(tmpFile3)
		cmd = ("hash.annotate.pl --hashFile " + tmpFile2 +
								" --scanFile " + ENST_ENSG_HGNC_FILE +
								" --matchColHash 7 " +
								" --matchColScan 3 " +
								" --annotateColHash 4 " +
								" --annotateColScan 2 " +
								" > " + tmpFile3)
		docall(cmd)
		if os.path.getsize(tmpFile3) <= 0:
			die("Failed!")

		# Add a final column which is the ENSG gene ID, but then reorder the columns
		# so that this is the 4th column of the output.
		tmpFileAnn = args.output + ".tpm.split.fixed.ann.txt"
		tempFiles.append(tmpFileAnn)
		cmd = ("hash.annotate.pl --hashFile " + tmpFile3 +
								" --scanFile " + ENST_ENSG_HGNC_FILE +
								" --matchColHash 4 " +
								" --matchColScan 2 " +
								" --annotateColHash 8 " +
								" --annotateColScan 1 " +
								" | awk 'BEGIN {OFS=\"\t\";FS=\"\t\"} {print $1,$2,$3,$8,$4,$5,$7,$6}' " +
								" | sortByChrPos.pl -f stdin -chrcol 1 -poscol 2 > " + tmpFileAnn)
		docall(cmd)
		if os.path.getsize(tmpFileAnn) <= 0:
			die("Failed!")
		
		tpmFileFinal = args.output + ".mintpm.%.1f.txt" % args.threshold
		cmd = "awk 'BEGIN {{OFS=\"\t\";FS=\"\t\"}} {{if ($4 ~ /ENSG/ && $8 >= {0}) {{print $1,$2,$4,$8}}}}' {1} | stripChrCol.pl > {2}".format(args.threshold, tmpFileAnn, tpmFileFinal)
		docall(cmd)
		if os.path.getsize(tpmFileFinal) <= 0:
			die("Failed!")
		
		if args.proteincoding:
			tempFiles.append(tpmFileFinal)
			tpmFileFinalPC = args.output + ".mintpm.%.1f.pc.txt" % args.threshold
			cmd = ("hashJoin.pl --scanFile " + tpmFileFinal + " --hashFile " + GENCODE_PC_IDS_FILE + " --colHashFile 1 --colScanFile 3 " +
					" > " + tpmFileFinalPC)
			docall(cmd)
			if os.path.getsize(tpmFileFinalPC) <= 0:
				die("Failed!")
		
		if args.debug:
			sys.stderr.write("Done!\n\n")
		
	finally:
		if not args.debug:
			for fname in tempFiles:
				os.remove(fname)



import gzip
def checkGzFile(infile, mode='r'):
	if re.search(".gz$", infile.name) is not None:
		close(infile)
		if re.search('b', mode) is None:
			mode += 'b'
		return gzip.open(fname, mode)
	return infile

def die(msg):
	sys.stderr.write(msg + "\n")
	exit(1)

def docall(cmd):
	if args.debug:
		sys.stderr.write(cmd + "\n")
	retval = 0
	if (not args.norun):
		retval = subprocess.call(cmd, shell=True, executable='/bin/bash')
	if (retval != 0):
		die("Error running command.")

def docheckoutput(cmd):
	if args.debug:
		sys.stderr.write(cmd + "\n")
	retval = ""
	if (not args.norun):
		retval = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	return retval



main()
