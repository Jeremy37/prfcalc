#!/nfs/users/nfs_j/js29/anaconda/bin/python
##
## getFantomEnh.py
## Created 2015-8-20 by Jeremy Schwartzentruber
##
## Given the FANTOM file of enhancers and their TPM (transcripts per million)
## outputs a file with the enhancer TPM for an "epigenome", i.e. a combined set
## of columns with weights.
import argparse
import sys
import subprocess
import re
import os
import os.path

parser = argparse.ArgumentParser(description="Processes the FANTOM enhancer TPM file to output a bed file with positions and TPMs.")
parser.add_argument("--fantomfile", type=str, required=True, metavar='FILE', help="human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz")
parser.add_argument("--cols", type=str, required=True, metavar='STR', help="Comma-separated list of integer column indexes indicating which columns to combine (average) for TSS TPM levels, e.g. '34,35,36'")
parser.add_argument("--weights", type=str, required=True, metavar='STR', help="Comma-separated list of floating-point weights for columns, in the same order as --cols, e.g. '1.234,0.45,1.78'")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")
args = parser.parse_args()

sys.stderr.write("getFantomEnh.py --fantomfile {0} --cols {1}{2}{3}\n".format(
					args.fantomfile,
					args.cols,
					" --weights " + args.weights if args.weights else "",
					" --debug" if args.debug else ""))

def main():
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
	cmd += "awk 'BEGIN {{OFS=\"\\t\"; FS=\"\\t\"}} {{ if ($1 !~ /^#|^0/) {{ print $1,({0}) }} }}'".format(awkColStr)
	cmd += " | perl -ne '@cols=split(/\\t/); ($chr,$start,$end)=split(/:|-/,$cols[0]); if ($cols[1] > 0.0) {print join(\"\\t\",$chr,$start,$end,$cols[1]);}'"
	docall(cmd)
	
	if args.debug:
		sys.stderr.write("Done!\n\n")
		



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
	retval = subprocess.call(cmd, shell=True)
	if (retval != 0):
		die("Error running command.")

def docheckoutput(cmd):
	if args.debug:
		sys.stderr.write(cmd + "\n")
	retval = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	return retval



main()
