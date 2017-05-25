#!/nfs/users/nfs_j/js29/anaconda/bin/python
## Here we process the FANTOM TSS and enhancer TPM files to produce 2 output files per
# Roadmap epigenome, one for TSSes and one for enhancers. We use a custom input file that
# indicates which FANTOM tissues we mapped to which Roadmap epigenomes, and with what
# weights.
import argparse
import sys
import subprocess
import re
import os
import os.path

parser = argparse.ArgumentParser(description="Processes the FANTOM TSS and enhancer TPM files to produce 2 output files per Roadmap epigenome, one for TSSes and one for enhancers.")
parser.add_argument("--tpm", type=float, default=0.0, required=True, metavar='FLOAT', help="Minimum TPM expression for TSSes")
args = parser.parse_args()

FANTOM_TPM_FILE="./input/hg19.cage_peak_tpm_ann.osc.txt.gz"
FANTOM_ENH_FILE="./input/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz"

def main():
	with open("epigenome.tsscols.txt", 'r') as f:
		for line in f:
			# Skip commented lines
			if re.search("^\s*#|^\s*$", line) is not None:
				continue
			lineVals = line.strip().split("\t")
			if len(lineVals) < 3:
				die("Only {0} columns found, but 3 expected. Line:\n{1}".format(len(lineVals), line))
			[epigenomeID, cols, weights] = line.strip().split("\t")
			cmd = "./getFantomGeneTSSes.py --fantomfile {0} --threshold {1} --cols {2} --weights {3} --genes --output tss/{4}.fantomTSS".format(FANTOM_TPM_FILE, args.tpm, cols, weights, epigenomeID)
			docall(cmd)

	with open("epigenome.enhcols.txt", 'r') as f:
		for line in f:
			# Skip commented lines
			if re.search("^\s*#|^\s*$", line) is not None:
				continue
			lineVals = line.strip().split("\t")
			if len(lineVals) < 3:
				die("Only {0} columns found, but 3 expected. Line:\n{1}".format(len(lineVals), line))
			[epigenomeID, cols, weights] = line.strip().split("\t")
			cmd = "getFantomEnh.py --fantomfile {0} --cols {1} --weights {2} --debug > enh/{3}.fantomEnh.bed".format(FANTOM_ENH_FILE, cols, weights, epigenomeID)
			docall(cmd)

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
