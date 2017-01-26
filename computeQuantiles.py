#!/nfs/users/nfs_j/js29/anaconda/bin/python
##
## computeQuantiles.py
## Created 2015-7-17 by Jeremy Schwartzentruber
##
import argparse
import sys
import numpy as np
import pandas as pd
import scipy.stats
import math

parser = argparse.ArgumentParser(description="Given a set of values, does a quantile normal transformation of the values, and then outputs a table (with specified number of bins) that can be used to do an approximate transform on the same scale using another set of values. NA values are ignored.")

parser.add_argument("--input", type=file, metavar='FILE',
					help="file of values, one per line")
parser.add_argument("--numbins", type=int, required=True, metavar='FILE',
					help="number of bins to output for quantiles")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")

#argslist = ["--input", "test.effect-snp.nmotifs.txt", '--numbins', '1000', '--debug']
#args = parser.parse_args(argslist)
args = parser.parse_args()

inputFile = args.input
if not inputFile or inputFile == "stdin" or inputFile == "-":
	inputFile = sys.stdin

dataDF = pd.read_table(inputFile, header=None, dtype=float)
if len(dataDF.columns) > 1:
	sys.stderr.write("Input file has {0} columns but expected only one.".format(len(dataDF.columns)))
	exit(1)

# We only consider values that are non NA
numValues = len(dataDF.index)
notNA = [not math.isnan(x) for x in dataDF.iloc[:,0]]
n = sum(notNA)
if args.debug:
	sys.stderr.write("Found {0} NA values out of {1} total values input.\n".format(numValues-n, numValues))

# Set up a dataframe with the values as well as an integer index from 0:n that leaves
# out NA values
datavecNoNAs = pd.DataFrame(dataDF[notNA].values, index=np.arange(n), columns=['val'])

# Get the 0-based rank of each value. There may be duplicate ranks in the case of tied values.
ranks = scipy.stats.rankdata(datavecNoNAs.val, method='average').astype(int) - 1

# Get normal distribution values at evenly-spaced quantiles (i.e. if there are n non-NA
# input values, we want n unique normal Z scores at the n spaced quantiles).
# We start at n/(n+1) and go down to 1/(n+1) because rankdata() gives low ranks to the
# highest values, whereas we want the smallest value to have rank 0.
u = scipy.stats.norm.isf(np.linspace(float(n)/(n+1), float(1)/(n+1), n))

# Assign the normal Z score based on the quantile for the rank of each value
# There may be duplicated scores when values are tied
datavecNoNAs['qnorm'] = u[ranks]
datavecNoNAs['rank'] = ranks

if args.debug:
	datavecNoNAs.to_csv("computQuantiles.datavecNoNAs.txt", index=False, na_rep="NaN", sep='\t')

# First sort the dataframe by rank. Then we can just output the rows at the appropriate
# indices specified by binIndexes.
datavecNoNAs.sort(columns='rank', inplace=True)
dfNodup = datavecNoNAs.drop_duplicates()

numUnique = len(dfNodup.index)

# args.numbins specifies the number of bins we should output. Usually this will be less
# than the number of input values, but this doesn't matter. This gives us an index into
# the normal quantiles that we want to output. E.g. if we had only 3 data points as input
# our quantiles would be 0.25,0.5,0.75, and the Z scores would be -0.67, 0, 0.67. If we
# then had 5 bins, binIndexes would be [0,0,1,1,2]
numbins = args.numbins if args.numbins <= numUnique else numUnique
binIndexes = np.linspace(0, numUnique-1, numbins, dtype=int)

if args.debug:
	sys.stderr.write("There are {0} unique annotation values.\n".format(numUnique))
	sys.stderr.write("Outputting {0} non-duplicate bins out of total {1} bins requested.\n".format(len(binIndexes), args.numbins))
	
sys.stdout.write(dfNodup.iloc[binIndexes].drop_duplicates().to_csv(header=False, index=False, na_rep="NaN", sep='\t'))

if args.debug:
	dfNodup.to_csv("computQuantiles.dfNodup.txt", index=False, na_rep="NaN", sep='\t')
