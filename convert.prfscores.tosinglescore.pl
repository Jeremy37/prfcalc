#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

sub usage($);
my $usage = "$0 - Convert PRF scores output by prfcalc into a format with one SNP-gene score per line.
    usage:
       $0 --prf FILE [--genes FILE]
       \n";

my ($prfScoresFile, $genesFile);
GetOptions(	'prf=s' => \$prfScoresFile);

(!defined $prfScoresFile) and die("Missing parameter --prf.");

&checkFileExists($prfScoresFile) unless $prfScoresFile eq 'stdin';

*IN = &open_for_read2($prfScoresFile);
while(<IN>)
{
	chomp;
	my ($chr, $pos, $name, $XStr, $geneIDStr) = split("\t", $_, -1);
	my @geneIDs = split(',', $geneIDStr);
	my @Xs = split(',', $XStr);
	
	my $score = -20;
	my $geneID = "";
	if (scalar(@geneIDs)) {
		# Get the index of the max X value
		my $idxMax = 0;
		$Xs[$idxMax] >= $Xs[$_] or $idxMax = $_ for 1 .. $#Xs;
		$score = $Xs[$idxMax];
		$geneID = $geneIDs[$idxMax];
	}
	print join("\t", $chr, $pos, $name, $score, $geneID)."\n";
}
close(IN) unless $prfScoresFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;
