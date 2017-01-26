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
GetOptions(	'prf=s' => \$prfScoresFile,
			'genes=s' => \$genesFile);

(!defined $prfScoresFile) and die("Missing parameter --prf.");

&checkFileExists($prfScoresFile) unless $prfScoresFile eq 'stdin';
my %genesHash = ();
if ($genesFile) {
	&checkFileExists($genesFile);
	*GENES_IN = &open_for_read2($genesFile);
	while(<GENES_IN>)
	{
		chomp;
		$genesHash{$_} = 1;
	}
	close(IN);
	printf STDERR "$0: Read %d genes in file $genesFile\n", scalar(keys %genesHash);
}


*IN = &open_for_read2($prfScoresFile);
my $header = <IN>;
#print $header;
while(<IN>)
{
	chomp;
	my ($chr, $pos, $name, $XStr, $geneIDStr) = split("\t", $_, -1);
	my @geneIDs = split(',', $geneIDStr);
	my @Xs = split(',', $XStr);
	my $i = 0;
	for my $geneID (@geneIDs) {
		if (!$genesFile or exists $genesHash{$geneID}) {
			print join("\t", $chr, $pos, $name, $Xs[$i], $geneID)."\n";
		}
		$i += 1;
	}
}
close(IN) unless $prfScoresFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;
