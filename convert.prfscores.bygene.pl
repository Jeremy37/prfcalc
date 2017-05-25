#!/usr/bin/perl -w

use strict;
use Getopt::Long;
sub open_for_read;
sub checkFileExists;

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
	*GENES_IN = &open_for_read($genesFile);
	while(<GENES_IN>)
	{
		chomp;
		$genesHash{$_} = 1;
	}
	close(IN);
	printf STDERR "$0: Read %d genes in file $genesFile\n", scalar(keys %genesHash);
}


*IN = &open_for_read($prfScoresFile);
my $header = <IN>;
#print $header;
while(<IN>)
{
	chomp;
	my ($chr, $pos, $name, $maxX, $maxGeneID, $XStr, $geneIDStr) = split("\t", $_, -1);
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


sub checkFileExists {
    die "$0: checkIsFile requires 1 arg - file" unless @_==1;
    my ($file) = @_;
    die "$0: checkIsFile: File $file does not exist" unless -e $file;
    die "$0: checkIsFile: File $file is not a regular file" unless -f $file;
}

sub open_for_read {
    die("open_for_read: Insufficient arguments to function open_for_read") unless @_ > 0;
    my $fname = shift;
    my $stat = shift if(@_);
    local *READ_FILE;
    if($fname eq 'stdin') {
	*READ_FILE = *STDIN;
    } else {
	if($fname =~ /.gz$/) {
	    die "Cannot open file $fname for read. Reason :$!.\nExiting...\n" unless -e $fname;
	    open(READ_FILE,"gunzip -c $fname |") || die "Cannot open file $fname for read. Reason :$!.\nExiting...\n";
	} else {
	    open(READ_FILE,"$fname") || die "Cannot open file $fname for read. Reason :$!.\nExiting...\n";
	}
	print STDERR "$0: Open file $fname for read\n" unless $stat && $stat eq 'q';
    }
    return *READ_FILE;
}

