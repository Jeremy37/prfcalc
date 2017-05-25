#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Cwd qw(abs_path);
use File::Basename;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Reads a list of gene symbols and uses the HGNC file to convert any old symbols to the latest.
    usage:
       $0 --file FILE --hgncfile FILE
	   \n";

## Read options from CL
my ($infile, $hgncfile);
GetOptions(	"verbose=i" => \$Common::verbosity,
			'file=s' => \$infile,
			'hgncfile=s' => \$hgncfile);

my $scriptPath = dirname(abs_path($0));
!$hgncfile and $hgncfile = $scriptPath."/gene_with_protein_product.txt";

## Check CL args
&checkFileExists($infile) unless $infile eq 'stdin';
&checkFileExists($hgncfile);

my %approvedSymbols;
my %previousSymbols;

# First read the whole hgnc file into a hash
*HGNC = &open_for_read2($hgncfile);
while(<HGNC>) {
	my $line = $_;
	chomp();
	my @l = split(/\t/, $_, -1);
	my $approvedSym = $l[1];
	if (!$approvedSym) { next; }
	$approvedSymbols{uc($approvedSym)} = 1;

	# Add each previous symbol or synonym in a hash with the value
	# being the approved symbol
	my @prevSymbols = split(', ', $l[4]);
	my @symonyms = split(', ', $l[6]);
	foreach my $sym (@prevSymbols) {
		$previousSymbols{uc($sym)} = $approvedSym;
	}
	foreach my $sym (@symonyms) {
		$previousSymbols{uc($sym)} = $approvedSym;
	}
}
close(HGNC);


*IN = &open_for_read2($infile);
while(<IN>) {
	chomp();
	my $sym = uc($_);
	if (!exists $approvedSymbols{$sym} and exists $previousSymbols{$sym}) {
		print $previousSymbols{$sym}."\n";
	} else {
		print "$_\n";
	}
}
close(IN) unless $infile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;


sub usage {
    print STDERR $usage;
}
