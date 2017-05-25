#!/usr/bin/perl -w

use strict;
use Getopt::Long;
sub checkFileExists;
sub open_for_read;

sub usage($);
my $usage = "$0 - Do a hash join on two files
    usage:
       $0 --hashFile FILE --scanFile FILE --colHashFile INT --colScanFile INT [--sep '\\t']
    options:
       --invert=XXX - print lines from --scanFile that do NOT match anything in hashFile
       \n";

## Read options from CL
my $colHashFile = 1;
my $colScanFile = 1;

my ($hashFile, $scanFile, $outerjoin, $hasHeader, $invert, $caseinsensitive, $sep);
GetOptions(	'hashFile=s' => \$hashFile,
			'scanFile=s' => \$scanFile,
			'colHashFile=i' => \$colHashFile,
			'colScanFile=i' => \$colScanFile,
			'outerjoin' => \$outerjoin,
			'header' => \$hasHeader,
			'caseinsensitive' => \$caseinsensitive,
			'invert|v' => \$invert,
			'sep=s' => \$sep);

## Check CL args
(!defined $hashFile) and &usage("Missing parameter --hashFile.");
(!defined $colHashFile or !defined $colScanFile) and &usage("Missing one of the parameters colHashFile, colScanFile.");
my $inSep = '\t';
my $outSep = "\t";
(defined $sep) and $inSep = $sep;
(defined $sep) and $outSep = $sep;
($outerjoin and $invert) and &usage("Parameters --outerjoin and --invert are incompatible.");

($colHashFile < 1 or $colScanFile < 1) and &usage("Error: column parameters should be 1-based.");
$colHashFile = $colHashFile - 1;
$colScanFile = $colScanFile - 1;

&checkFileExists($hashFile) unless $hashFile eq 'stdin';
&checkFileExists($scanFile);
my %hashFile1;

## Read input files
*IN = &open_for_read($hashFile);
my $headerHashFile;
if ($hasHeader) {
	$headerHashFile = <IN>;
	chomp($headerHashFile);
}
while(<IN>)
{
	chomp;
	my @l = split($inSep, $_, -1);
	die "$0: Column $colHashFile does not exist in file $hashFile, line $." unless $#l >= $colHashFile;
	($colHashFile > $#l) and next; # Skip lines with blank key value
	my $key = $l[$colHashFile];
	($caseinsensitive) and $key = uc($l[$colHashFile]);
	# Store the first line that has this key
	if (not exists $hashFile1{$key}) {
		$hashFile1{$key} = $_;
	}
}
close(IN) unless $hashFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;
printf STDERR "$0: Read %d records in file $hashFile\n", scalar(keys %hashFile1);


*IN = &open_for_read($scanFile);
if ($hasHeader) {
	my $headerScanFile = <IN>;
	chomp($headerScanFile);
	my @lineHashFile = split($inSep, $headerHashFile, -1);
	my @lineScanFile = split($inSep, $headerScanFile, -1);
	splice @lineHashFile, $colHashFile, 1;
	print join($outSep, @lineScanFile, @lineHashFile)."\n";
}
while(<IN>)
{
	chomp;
	($_ =~ /^#/) and next;
	my @l = split($inSep, $_, -1);
	($colScanFile > $#l) and print "$_\n" and next; # Just print lines with blank key value

	my $key = $l[$colScanFile];
	($caseinsensitive) and $key = uc($l[$colScanFile]);
	if (not $invert) {
		if (exists $hashFile1{$key}) {
			my @lineFile1 = split($inSep, $hashFile1{$key}, -1);
			splice @lineFile1, $colHashFile, 1;
			print join($outSep, @l, @lineFile1)."\n";
		} elsif ($outerjoin) {
			print "$_\n";
		}
	} else {
		if (not exists $hashFile1{$key}) {
			print "$_\n";
		}
	}
}
close(IN) unless $hashFile;



## Functions 

sub usage($) {
	print $_[0]."\n";
	print STDERR $usage;
	exit;
}

sub errAbort {
    ## Optional args:
    ## 1: Error message
    ## 2: 'q' suppresses usage message
    my $arg = shift;
    if(defined $arg) {
	print STDERR "Error: $arg\n";
    }
    undef $arg;
    $arg = shift(@_);
    if(defined $arg) {
	if($arg ne 'q') {
	    &usage();
	}
    } else {
	&usage();	
    }
    exit(1);
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

sub checkFileExists {
    die "$0: checkIsFile requires 1 arg - file" unless @_==1;
    my ($file) = @_;
    die "$0: checkIsFile: File $file does not exist" unless -e $file;
    die "$0: checkIsFile: File $file is not a regular file" unless -f $file;
}
