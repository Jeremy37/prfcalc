#!/usr/bin/perl -w

use strict;
use Getopt::Long;
sub open_for_read;
sub checkFileExists;
sub usage($);
my $usage = "$0 - Do a hash join on two files
    usage:
       $0 --hashFile FILE --scanFile FILE --matchColHash INT --matchColScan INT --annotateColHash INT --annotateColScan INT --annotationSep STR [--sep '\\t']
    options:
       --verbose=XXX - set verbosity level
       \n";

## Read options from CL

my ($hashFile, $scanFile, $matchColHash, $matchColScan, $annotateColHash, $annotateColScan, $annotationSep, $sep, $verbosity);
$verbosity = 0;

GetOptions(	"verbose|v=i" => \$verbosity,
			'hashFile=s' => \$hashFile,
			'scanFile=s' => \$scanFile,
			'matchColHash=s' => \$matchColHash,
			'matchColScan=s' => \$matchColScan,
			'annotateColHash=s' => \$annotateColHash,
			'annotateColScan=s' => \$annotateColScan,
			'annotationSep=s' => \$annotationSep,
			'sep=s' => \$sep);

## Check CL args
(!defined $hashFile) and &usage("Missing parameter --hashFile.");
(!defined $matchColHash or !defined $matchColScan) and &usage("Missing one of the parameters matchColHash, matchColScan.");
(!defined $annotateColScan) and &usage("Missing parameter annotateColScan.");
(defined $sep) or $sep = "\t";

my @colsHashFile = split(/,/, $matchColHash);
my $maxMatchColHash = 0;
for (my $i = 0; $i <= $#colsHashFile; $i++) {
	$colsHashFile[$i] -= 1;
	($colsHashFile[$i] < 0) and &usage("Error: column parameters should be 1-based.");
	($colsHashFile[$i] > $maxMatchColHash) and $maxMatchColHash = $colsHashFile[$i];
}
my @colsScanFile = split(/,/, $matchColScan);
my $maxMatchColScan = 0;
for (my $i = 0; $i <= $#colsScanFile; $i++) {
	$colsScanFile[$i] -= 1;
	($colsScanFile[$i] < 0) and &usage("Error: column parameters should be 1-based.");
	($colsScanFile[$i] > $maxMatchColScan) and $maxMatchColScan = $colsScanFile[$i];
}

my @annotateColsScanFile = split(/,/, $annotateColScan);
my $maxAnnColScan = 0;
for (my $i = 0; $i <= $#annotateColsScanFile; $i++) {
	$annotateColsScanFile[$i] -= 1;
	($annotateColsScanFile[$i] < 0) and &usage("Error: column parameters should be 1-based.");
	($annotateColsScanFile[$i] > $maxAnnColScan) and $maxAnnColScan = $annotateColsScanFile[$i];
}

my @annotateColsHashFile;
if (defined $annotateColHash) {
	@annotateColsHashFile = split(/,/, $annotateColHash);

	if (scalar(@annotateColsHashFile) != scalar(@annotateColsScanFile)) {
		die "$0: if any columns are specified with --annotateColHash, the number of them must match the number specified with --annotateColScan.";
	}

	my $maxAnnColHashFile = 0;
	for (my $i = 0; $i <= $#annotateColsHashFile; $i++) {
		$annotateColsHashFile[$i] -= 1;
		($annotateColsHashFile[$i] < 0) and &usage("Error: column parameters should be 1-based.");
		($annotateColsHashFile[$i] > $maxAnnColHashFile) and $maxAnnColHashFile = $annotateColsHashFile[$i];
	}
}
if ($verbosity > 0) {
	print STDERR "$0: colsHashFile: ".join(",", @colsHashFile)."\n";
	print STDERR "$0: colsScanFile: ".join(",", @colsScanFile)."\n";
	(defined $annotateColHash) and print STDERR "$0: annotateColsHashFile: ".join(",", @annotateColsHashFile)."\n";
	print STDERR "$0: annotateColsScanFile: ".join(",", @annotateColsScanFile)."\n";
}

&checkFileExists($hashFile) unless $hashFile eq 'stdin';
&checkFileExists($scanFile);
my %hashFile1;

## Read input files
my $numHashFileCols = 0;
*IN = &open_for_read($hashFile);
my @hashFileLines = <IN>;
for (my $i = 0; $i <= $#hashFileLines; $i++)
{
	my $line = $hashFileLines[$i];
	($line =~ /^#/ or $line =~ /^\s*$/) and next;
	chomp($line);
    my @l = split($sep, $line, -1);
    die "$0: Column $maxMatchColHash does not exist in file $hashFile, line $.\n$line" unless @l >= $maxMatchColHash;
    my $numFields = scalar(@l);
    ($numFields < $numHashFileCols) and print STDERR "$0: WARNING:Found line with fewer columns ($numFields) than some other lines ($numHashFileCols)!\n";
    ($numFields > $numHashFileCols) and $numHashFileCols = $numFields;
    
	my $key = join('\t', @l[@colsHashFile]);
    $hashFile1{$key} = undef;
}
close(IN) unless $hashFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

if (!defined $annotateColHash) {
	# The columns to replace (annotateColsHashFile) were not specified, so the annotation
	# columns are added after all current columns.
	for (my $i = 0; $i <= $#annotateColsScanFile; $i++) {
		$annotateColsHashFile[$i] = $numHashFileCols + $i;
	}
	($verbosity > 0) and print STDERR "$0: annotateColsHashFile: ".join(",", @annotateColsHashFile)."\n";
}

printf STDERR "$0: Read %d unique records in file $hashFile\n", scalar(keys %hashFile1);
my $numMatches = 0;
my $matchSep = '@#';
my $valSep = ',;';
*IN = &open_for_read($scanFile);
while(<IN>)
{
	($_ =~ /^#/) and next;
	chomp();
	my @l = split($sep, $_, -1);
	my $key = join('\t', @l[@colsScanFile]);
	my $allempty = ($key =~ /^\s*$/);
	if (not $allempty and exists $hashFile1{$key}) {
		$numMatches = $numMatches + 1;
		
		my $newMatchVals = join($valSep, @l[@annotateColsScanFile]);
		if (defined $annotationSep and defined $hashFile1{$key}) {
			$hashFile1{$key} = join($matchSep, $hashFile1{$key}, $newMatchVals);
		} else {
			$hashFile1{$key} = $newMatchVals;
		}
	}
}

($verbosity > 0) and print STDERR "$0: Matches to hash file lines: $numMatches\n";
($verbosity > 0) and print STDERR "$0: Outputting data\n";

my $numDataLines = 0;
my $numAnnotatedLines = 0;
my $numAnnotationsWritten = 0;
for (my $i = 0; $i <= $#hashFileLines; $i++)
{
	($hashFileLines[$i] =~ /^#/) and $hashFileLines[$i] and next;
	$numDataLines++;
	chomp($hashFileLines[$i]);
    my @l = split($sep, $hashFileLines[$i], -1);
	my $key = join('\t', @l[@colsHashFile]);
	(!exists $hashFile1{$key}) and die "ERROR: All lines from input hash file should be present in hash table!\n";
	
	my @annotationsToWrite;
	if (defined $hashFile1{$key}) {
		$numAnnotatedLines++;
		# $hashFile1{$key} has data of the form "col1valMatch1,col2valMatch1,col3valMatch1@#col1valMatch2,col2valMatch2,col3valMatch2".
		# We want to rearrange it to get an array with "col1valMatch1,col1valMatch2", "col2valMatch1,col2valMatch2", etc.
		my @scanColMatches = split($matchSep, $hashFile1{$key}, -1);
		# Need to collect matches for a particular column together.
		for (my $j = 0; $j <= $#scanColMatches; $j++) {
			my @colVals = split($valSep, $scanColMatches[$j], -1);
			for (my $k = 0; $k <= $#colVals; $k++) {
				if (defined $annotationSep and defined $annotationsToWrite[$k]) {
					$annotationsToWrite[$k] .= $annotationSep.$colVals[$k];
				} else {
					if (!defined $colVals[$k]) {
						print STDERR "WTF?? k=$k, j=$j, hashFileVal=$hashFile1{$key}\n";
						exit(1);
					}
					$annotationsToWrite[$k] = $colVals[$k];
				}
			}
		}
	}
	
	# Copy in values from the scan file into the hash file line
	for (my $k = 0; $k <= $#annotateColsHashFile; $k++) {
		my $col = $annotateColsHashFile[$k];
		($col > scalar(@l)) and die "$0: ERROR:Hash file column to annotate ($col) is beyond end of the line (length ".scalar(@l).") \n";
		if (defined $annotationsToWrite[$k]) {
			$numAnnotationsWritten++;
			$l[$col] = $annotationsToWrite[$k];
		} elsif (!defined $l[$col]) {
			# Neither the annotation/scan file nor the hash file has a value at
			# this column. Write an empty value so that the line output has the
			# correct number of fields.
			$l[$col] = "";
		}
	}
    print join($sep, @l)."\n";
}
close(IN) unless $hashFile;

($verbosity > 0) and print STDERR "$0: Annotated $numAnnotatedLines of $numDataLines lines\n";
($verbosity > 0) and print STDERR "$0: Added $numAnnotationsWritten annotations\n";



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

