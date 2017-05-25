#!/usr/bin/perl -w

use strict;
use Getopt::Long;
sub checkFileExists;
sub open_for_read;

my $usage = "$0 - Code to take TSS info from FANTOM if available, or Ensembl otherwise.
    usage:
       $0 --fantom FILE --ensembl FILE
    options:
       -XXX\n";

## Read options from CL
my ($fantomFile, $ensemblFile);
GetOptions(	'fantom=s' => \$fantomFile,
			'ensembl=s' => \$ensemblFile);

## Check CL args
&checkFileExists($fantomFile);
&checkFileExists($ensemblFile);

my %fantomGenesSeen;
*FANTOM = &open_for_read($fantomFile);
while(<FANTOM>) {
	chomp();
	my @l = split;
	$fantomGenesSeen{$l[2]} = 1;
	if ($l[0] =~ /chr(.+)/) {
		$l[0] = $1;
	}
	print join("\t", @l)."\n";
}
close(FANTOM);

*ENSEMBL = &open_for_read($ensemblFile);
while(<ENSEMBL>) {
	my $line = $_;
	chomp();
	my @l = split;
	
	(!exists $fantomGenesSeen{$l[2]}) and print $line;
}
close(ENSEMBL);


sub usage {
    print STDERR $usage;
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

