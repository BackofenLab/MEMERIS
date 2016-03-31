#!/usr/bin/perl

# call: perl -w ReplaceUbyT.perl inputfile outputfile\n";
# inputfile has to have the following format:
#		>name
#		sequence
#		>name2
#		sequence2
#	   ...


use strict;
use warnings;
use Cwd;


if (@ARGV < 2) {
	print "usage: perl -w ReplaceUbyT.perl inputfile outputfile \n";
	print "this perl file replaces all U's by T's in the input sequences\n";
	exit(0);
}


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];


# open input and output file
open(INPUT, $inputfile);
open(OUTPUT, ">$outputfile");

# read all sequences
while (my $line1 = <INPUT>) {
	chomp($line1); 					# fasta header
	my $line2 = <INPUT>; 			# sequence
	chomp($line2); 

	my $seq = uc $line2;				# sequence
	$seq =~ tr/U/T/;					# T instead of U
	
	print OUTPUT "$line1\n$seq\n";			# output

}	

close(INPUT);	
close(OUTPUT);




