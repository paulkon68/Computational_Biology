#!/usr/bin/perl
use warnings;
use strict;
use Scalar::Util qw(reftype);

my $seqs=$ARGV[0];  # take the filename as the first argument from the command line.
my $output=$ARGV[1];
my %PFM;  # create a hash for the PFM
my $base;
my $length;



open(INPUT,"<",$seqs) or die "$!\n";  # open the file with filename contained in $seq for reading
open(OUTPUT,">","$output") or die "$!\n";  # open the file "PFM.txt" for writing

# we check if the sequences contain valid letters.
while(my $line=<INPUT>){ 
	chomp $line;  # remove the '\n' character
	if($line!~/^[ATCG]+$/) # check if all lines contain A, C, T, or G
	{ 	
		print "No valid letters in $seqs \n";
		exit;
	}
}
close(INPUT);

open(INPUT,"<",$seqs) or die "$!\n";

# Calculate each nucleotide's frequency in each position
while(my $line=<INPUT>){ 

	chomp $line;
	$length=length($line);

	for (my $i=0; $i<length($line); $i++)
	{ 	
		$base=substr($line,$i,1);
		$PFM{$base}{$i+1}++;
	}
}

# Print the Position Frequency Matrix

foreach my $base (keys %PFM){
	print OUTPUT $base . "\t";  # print the key of the initial hash, which is the base
	for (my $i = 1; $i <= $length; $i++){
		if(!defined($PFM{$base}{$i})){
			$PFM{$base}{$i} = 0; # initialize the secondary hash values if they are not defined, which are the 
								 # frequencies of the base in the different positions in the PFM
		}
		print OUTPUT $PFM{$base}{$i} . "\t"; # print the lines of PFM
	}
	print OUTPUT "\n";
}

# Find the consensus sequence and calculate the score

my @consensus;
my @consensus2;
my $max = 0;
my $score = 0;
print "\n";

for(my $i = 1; $i <= $length; $i++){

	$max = 0;
	# print $PFM{'A'}{$i}, "\t", $PFM{'T'}{$i}, "\t", $PFM{'C'}{$i}, "\t", $PFM{'G'}{$i}, "\n";
	# Find for each position the most frequent nucleotide and save its value to "$max" variable
	if ($PFM{'A'}{$i} >= $max){
		$max = $PFM{'A'}{$i};
		$consensus[$i-1] = 'A'; # save the most frequent nucleotide found, in the respective position in the array "$consensus"
	}
	if ($PFM{'T'}{$i} >= $max){
		$max = $PFM{'T'}{$i};
		$consensus[$i-1] = 'T';
	}
	if ($PFM{'C'}{$i} >= $max){
		$max = $PFM{'C'}{$i};
		$consensus[$i-1] = 'C';
	}
	if ($PFM{'G'}{$i} >= $max){
		$max = $PFM{'G'}{$i};
		$consensus[$i-1] = 'G';
	}

	# Check if two or more nucleotides have the same frequency in each position
	# If that's true, save the nucleotides in the multidimensional array "@consensus2"
	if ($PFM{'A'}{$i} == $max){
		push(@{$consensus2[$i-1]}, 'A');
	}
	if ($PFM{'T'}{$i} == $max){
		push(@{$consensus2[$i-1]}, 'T');
	}
	if ($PFM{'C'}{$i} == $max){
		push(@{$consensus2[$i-1]}, 'C');
	}
	if ($PFM{'G'}{$i} == $max){
		push(@{$consensus2[$i-1]}, 'G');
	}

	$score+=$max;  # score of the consensus sequence
}

# Print the consensus sequence to the file
print OUTPUT @consensus, "\n";

# print the consensus sequence, but with all of the possible bases for each position
for (my $i = 0; $i < $length; $i++){

	print OUTPUT ${$consensus2[$i]}[0];  

	if (defined ${$consensus2[$i]}[1]){
		print OUTPUT "|" . ${$consensus2[$i]}[1];
	}
	if (defined ${$consensus2[$i]}[2]){
		print OUTPUT "|" . ${$consensus2[$i]}[2];
	}
	if (defined ${$consensus2[$i]}[3]){
		print OUTPUT "|" . ${$consensus2[$i]}[3];
	}
	
}

# Print the Score of the consensus sequence to the output file
print OUTPUT "\nScore= ", $score, "\n";

# Close the two file handlers
close(INPUT);
close(OUTPUT);