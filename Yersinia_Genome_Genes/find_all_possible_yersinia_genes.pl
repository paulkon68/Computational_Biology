# This script is designed to find all genes of the yersinia genome (yersinia_genome.fasta) and then save them in fasta format (yersinia_genes.fasta).

#! /usr/bin/perl

use strict;
use warnings;


# Check if the file exists
unless ( -e "yersinia_genome.fasta") {
	print "File \"yersinia_genome.fasta\" doesn\'t seem to exist!!\n";
	exit;
}

# Check if we can open the file correctly
unless ( open(FASTA, "<", "yersinia_genome.fasta") ) {
	print "Cannot open file \"yersinia_genome.fasta\"\n\n";
	exit;
}

# Check if we can open the file correctly
unless ( open(OUT, ">", "yersinia_genes.fasta") ) {
	print "Cannot open file \"yersinia_genes.fasta\"\n\n";
	exit;
}

my $fasta_seq = "";
my $fasta_header = <FASTA>; # The header of the fasta file: "yersinia_genome.fasta"
my @fasta_seq = <FASTA>; # the yersinia genome

close(FASTA); # close the input file (yersinia_genome.fasta)

$fasta_seq = join( '', @fasta_seq); # concatenate all the lines of the sequence of the genome into a string.
$fasta_seq =~ s/\s//g; # Remove whitespaces
my $gene_counter = 0; # Initialize the gene counter

# The regex is trying to match the 8-letter Shine-Dalgarno sequence ([TA][AC]AGGA[GA][GA]), 
# followed by 4-10 bases downstream before the initiation codon (ATG), including the initiation codon ATG.
# We should make the regex non-greedy in order to find the minimum number of four to ten bases before ATG codon.
# Otherwise, ATG codon could be found on those 4-10 bases thus giving us wrong gene matching.
while ($fasta_seq =~ /([TA][AC]AGGA[GA][GA][ATCG]{4,10}?)(ATG)/g) {
	
 	# Whenever a Shine-Dalgarno sequence is matched, store the position (index)
 	# of the first matched character in the sequence
 	# to remember where the potential gene begins.
 	my $start_matching_index = $-[1]; # index of the first matching character.
    my $start_searching_index = $-[2]; # index of the initiation codon (ATG).
   
 	

	my $gene = ""; # this variable will hold the different genes that we find.
	$gene .= $&; # assign to the gene variable the matching regex sequence (Shine-Dalgarno sequence to ATG initiation codon).

	# Start from the codon after the initiation codon (ATG) ($start_searching_index+3)
 	# Iterate until the end of the genome of yersinia (length($fasta_seq))
 	# Each time begin from the start of the next triplet ($pos+=3)
	for (my $pos = $start_searching_index+3; $pos < length($fasta_seq); $pos+=3) {
		
		# In the sequence, go to $pos (start of the current triplet)
		# and subset 3 characters starting from there.
		my $codon = substr($fasta_seq, $pos, 3);
		
		$gene .= $codon; # append each codon to the $gene variable.

		# If the current codon is a stop-codon, we reached the end of a potential gene!
		if ($codon eq "TAA" | $codon eq "TAG" | $codon eq "TGA") {			
			
			# The stop codon must NOT be found right after ATG, so only count this
			# as a gene.
			# NOTE: the first codon that is being assigned to $codon variable is the 
			# codon right after the ATG initiation codon.
			if ($pos != $start_searching_index+3) {

				# We need the position of last character of the stop codon
				# so add 2 to the current position -> start of the stop codon
				my $gene_end = $pos + 2;
				
				# Print the asked information
				# Remember to add 1 to the stored positions because Perl is zero-based!

				$gene_counter +=  1; # increment the gene counter

				my $gene_length = ($gene_end-$start_matching_index)+1; # compute the length of the gene.
				print OUT ">$gene_counter|+|", $start_matching_index+1, "|", $gene_end+1, "|$gene_length\n";
				print OUT "$gene\n";
				
			}
			last;
		}
	}
}

# The search in the reverse sequence is identical, 
# we just need to be be careful with the coordinates we print
# because we need to report them in respect to the forward strand

# Calculate the reverse of the fasta sequence
my $reverse_seq = reverse($fasta_seq);

# Substitute each base with its complementary one
$reverse_seq =~ tr/ACTG/TGAC/;

# The regex is trying to match the 8-letter Shine-Dalgarno sequence ([TA][AC]AGGA[GA][GA]), 
# followed by 4-10 bases downstream before the initiation codon (ATG), including the initiation codon ATG.
while ($reverse_seq =~ /([TA][AC]AGGA[GA][GA][ATCG]{4,10}?)(ATG)/g) {

	# Whenever a Shine-Dalgarno sequence is matched, store the position (index)
 	# of the first matched character in the sequence
 	# to remember where the potential gene begins.
	my $start_matching_index = $-[1]; # index of the first matching character.
	my $start_searching_index = $-[2]; # index of the initiation codon (ATG).

	my $gene = ""; # this variable will hold the different genes that we find.
	$gene .= $&; # assign to the gene variable the matching regex sequence (Shine-Dalgarno sequence to ATG initiation codon).

	# Start from the codon after the initiation codon (ATG) ($start_searching_index+3)
 	# Iterate until the end of the genome of yersinia (length($fasta_seq))
 	# Each time begin from the start of the next triplet ($pos+=3)
	for (my $pos = $start_searching_index+3; $pos < length($reverse_seq); $pos+=3) {		

		# In the sequence, go to $pos (start of the current triplet)
		# and subset 3 characters starting from there.
		my $codon = substr($reverse_seq, $pos, 3);

		$gene .= $codon; # append each codon to the $gene variable.

		# If the current codon is a stop-codon, we reached the end of the potential gene!
		if ($codon eq "TAA" | $codon eq "TAG" | $codon eq "TGA") {

			# The stop codon must NOT be found right after ATG, so only count this
			# as a gene.
			# NOTE: the first codon that is being assigned to $codon variable is the 
			# codon right after the ATG initiation codon.
			if ($pos != $start_searching_index+3) {
				
				# The start of the potential gene in reverse can be found if we subtract
				# the position of its start as calculated in reverse from the
				# total length of the sequence
				my $rev_start = length($reverse_seq) - $start_matching_index;

				my $rev_gene_end = length($reverse_seq) - ($pos+2); # same as above goes for the end of the potential gene.

				$gene_counter +=  1; # increment the gene counter

				my $gene_length = $rev_start-$rev_gene_end+1; # the length of the gene
				
				print OUT ">$gene_counter|-|", $rev_start, "|", $rev_gene_end, "|$gene_length\n";
				print OUT "$gene\n";

				# IMPORTANT NOTE FOR THE REVERSE STRAND
				# When we subtract, we don't need to add 1 anymore to the coordinates because
				# of the zero-based system of Perl
				# That is because when we calculate the length, we do it in one-based system
				# (e.g. the length of 0-2 is 3)
			}
			last;
		}
	}
}

close(OUT); # Close the output file (yersinia_genes.fasta)
