#! /usr/bin/perl


use strict;
use warnings;


# open the "human_exons_prCoding_exercise_set.bed" file as input file.
open(BEDFILE, "<", "human_exons_prCoding_exercise_set.bed") or die "$!\n"; 

# open the "first_exons_coordinates.bed" file as the first output file.
open(WRITEFILE_0, ">", "first_exons_coordinates.bed") or die "$!\n";

# open the "exonic_length_per_transcript.txt" file as the second output file.
open(WRITEFILE_1, ">", "exonic_length_per_transcript.txt") or die "$!\n";


# Create a hash to store the start and the end of each transcript, belonging to a specific transcript id.
my %lengths;

# read the first line of the bed file.
my $bedline = <BEDFILE>; 

# Remove trailing \n (or \r\n for Windows) at the end of the line
chomp $bedline; 

# Split the line in tabs (\t) and store the values in an array
my @first_entry = split("\t", $bedline);

# pick the first element of the array '@first_entry' that contains the chromosome name.
my $chr_number = $first_entry[0];

# pick the second element of the array '@first_entry' that contains the start coordinate of the exon in the chromosome.
my $chrom_Start = $first_entry[1];

# pick the third element of the array '@first_entry' that contains the end coordinate of the exon in the chromosome.
my $chrom_End = $first_entry[2];

# split the 4th column of the file (3 because Perl is 0-based, 0,1,2,3...) which stores the Ensembl IDs in @ 
# and save the values in an array.
my @ensembl_ids = split("@", $first_entry[3]);

# pick the first element of the array '@ensembl_ids' (1 due to 0-based), which is the Ensembl Gene ID.
my $gene_id = $ensembl_ids[0];

# pick the second element of the array '@ensembl_ids', which is the transcript ID.
my $transcript_id = $ensembl_ids[1];

# pick the fifth element of the array '@first_entry', which is the score.
my $score = $first_entry[4];

# pick the sixth element of the array '@first_entry', which is the strand.
my $strand = $first_entry[5];

# Because 1 gene can have more than 1 transcript,
# the value of the hash for each transcript ID (key)
# will be an array with all its start and end coordinates, in the following pattern:
# start_coordinate,end_coordinate,start_coordinate,end_coordinate,... .
	
# Thus we push this transcript_id into the 'lengths' hash at the key 'transcript_id'
# But in order to push, we tell Perl that there is an array at that slot of the hash
# by dereferencing that position as an array! -> @{hash{key}}
# Internally, Perl also creates the array if it doesn't exist 
push @{$lengths{$transcript_id}}, $chrom_Start;
push @{$lengths{$transcript_id}}, $chrom_End;

# Write in the output file WRITEFILE_0 as asked
print WRITEFILE_0 $chr_number, "\t", @{$lengths{$transcript_id}}[0], "\t", @{$lengths{$transcript_id}}[1], "\t", $first_entry[3], "\t", $score, "\t", $strand, "\n";

# For each line in BEDFILE (while it has lines...)
while ($bedline = <BEDFILE>) {

	# Remove trailing \n (or \r\n for Windows) at the end of the line
	chomp $bedline;
	
	# Split the line in tabs (\t) and store the values in an array
	my @columns = split("\t", $bedline);

	# pick the first element of the array that contains the chromosome name.
    $chr_number = $columns[0];

	# pick the second element of the array that contains the start coordinate of the exon in the chromosome.
    $chrom_Start = $columns[1];

	# pick the third element of the array that contains the end coordinate of the exon in the chromosome.
    $chrom_End = $columns[2];

	# split the 4th column of the file (3 because Perl is 0-based, 0,1,2,3...) which stores the Ensembl IDs in @ 
	# and save the values in an array.
    @ensembl_ids = split("@", $columns[3]);
	
	# pick the second element of the array '@ensembl_ids', which is the transcript ID.
    $transcript_id = $ensembl_ids[1];

	# pick the fifth element of the array , which is the score.
    $score = $columns[4];

	# pick the sixth element of the array, which is the strand.
    $strand = $columns[5];

	# if the Ensembl Gene ID is equal to the previous Ensembl Gene ID, that means that there are still exons with the same 
	# gene id. Save all the start and end coordinates as an array to the 'length' hash.
	# We will use this hash to compute the total length of each exonic region.
    if ($ensembl_ids[0] eq $gene_id){
        push @{$lengths{$transcript_id}}, $chrom_Start;
        push @{$lengths{$transcript_id}}, $chrom_End;
    }
	# Alternatively, (means that we found an exon with different gene id than the previous)
    else{

		# change the 'gene_id' variable to this new first exon of the new transcript
        $gene_id = $ensembl_ids[0];

		# save the start and end coordinates to the 'lengths' hash with key: '$transcript_id'
        push @{$lengths{$transcript_id}}, $chrom_Start;
        push @{$lengths{$transcript_id}}, $chrom_End;

		# print to the 'first_exons_coordinates.bed' file as asked. The first exon of each transcipt will be saved.
        print WRITEFILE_0 $chr_number, "\t", $chrom_Start, "\t", $chrom_End, "\t", $columns[3], "\t", $score, "\t", $strand, "\n";
    }
	
	
}
# For each key (i.e. transcript ID) in the 'lengths' hash
foreach my $key (keys %lengths) {

	# this variable will hold the length of each exon for a specific transcript.
    my $length_of_exonic_region = 0;

	# save the $key which is the transcript_id in the 'exonic_length_per_transcript.txt' file
    print WRITEFILE_1 $key;

	# compute the lengths of the exons of each of the transcripts.
	for (my $i=0; $i < scalar(@{$lengths{$key}})-1; $i+=2){
        $length_of_exonic_region += (@{$lengths{$key}}[$i+1] - @{$lengths{$key}}[$i]);
    }

	# print the lengths to the 'exonic_length_per_transcript.txt' file.
    print WRITEFILE_1 "\t", $length_of_exonic_region, "\n";
}

# Close the input and output file handlers.
close(BEDFILE);
close(WRITEFILE_0);
close(WRITEFILE_1);
