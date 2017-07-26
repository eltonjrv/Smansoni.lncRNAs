#!/usr/bin/perl
# Programmer: Elton Vasconcelos (24/Fev/2015)
# Script that reports masked sequences by an user-defined coverage percentage threshold
# Usage: perl RM-cov_cutoff.pl [RM_file.out.gff|.cat] [seq.fasta] [fraction of coverage, fraction from 0.x to 1] >outfile.gff
####################################################################
# ATTENTION-1: Use the option "-gff" on your RepeatMasker cmd line to generate a gff output
# If you didn't use "-gff" on your RepeatMasker execution, use the *.cat file as input for this script instead
# ATTENTION-2: the headline of the fasta file must end with the regexp '\d+ bp$'
# run readseq.jar to generate this kind of fasta file, but make sure the name of the sequences haven't changed
# ATTENTION-3: The cat file must be sorted by field 5 (gene or contig ID)
############################################################################
######## Working with the fasta file to make the hash table ################
############################################################################
open (SEQFILE, "$ARGV[1]") or die ("Can't open $ARGV[1]\n");
my $seq_line = <SEQFILE>;
chomp($seq_line);
my $seqID;
my $seq_size;
my %hash;

while ($seq_line ne "") {
    if ($seq_line =~ m/^>.*/) {
        $seqID = $&;                #catching only the seqID ID
        $seqID =~  s/^>//g;
        $seqID =~  s/ .*//g;        #Keeping only the first word of the sequence (after ">" and before the first blank space)
        $seq_line =~ m/\d+ bp$/;
        $seq_size = $&;	            #armazenando o tamanho (pb) da sbjct
        $seq_size =~ s/\,//g;         #removendo a virgula
        $seq_size =~ s/ bp//g;          #armazenando apenas o tamanho (pb) da sbjct
        $hash{$seqID} = $seq_size;      #Creating a hash containing each seqID size. It will be used on line 73
        $seq_line = <SEQFILE>;
        chomp($seq_line);
    }
    else {        #Just in the case the FASTA file is not starting with a headline (^>.+)
        $seq_line = <SEQFILE>;
        chomp($seq_line);
    }
}
#my @test = keys %hash;
######################################################
if ($ARGV[0] =~ m/\.gff$/) {    #input file is in gff format

    open (INFILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
    my $line = <INFILE>;
    chomp($line);
    my @array;
    my $masked_region;              # region on the sequence that was masked by RepeatMasker (cols 4 and 5 from the gff)
    my $cov_frac = $ARGV[2];	    # coverage fraction cutoff
    while ($line ne "") {
        if ($line =~ m/^\#/) {
            $line = <INFILE>;
            chomp($line);

        }

        else {
            @array = split(/\t/, $line);
            $masked_region = ($array[4] - $array[3]) + 1;

            if ($masked_region >= $hash{$array[0]} * $cov_frac) {
                print("$line\n");

            }
            $line = <INFILE>;
            chomp($line);
        }
    }
}

elsif ($ARGV[0] =~ m/\.cat$/) {    #input file is in cat format
    open (INFILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
    my $line = <INFILE>;
    chomp($line);
    my @array;
    my $masked_region;              # region on the sequence that was masked by RepeatMasker (cols 6 and 7 from the cat)
    #my $masked_region_blocks = 0;   # sum of all the masked blocks in the same sequence

    my $cov_frac = $ARGV[2];	    # coverage fraction cutoff
    my $geneID;
    my $percent_mask;
    while ($line ne "") {
        if ($line !~ m/^\d/) {	    # Lines not starting with digits are not results
            $line = <INFILE>;
            chomp($line);

        }

        else {
            @array = split(/ /, $line);    #cat file is not a tabular one. The fields are separated by blank space
            $geneID = $array[4];
            my $coord1 = $array[5];
            my $coord2 = $array[6];
	    $masked_region = ($coord2 - $coord1) + 1; 
            while ($array[4] eq $geneID) {
                $line = <INFILE>;
	        chomp($line);
        	@array = split(/ /, $line);  
		if ($array[4] eq $geneID && $array[5] <= $coord2 && $array[6] > $coord2) {
                	$masked_region = $masked_region + ($array[6] - $coord2);
			$coord1 = $array[5];
			$coord2 = $array[6];
			#print("$geneID\t$masked_region\n");
 		}
		elsif ($array[4] eq $geneID && $array[5] >= $coord2) {
        	        $masked_region = $masked_region + ($array[6] - $array[5] + 1);
			$coord1 = $array[5];
			$coord2 = $array[6];
			#print("$geneID\t$masked_region\n")
     	       	}
	    }
            if ($masked_region >= $hash{$geneID} * $cov_frac) {
                $percent_mask = ($masked_region / $hash{$geneID}) * 100;
                print("$geneID\t$percent_mask\n");

            }
        }
    }
}

else {
    print STDERR ("*** Error: The input file must be in either .gff or .cat format ***\n");
}
    

    
