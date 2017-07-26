#!/usr/bin/perl
# Programmer: Elton Vasconcelos (10/08/2015)
# Script that searches for canonical splice donor (GT) and acceptor (AG) within introns' boundaries
# Usage: perl catch_canonical_spliceSites.pl [infile.bed] [whole_genome.fasta] >outfile.bed
############################################################
# NOTE-1: The chromosomes' IDs must be the same on both .bed and .fasta input files.

########################################################
######## Working with the genome.fasta file ############
########################################################
open (FILE2, "$ARGV[1]") or die ("Can't open file $ARGV[1]\n");
my $line2 = <FILE2>;
chomp($line2);
my $seq;
my $chr;
my %hash;

while ($line2 ne "") {
    if ($line2 =~ m/^>[\w\d\._]+ /) {
        $chr = $&;                #catching only the chr ID (i.e. LmjF.01)
        $chr =~  s/^>//;
        $chr =~  s/ //;
        $line2 = <FILE2>;
        chomp($line2);
        until ($line2 =~ m/^>/ || $line2 eq "") {             
            $seq .= $line2;            #putting the whole chromosome sequence in one single line
            $line2 = <FILE2>;
            chomp($line2);
        }
        $hash{$chr} = $seq;            #Creating a hash containing each chromosome in the file
        $seq = "";
    }
    else {        #Just in the case the FASTA file is not starting with a headline (^>.+)
        $line2 = <FILE2>;
        chomp($line2);
    }
}

##############################################################################
########################## Working on the bed file ###########################
##############################################################################

open (INFILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my $line = <INFILE>;
chomp($line);
my (@array, @col11, @col12, $m, $splice_donor, $splice_acceptor, $intron_coord1, $intron_coord2);
#my $cutoff = $ARGV[1];    #$ARGV[1] must be an integer (>= 1). That's the minimum intron length

while ($line ne "") {
    $m = 0;
    @array = split(/\t/, $line);
    $array[10] =~ s/\,$//g;    #only in case these columns are ending with ","
    $array[11] =~ s/\,$//g;
    @col11 = split(/\,/, $array[10]);
    @col12 = split(/\,/, $array[11]);
    if ($line =~ m/\t\+\t/) {           ### working on plus strand
        for ($i = 0; $i < $array[9] - 1; $i++) {         #$array[9] is the number of exons. So, $array[9] - 1 is the number of introns
            $intron_coord1 = $array[1] + $col12[$i] + $col11[$i];
            $intron_coord2 = $array[1] + $col12[$i+1] - 1;
            $splice_donor = substr($hash{$array[0]}, $intron_coord1, 2);         #catching 2 nt from the intron start
            $splice_acceptor = substr($hash{$array[0]}, $intron_coord2 - 1, 2);  #catching the last 2 nt of the intron
            print STDERR ("$array[0]\t$array[3]\t$array[5]\t$splice_donor\t$splice_acceptor\n");
            if ($splice_donor eq "GT" && $splice_acceptor eq "AG") {
                $m++; 
            }
	    elsif ($splice_donor eq "gt" &&  $splice_acceptor eq "ag") {
		$m++;
	    }
        }
        if ($m == $array[9] - 1) {       # All the introns have cannonical splice sites
            print("$line\n");
        }
    }
    elsif ($line =~ m/\t\-\t/) {        ### working on minus strand
        for ($i = 0; $i < $array[9] - 1; $i++) {    #$array[9] is the number of exons. So, $array[9] - 1 is the number of introns
            $intron_coord2 = $array[1] + $col12[$i] + $col11[$i];
            $intron_coord1 = $array[1] + $col12[$i+1] - 1;
            $splice_acceptor = substr($hash{$array[0]}, $intron_coord2, 2);         #catching 2 nt from the intron start
            $splice_donor = substr($hash{$array[0]}, $intron_coord1 - 1, 2);  #catching the last 2 nt of the intron
            print STDERR ("$array[0]\t$array[3]\t$array[5]\t$splice_donor\t$splice_acceptor\n");
            if ($splice_donor eq "AC" && $splice_acceptor eq "CT" ) {
                $m++; 
            }
	    elsif ($splice_donor eq "ac" && $splice_acceptor eq "ct") {
		$m++;
	    }
        }
        if ($m == $array[9] - 1) {       # All the introns have cannonical splice sites
            print("$line\n");
        }
    }     
    $line = <INFILE>;
    chomp($line);
}
