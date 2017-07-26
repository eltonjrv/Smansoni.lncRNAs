#!/usr/bin/perl
# Programmer: Elton Vasconcelos (15/04/2015)
# Script that selects, on a bed input file, only genes that bear at least one intron with length greater then the user-defined cutoff
# Usage: perl intronSizes_cutoff_onBed.pl [infile.bed] [intron_length_cutoff_integer]  >outfile.bed
############################################################
# The bed output file will contain only the genes that have passed through the cutoff

open (INFILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my $line = <INFILE>;
chomp($line);
my (@array, @col11, @col12, $m);
my $cutoff = $ARGV[1];    #$ARGV[1] must be an integer (>= 1). That's the minimum intron length

while ($line ne "") {
    $m = 0;
    @array = split(/\t/, $line);
    $array[10] =~ s/\,$//g;    #only in case these columns are ending with ","
    $array[11] =~ s/\,$//g;
    @col11 = split(/\,/, $array[10]);
    @col12 = split(/\,/, $array[11]);
    for ($i = 0; $i < @col11; $i++) {
        if ($col12[$i+1] - ($col12[$i] + $col11[$i]) >= $cutoff) {
           $m = 1; 
        }        
    }
    if ($m == 1) {
        print("$line\n");
    }     
    $line = <INFILE>;
    chomp($line);
}