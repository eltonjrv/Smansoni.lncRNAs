#!/usr/bin/perl
# Programmer: Elton Vasconcelos (last update: Sep/2016).
# Pipeline for identification of novel putative lincRNAs that might be present in a transcriptome assembly (e.g. trinity contigs or tuxedo TCONS).
# This code was written by me as part of my postdoctoral research grant entitled "Identification and characterization of regulatory Long Non-coding RNAs on the Schistosoma mansoni genome through NGS strategies and systems approach", funded by Sao Paulo Research Foundation (FAPESP: 14/24560-8).
# If you use this whole pipeline or part of it, please cite this github page as well as its embedded bioinformatics tools individually.
# Our manuscript is currently under review by Scientific Reports (a NPG Journal).
# Five input files must be provided at the cmd-line: contigs.fasta, contigs.bed, annotated_genes.bed, repeats_library.fasta, ref-genome.fasta
###############################################################################################################################
# Usage: nohup perl pipe-lincRNAs.pl [infile.fasta] [infile.bed] [annotated_genes.bed] [repeats_library.fasta] [ref-genome.fasta] >nohup-pipe-lincRNAs.out 2>nohup-pipe-lincRNAs.err &
################################################################################################################################
# IMPORTANT NOTES:
# Note 1: 'git clone' the "lncRNA-pipeTools" branch from this repository to your "home" folder at your workstation.
# Note 2: Please have bedtools (http://bedtools.readthedocs.io/en/latest/) and EMBOSS suite (http://emboss.sourceforge.net/) installed and placed at your environment variables.
# Note 3: The chromosome IDs must be identical in all input files provided.
#################################################################################################################
# This program and its embedded tools are free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details at http://www.gnu.org/licenses/.

if (@ARGV != 5) {
die ("**Error**\nThe cmd line must contain 5 arguments:\n\$perl pipe-ASlncRNAs.pl [infile.fasta] [infile.bed] [annotated_genes.bed] [repeats_library.fasta] [ref-genome.fasta]\nRead script's initial commented lines for a better explanation\n");
}

my $base_fasta = $ARGV[0];
$base_fasta =~ s/\.\w+$//g;
my $base_bed = $ARGV[1];
$base_bed =~ s/\.\w+$//g;

### Running readseq.jar to add the sequence length at the end of each headline (>.+ \d+ bp$)
`java -jar ~/lncRNA-pipeTools/readseq.jar -inform fasta -f fasta -o $base_fasta.fasta2 $ARGV[0]`;
`perl -pi -e \'s/\\\|c/_c/g\' $base_fasta.fasta2`;     #If it is a trinity assembly, the contig IDs have a "|" character that must be replaced by "_"
`perl -pi -e \'s/\\\|c/_c/g\' $ARGV[1]`;

### bedtools intersect: catching all features that do not overlap to annotated genes
`bedtools intersect -v -a $ARGV[1] -b $ARGV[2] >$base_bed-nonOverlapped2PCs.bed`;
`cut -f 4 $base_bed-nonOverlapped2PCs.bed >$base_bed-nonOverlapped2PCs.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-nonOverlapped2PCs.nam -seq $base_fasta.fasta2 >$base_fasta-nonOverlapped2PCs.fasta`;

### RepeatMasker: Eliminating transposable elements and low complexity repeats from the dataset
`~/lncRNA-pipeTools/RepeatMasker/RepeatMasker -s -lib $ARGV[3] -x -gff -gc 35 -dir . -pa 8 $base_fasta-nonOverlapped2PCs.fasta`;
`perl ~/lncRNA-pipeTools/perl-scripts/myIQUSP-scripts/RM-cov_cutoff.pl $base_fasta-nonOverlapped2PCs.fasta.cat $base_fasta-nonOverlapped2PCs.fasta 0.5 >$base_fasta-masked_gt50percent-Blocks.tab`;
`cut -f 1 $base_fasta-masked_gt50percent-Blocks.tab >$base_fasta-masked_gt50percent-Blocks.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -excl $base_fasta-masked_gt50percent-Blocks.nam -seq $base_fasta-nonOverlapped2PCs.fasta >$base_fasta-nonOverlapped2PCs-noRepeats.fasta`;
`grep \'>\' $base_fasta-nonOverlapped2PCs-noRepeats.fasta | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-nonOverlapped2PCs-noRepeats.nam`; 
`cat $base_bed-nonOverlapped2PCs-noRepeats.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-nonOverlapped2PCs-noRepeats.bed`;

#### Excluding Ribosomal RNAs
`perl ~/lncRNA-pipeTools/ribopicker-standalone-0.4.3/ribopicker.pl -i 70 -c 50 -out_dir ./$base_fasta-RiboPickerOUT -f $base_fasta-nonOverlapped2PCs-noRepeats.fasta -dbs rrnadb`;
`grep -P \'^>\' $base_fasta-RiboPickerOUT/*nonrrna.fa | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-nonOverlapped2PCs-noRepeats-nonrrna.nam -seq $base_fasta-nonOverlapped2PCs-noRepeats.fasta >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna.fasta`;
`cat $base_bed-nonOverlapped2PCs-noRepeats-nonrrna.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna.bed`;

### Catching spliced only (at least one intron greater than 30 bp)
`grep -v -P \'\\\t0,*\$\' $base_bed-nonOverlapped2PCs-noRepeats-nonrrna.bed >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced.bed`;
`perl ~/lncRNA-pipeTools/perl-scripts/myIQUSP-scripts/intronSizes_cutoff_onBed.pl $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced.bed 30 >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.bed`;
`cut -f 4 $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.bed >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.nam -seq $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna.fasta >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.fasta`;
### Catching the ones with canonical splice sites only (GT - AG)
`perl ~/lncRNA-pipeTools/perl-scripts/myIQUSP-scripts/catch_canonical_spliceSites.pl $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.bed $ARGV[4] >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed`;
`cut -f 4 $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed >$base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam -seq $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30.fasta >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta`;

### Getorf
`getorf -noreverse -minsize 75 -find 0 -sequence $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta -outseq getorf_out.fasta`;
`perl ~/lncRNA-pipeTools/perl-scripts/myIQUSP-scripts/getorf-byPercentage.pl getorf_out.fasta 25 >getorf_out-gt25aaAND25cov.fasta`;
`grep \'>\' getorf_out-gt25aaAND25cov.fasta | sed \'s/ .*//g\' | sed \'s/_[0-9]*\$//g\' | sed \'s/^>//g\' | sort -u >withORFsgt25aaAND25percentCov.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outmft fasta -excl withORFsgt25aaAND25percentCov.nam -seq $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;

### CPC
`~/lncRNA-pipeTools/cpc-0.9-r2/bin/run_predict.sh $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta cpc-results.tab ./ cpc-evd`;
`grep -P \'\\\tcoding\\\t\' cpc-results.tab | cut -f 1 | sort -u >cpc-coding.nam`;

### TransDecoder
`~/lncRNA-pipeTools/TransDecoder-2.0.1/TransDecoder.LongOrfs -S -t $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`~/lncRNA-pipeTools/TransDecoder-2.0.1/TransDecoder.Predict -t $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`grep -P \'\\\tCDS\'  $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta.transdecoder.gff3 | cut -f 1 | sort -u >transDecoder-ORFs.nam`;

### Removing CPC and transDecoder coding predictions
`cat cpc-coding.nam transDecoder-ORFs.nam | sort -u >cpc-transDecoder-2remove.nam`;
`perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -excl cpc-transDecoder-2remove.nam -seq $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta`;
`grep '>' $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta | sed \'s/>//g\' | sed \'s/ .*//g\' >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam`;
`cat $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam |  xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.bed`;

### InterproScan
`~/lncRNA-pipeTools/perl-scripts/split-FASTA.pl $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta 2171`;
`nohup nice interproscan.sh -appl Panther,Pfam -t n -i $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta--file1.fasta -b lincRNAs-iprscan-part1 -T temp01 -goterms -iprlookup &`;
`nohup nice interproscan.sh -appl Panther,Pfam -t n -i $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta--file2.fasta -b lincRNAs-iprscan-part2 -T temp02 -goterms -iprlookup &`;
`nohup nice interproscan.sh -appl Panther,Pfam -t n -i $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta--file3.fasta -b lincRNAs-iprscan-part3 -T temp03 -goterms -iprlookup &`;
`nohup nice interproscan.sh -appl Panther,Pfam -t n -i $base_fasta-nonOverlapped2PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta--file4.fasta -b lincRNAs-iprscan-part4 -T temp04 -goterms -iprlookup &`;


