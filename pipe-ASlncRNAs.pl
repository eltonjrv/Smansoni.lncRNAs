#!/usr/bin/perl
# Programmer: Elton Vasconcelos (last update: Sep/2016)
# Pipeline for identification of antisense lncRNAs that might be present in a transcriptome assembly (e.g. trinity contigs or tuxedo TCONS).
# This code was written by me as part of my postdoctoral research grant entitled "Identification and characterization of regulatory Long Non-coding RNAs on the Schistosoma mansoni genome through NGS strategies and systems approach", funded by Sao Paulo Research Foundation (FAPESP: 14/24560-8).
# If you use this whole pipeline or part of it, please cite Vasconcelos et. al., 2017 (https://www.nature.com/articles/s41598-017-10853-6) as well as its embedded bioinformatics tools individually.
# Five input files must be provided at the cmd-line: contigs.fasta, contigs.bed, annotated_genes.bed, repeats_library.fasta, ref-genome.fasta
###############################################################################################################################
# Usage: nohup perl pipe-ASlncRNAs.pl [infile.fasta] [infile.bed] [annotated_genes.bed] [repeats_library.fasta] [ref-genome.fasta] >nohup-pipe-ASlncRNAs.out 2>nohup-pipe-ASlncRNAs.err &
###############################################################################################################################
# IMPORTANT NOTES:
# Note 1: Do "git clone https://github.com/eltonjrv/Smansoni.lncRNAs/", change the directory to "Smansoni.lncRNAs", and then run this pipeline within that directory.
# Note 2: Please take a look at "software2install.txt" file and make sure you have all those tools installed at your workstation. If they are not set on your environment variables, please edit the lines below where they are invoked, typing the program full PATH.
# Note 3: Also take a look at the "Assembly of RNA-Seq reads" topic from the Supplementary Methods of our publication (if you have no experience on transcriptome assemblies) in order to generate both "contigs.fasta" and "contigs.bed" input files.
# Note 4: Accessories ad-hoc scripts written to support this pipeline are placed within "lncRNA-pipeTools/perl-scripts" dir.
# Note 5: One must edit line 9 from "lncRNA-pipeTools/perl-scripts/seqtools.pl" to correctly point to his/her BioPerl full PATH, as well as line 31 from "lncRNA-pipeTools/perl-scripts/seqs1.pl", replacing the seqtools.pl right location.
# Note 6: The chromosome IDs must be identical on the two bed input files as well as on the ref-genome.fasta provided in the command line
# Note 7: All the FASTA and BED output files generated during the execution of this pipeline are named intuitively, indicating the filtering steps that have been performed.
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
`java -jar lncRNA-pipeTools/readseq.jar -inform fasta -f fasta -o $base_fasta.fasta2 $ARGV[0]`;
`perl -pi -e \'s/\\\|c/_c/g\' $base_fasta.fasta2`;     #If it is a trinity assembly, the contig IDs have a "|" character that must be replaced by "_"
`perl -pi -e \'s/\\\|c/_c/g\' $ARGV[1]`;

### bedtools intersect: catching all features that overlap to annotated genes on the opposite strand with 33% of their (features being searched) coverage
# Make sure your annotated_genes.bed input file contains protein-coding genes only.
`bedtools intersect -S -loj -f 0.33 -a $ARGV[1] -b $ARGV[2] >$base_bed-AS_PCs.bedloj`;
`grep 'Smp_' $base_bed-AS_PCs.bedloj | cut -f 1,2,3,4,5,6,7,8,9,10,11,12 | uniq >$base_bed-AS_PCs.bed`;
`cut -f 4 $base_bed-AS_PCs.bed >$base_bed-AS_PCs.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_PCs.nam -seq $base_fasta.fasta2 >$base_fasta-AS_PCs.fasta`;

### RepeatMasker: Eliminating transposable elements and low complexity repeats from the dataset
`RepeatMasker -s -lib $ARGV[3] -x -gff -gc 35 -dir . -pa 8 $base_fasta-AS_PCs.fasta`;
`perl lncRNA-pipeTools/perl-scripts/RM-cov_cutoff.pl $base_fasta-AS_PCs.fasta.cat $base_fasta-AS_PCs.fasta 0.5 >$base_fasta-masked_gt50percent-Blocks.tab`;
`cut -f 1 $base_fasta-masked_gt50percent-Blocks.tab >$base_fasta-masked_gt50percent-Blocks.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -excl $base_fasta-masked_gt50percent-Blocks.nam -seq $base_fasta-AS_PCs.fasta >$base_fasta-AS_PCs-noRepeats.fasta`;
`grep \'>\' $base_fasta-AS_PCs-noRepeats.fasta | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-AS_PCs-noRepeats.nam`; 
`cat $base_bed-AS_PCs-noRepeats.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-AS_PCs-noRepeats.bed`;

#### Excluding Ribosomal RNAs
`perl ribopicker-standalone-0.4.3/ribopicker.pl -i 70 -c 50 -out_dir ./$base_fasta-RiboPickerOUT -f $base_fasta-AS_PCs-noRepeats.fasta -dbs rrnadb`;
`grep -P \'^>\' $base_fasta-RiboPickerOUT/*nonrrna.fa | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-AS_PCs-noRepeats-nonrrna.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_PCs-noRepeats-nonrrna.nam -seq $base_fasta-AS_PCs-noRepeats.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna.fasta`;
`cat $base_bed-AS_PCs-noRepeats-nonrrna.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-AS_PCs-noRepeats-nonrrna.bed`;

### Catching spliced only (at least one intron greater than 30 bp)
`grep -v -P \'\\\t0,*\$\' $base_bed-AS_PCs-noRepeats-nonrrna.bed >$base_bed-AS_PCs-noRepeats-nonrrna-spliced.bed`;
`perl lncRNA-pipeTools/perl-scripts/intronSizes_cutoff_onBed.pl $base_bed-AS_PCs-noRepeats-nonrrna-spliced.bed 30 >$base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.bed`;
`cut -f 4 $base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.bed >$base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.nam -seq $base_fasta-AS_PCs-noRepeats-nonrrna.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.fasta`;
### Catching the ones with canonical splice sites only (GT - AG)
`perl lncRNA-pipeTools/perl-scripts/catch_canonical_spliceSites.pl $base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.bed $ARGV[4] >$base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed`;
`cut -f 4 $base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed >$base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam -seq $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta`;

### Getorf
`getorf -noreverse -minsize 75 -find 0 -sequence $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta -outseq getorf_out.fasta`;
`perl lncRNA-pipeTools/perl-scripts/getorf-byPercentage.pl getorf_out.fasta 25 >getorf_out-gt25aaAND25cov.fasta`;
`grep \'>\' getorf_out-gt25aaAND25cov.fasta | sed \'s/ .*//g\' | sed \'s/_[0-9]*\$//g\' | sed \'s/^>//g\' | sort -u >withORFsgt25aaAND25percentCov.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outmft fasta -excl withORFsgt25aaAND25percentCov.nam -seq $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;

### CPC
`cpc-0.9-r2/bin/run_predict.sh $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta cpc-results.tab ./ cpc-evd`;
`grep -P \'\\\tcoding\\\t\' cpc-results.tab | cut -f 1 | sort -u >cpc-coding.nam`;

### TransDecoder
`TransDecoder-2.0.1/TransDecoder.LongOrfs -S -t $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`TransDecoder-2.0.1/TransDecoder.Predict -t $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`grep -P \'\\\tCDS\'  $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta.transdecoder.gff3 | cut -f 1 | sort -u >transDecoder-ORFs.nam`;

### Removing CPC and transDecoder coding predictions
`cat cpc-coding.nam transDecoder-ORFs.nam | sort -u >cpc-transDecoder-2remove.nam`;
`perl lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -excl cpc-transDecoder-2remove.nam -seq $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta`;
`grep  '>' $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta | sed \'s/>//g\' | sed \'s/ .*//g\' >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam`;
`cat $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam |  xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.bed`;

### InterproScan for removing sequences showing hits against either Pfam or PANTHER protein domains databases
## ATTENTION: Depending on the size of the input fasta file, InterProScan may take several hours or even days running.
# That's why its execution is "commented" below, so the user can appropriately split his/her input fasta file and go along with the protein domains detection independently.
# The commands below are suggestions on how to split the FASTA input, run interproscan, merge the results and remove the domain-hitting sequences
# $ ~/lncRNA-pipeTools/perl-scripts/split-FASTA.pl $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta 1000
# $ for i in `ls *--file*.fasta`; do nohup nice interproscan.sh -appl Panther,Pfam -t n -i $i -b `echo $i | sed 's/.*--file/lincRNAs-iprscan-part/g' | sed 's/\.fasta//g'` -T temp`echo $i | sed 's/.*--file//g' | sed 's/\.fasta//g'` -goterms -iprlookup &; done
### After all interproscan runs have finished, run the following:
# $ cat lincRNAs-iprscan-part*gff3 | grep -P '\tPANTHER\t.*\t\+\t|\tPfam\t.*\t\+\t' >lincRNAs-iprscan-PANTHER_Pfam-hits.gff3;
# $ cut -f 1 lincRNAs-iprscan-PANTHER_Pfam-hits.gff3 | sort -u >lincRNAs-iprscan-PANTHER_Pfam-hits.nam
# $ perl ~/lncRNA-pipeTools/perl-scripts/seqs1.pl -outfmt fasta -excl lincRNAs-iprscan-PANTHER_Pfam-hits.nam -seq $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD-noIprScan.fasta
# $ grep '>' $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD-noIprScan.fasta | sed \'s/>//g\' | sed \'s/ .*//g\' >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD-noIprScan.nam
# $ cat $base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD-noIprScan.nam |  xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_fasta-AS_PCs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD-noIprScan.bed
