#!/usr/bin/perl
# Programmer: Elton Vasconcelos (12/09/2016)
# Pipeline for an automatic detection of antisense lncRNAs that might be present in a transcriptome assembly (e.g. trinity contigs)
# Five input files must be provided at the cmd-line: contigs.fasta, contigs.bed, annotated_genes.bed, repeats_library.fasta, ref-genome.fasta
# Usage: nohup perl pipe-ASlncRNAs.pl [infile.fasta] [infile.bed] [annotated_genes.bed] [repeats_library.fasta] [ref-genome.fasta] >nohup-pipe-ASlncRNAs.out 2>nohup-pipe-ASlncRNAs.err &
###########################################################################################################################
# Note-1: At the initial stage, it has to be run on verjo-server-01
# Note-2: The chromosome IDs must be the same in all input files
# Note-3: The flanking regions from the annotated genes (on the first bedtools window execution) should also be declared on the cmd-line
# Note-4: Will have to set the repeats_library.fasta as an optional argument, as well as the GC content given on RepeatMasker cmd line (line 28)

my $base_fasta = $ARGV[0];
$base_fasta =~ s/\.\w+$//g;
my $base_bed = $ARGV[1];
$base_bed =~ s/\.\w+$//g;

### Running readseq.jar to add the sequence length at the end of each headline (>.+ \d+ bp$)
`java -jar /home/elton/bioinformatics-tools/readseq.jar -inform fasta -f fasta -o $base_fasta.fasta2 $ARGV[0]`;
`perl -pi -e \'s/\\\|c/_c/g\' $base_fasta.fasta2`;     #If it is a trinity assembly, the contig IDs have a "|" character that must be replaced by "_"
`perl -pi -e \'s/\\\|c/_c/g\' $ARGV[1]`;

### bedtools intersect: catching all features that overlap to annotated genes on the opposite strand with 33% of their (features being searched) coverage
`bedtools intersect -S -loj -f 0.33 -a $ARGV[1] -b $ARGV[2] >$base_bed-AS_SMPs.bedloj`;
`grep 'Smp_' $base_bed-AS_SMPs.bedloj | cut -f 1,2,3,4,5,6,7,8,9,10,11,12 | uniq >$base_bed-AS_SMPs.bed`;
`cut -f 4 $base_bed-AS_SMPs.bed >$base_bed-AS_SMPs.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_SMPs.nam -seq $base_fasta.fasta2 >$base_fasta-AS_SMPs.fasta`;

### RepeatMasker: Eliminating transposable elements and low complexity repeats from the dataset
`nice /home/elton/bioinformatics-tools/RepeatMasker/RepeatMasker -s -lib $ARGV[3] -x -gff -gc 35 -dir . -pa 2 $base_fasta-AS_SMPs.fasta`;
`perl /home/elton/bioinformatics-tools/perl-scripts/myIQUSP-scripts/RM-cov_cutoff.pl $base_fasta-AS_SMPs.fasta.cat $base_fasta-AS_SMPs.fasta 0.5 >$base_fasta-masked_gt50percent-Blocks.tab`;
`cut -f 1 $base_fasta-masked_gt50percent-Blocks.tab >$base_fasta-masked_gt50percent-Blocks.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -excl $base_fasta-masked_gt50percent-Blocks.nam -seq $base_fasta-AS_SMPs.fasta >$base_fasta-AS_SMPs-noRepeats.fasta`;
`grep \'>\' $base_fasta-AS_SMPs-noRepeats.fasta | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-AS_SMPs-noRepeats.nam`; 
`cat $base_bed-AS_SMPs-noRepeats.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-AS_SMPs-noRepeats.bed`;

#### Excluding Ribosomal RNAs
`nice perl /home/elton/bioinformatics-tools/ribopicker-standalone-0.4.3/ribopicker.pl -i 70 -c 50 -out_dir ./$base_fasta-RiboPickerOUT -f $base_fasta-AS_SMPs-noRepeats.fasta -dbs rrnadb`;
`grep -P \'^>\' $base_fasta-RiboPickerOUT/*nonrrna.fa | sed \'s/^>//g\' | sed \'s/ .*//g\' >$base_bed-AS_SMPs-noRepeats-nonrrna.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_SMPs-noRepeats-nonrrna.nam -seq $base_fasta-AS_SMPs-noRepeats.fasta >$base_fasta-AS_SMPs-noRepeats-nonrrna.fasta`;
`cat $base_bed-AS_SMPs-noRepeats-nonrrna.nam | xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_bed-AS_SMPs-noRepeats-nonrrna.bed`;

### Catching spliced only (at least one intron greater than 30 bp)
`grep -v -P \'\\\t0,*\$\' $base_bed-AS_SMPs-noRepeats-nonrrna.bed >$base_bed-AS_SMPs-noRepeats-nonrrna-spliced.bed`;
`perl /home/elton/bioinformatics-tools/perl-scripts/myIQUSP-scripts/intronSizes_cutoff_onBed.pl $base_bed-AS_SMPs-noRepeats-nonrrna-spliced.bed 30 >$base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.bed`;
`cut -f 4 $base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.bed >$base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.nam -seq $base_fasta-AS_SMPs-noRepeats-nonrrna.fasta >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.fasta`;
### Catching the ones with canonical splice sites only (GT - AG)
`perl /home/elton/bioinformatics-tools/perl-scripts/myIQUSP-scripts/catch_canonical_spliceSites.pl $base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.bed $ARGV[4] >$base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed`;
`cut -f 4 $base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.bed >$base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -incl $base_bed-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.nam -seq $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30.fasta >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta`;

### Getorf
`getorf -noreverse -minsize 75 -find 0 -sequence $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta -outseq getorf_out.fasta`;
`perl /home/elton/bioinformatics-tools/perl-scripts/myIQUSP-scripts/getorf-byPercentage.pl getorf_out.fasta 25 >getorf_out-gt25aaAND25cov.fasta`;
`grep \'>\' getorf_out-gt25aaAND25cov.fasta | sed \'s/ .*//g\' | sed \'s/_[0-9]*\$//g\' | sed \'s/^>//g\' | sort -u >withORFsgt25aaAND25percentCov.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outmft fasta -excl withORFsgt25aaAND25percentCov.nam -seq $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice.fasta >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;

### CPC
`/home/elton/bioinformatics-tools/cpc-0.9-r2/bin/run_predict.sh $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta cpc-results.tab ./ cpc-evd`;
`grep -P \'\\\tcoding\\\t\' cpc-results.tab | cut -f 1 | sort -u >cpc-coding.nam`;

### TransDecoder
`/home/elton/bioinformatics-tools/TransDecoder-2.0.1/TransDecoder.LongOrfs -S -t $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`/home/elton/bioinformatics-tools/TransDecoder-2.0.1/TransDecoder.Predict -t $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta`;
`grep -P \'\\\tCDS\'  $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta.transdecoder.gff3 | cut -f 1 | sort -u >transDecoder-ORFs.nam`;

### Removing CPC and transDecoder coding predictions
`cat cpc-coding.nam transDecoder-ORFs.nam | sort -u >cpc-transDecoder-2remove.nam`;
`perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -excl cpc-transDecoder-2remove.nam -seq $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs.fasta >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta`;
`grep  '>' $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.fasta | sed \'s/>//g\' | sed \'s/ .*//g\' >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam`;
`cat $base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.nam |  xargs -i grep -P \'{}\\\t\' $ARGV[1] >$base_fasta-AS_SMPs-noRepeats-nonrrna-spliced-intron_gt30-canonicalSplice-noORFs-noCPC_TD.bed`;

