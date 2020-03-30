#! bin/bash

num=$(ls -f1 results/kallisto/kallisto/aliquot/*/abundance.tsv | wc -l)
upper=$(echo "$((5 * $num))")
myseq=$(seq 4 5 $upper | sed 's/^\|$//g' | paste -sd,)
myseq=$(echo "1,2,"$myseq)
paste results/kallisto/kallisto/aliquot/*/abundance.tsv | cut -f $myseq > results/kallisto/kallisto/final/transcript_count_matrix_all_samples.tsv

ls -f1 results/kallisto/kallisto/aliquot/*/abundance.tsv | sed 's/results\/kallisto\/kallisto\/aliquot\///g' | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > results/kallisto/kallisto/final/header.tsv
cat results/kallisto/kallisto/final/header.tsv results/kallisto/kallisto/final/transcript_count_matrix_all_samples.tsv | grep -v "est_counts" > results/kallisto/kallisto/final/count_tpm_matrix_all_samples.tsv2
mv results/kallisto/kallisto/final/count_tpm_matrix_all_samples.tsv2 results/kallisto/kallisto/final/transcript_count_matrix_all_samples.tsv
rm -f results/kallisto/kallisto/final/header.tsv