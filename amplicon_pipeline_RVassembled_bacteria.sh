#!/bin/bash

module load gparallel
module load Java/1.7.0_51
module load CDHIT/4.6.1c
module load R/2.15.1
module load Biopython/1.65


DIR=/mnt/research/germs/Schuyler/JRBP/Bacteria/germs_otus
RAWDAT_DIR=$DIR/reads

SCRIPTS=/mnt/home/smithsch/code
RDP=/mnt/research/rdp/public
PANDASEQ=$RDP/RDP_misc_tools/pandaseq/pandaseq
CHIMERA_DB=/mnt/research/germs/Schuyler/Databases/GreenGene/current_Bacteria_unaligned.fa
VSEARCH=/mnt/research/germs/Schuyler/Applications/bin/vsearch
CORES=2


mkdir -p $DIR/parallel_scripts/../cdhit_clustering/master_otus/../R/../combined_seqs/../mapped_samples/../otu_tables

cat $DIR/reads/* > $DIR/cdhit_clustering/combined_seqs/all_sequences.fa

cd-hit-est -d 0 -i $DIR/cdhit_clustering/combined_seqs/all_sequences.fa -o $DIR/cdhit_clustering/combined_seqs/all_seqs_cdhit.fasta -c 0.97 -M 200000 -T 20

python $SCRIPTS/rename_otus.py cdhit_clustering/combined_seqs/all_seqs_cdhit.fasta > cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta

for sample in $DIR/reads/*; 
	do name=`basename $sample`
	echo "cd-hit-est-2d -d 0 -c 0.97 -i $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta -i2 $sample -o $DIR/cdhit_clustering/mapped_samples/mapped_otus_$name -M 200000 -T $CORES"
done > $DIR/parallel_scripts/mapping_otus.sh
cat $DIR/parallel_scripts/mapping_otus.sh | parallel -j 10

for file in $DIR/cdhit_clustering/mapped_samples/*.clstr
	do name=`basename $file`
	echo "python $SCRIPTS/cdhit_otu_mapping_counts_only.py $file > $DIR/cdhit_clustering/otu_tables/$name.tbl"
done > $DIR/parallel_scripts/make_otu_tables.sh
cat $DIR/parallel_scripts/make_otu_tables.sh | parallel -j 20

cat $DIR/cdhit_clustering/otu_tables/* > $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt

Rscript $SCRIPTS/convert_otu_table_long_to_wide_format.R $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt $DIR/cdhit_clustering/R/cdhit_otu_table_wide.txt

java -Xmx24g -jar $RDP/RDPTools/classifier.jar classify -c 0.5 -f filterbyconf -o $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt -h $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf_hierarchy.txt $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta

cp $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt $DIR/cdhit_clustering/R/cdhit_taxa_table_w_repseq.txt

mv $DIR/cdhit_clustering/R/cdhit_taxa_table_w_repseq.txt $DIR/cdhit_clustering/R/bac_cdhit_taxa_table_w_repseq.txt
mv $DIR/cdhit_clustering/R/cdhit_otu_table_wide.txt $DIR/cdhit_clustering/R/bac_cdhit_otu_table_wide.txt 
sed -i 's/B//g' $DIR/cdhit_clustering/R/bac_cdhit_otu_table_wide.txt

