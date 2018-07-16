#!/bin/bash

module load gparallel
module load Java/1.7.0_51
module load CDHIT/4.6.1
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


mkdir -p $DIR/cdhit_clustering/master_otus/../R/../combined_seqs

python $SCRIPTS/renaming_seq_w_short_sample_name.py $DIR/cdhit_clustering/combined_seqs/sample_filename_map.txt $DIR/cdhit_clustering/combined_seqs/sequence_name_map.txt $RAWDAT_DIR/*.fasta > $DIR/cdhit_clustering/combined_seqs/all_sequences.fa

cd-hit-est -i $DIR/cdhit_clustering/combined_seqs/all_sequences.fa -o $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta -c 0.95 -M 200000 -T 20

python $SCRIPTS/cdhit_clstr_to_otu.py $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta.clstr > $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt

Rscript $SCRIPTS/convert_otu_table_long_to_wide_format.R $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt $DIR/cdhit_clustering/R/cdhit_otu_table_wide.txt

java -Xmx24g -jar $RDP/RDPTools/classifier.jar classify -c 0.5 -f filterbyconf -o $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt -h $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf_hierarchy.txt $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta

python $SCRIPTS/rep_seq_to_otu_mapping.py $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta.clstr > $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit_rep_seq_to_cluster.map

Rscript $SCRIPTS/renaming_taxa_rep_seq_to_otus.R $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit_rep_seq_to_cluster.map $DIR/cdhit_clustering/R/cdhit_taxa_table_w_repseq.txt

mv cdhit_taxa_table_w_repseq.txt bac_cdhit_taxa_table_w_repseq.txt
mv cdhit_otu_table_wide.txt bac_cdhit_otu_table_wide.txt 
sed -i 's/B//g' bac_cdhit_taxa_table_w_repseq.txt 
sed -i 's/B//g' bac_cdhit_otu_table_wide.txt
