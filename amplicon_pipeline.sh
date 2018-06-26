#!/bin/bash

module load gparallel
module load Java/1.7.0_51
module load CDHIT/4.6.1
module load R/2.15.1
module load Biopython/1.65


DIR=/mnt/home/smithsch/Yeast
RAWDAT_DIR=$DIR/original
MAPPING_FILE=$RAWDAT_DIR/map.txt

SCRIPTS=/mnt/home/smithsch/code
RDP=/mnt/research/rdp/public
PANDASEQ=$RDP/RDP_misc_tools/pandaseq/pandaseq
CHIMERA_DB=/mnt/scratch/smithsch/DataBase_GreenGene/current_Bacteria_unaligned.fa
VSEARCH=/mnt/research/germs/Schuyler/Applications/bin/vsearch
CORES=2

## assemble paired-ends. The below parameters work well with bacterial 16S. 
OVERLAP=10 ## minimal number of overlapped bases required for pair-end assembling. Not so critical if you set the length parameters (see below)
MINLENGTH=480 #16s: "250" ## minimal length of the assembled sequence
MAXLENGTH=520 #16s: "280" ## maximum length of the assembled sequence
Q=25 ## minimal read quality score.


ls $RAWDAT_DIR/*_R*.fastq* > $RAWDAT_DIR/raw_seq_list.txt
mkdir -p $DIR/parallel_scripts/../split_fastqs
while read SEQS; 
	do NAME=`basename $SEQS | cut -d "." -f 1`; 
	mkdir -p $DIR/split_fastqs/$NAME; 
	echo "split -l 1000000 $SEQS $DIR/split_fastqs/$NAME/"; 
done < $RAWDAT_DIR/raw_seq_list.txt > $DIR/parallel_scripts/split_fastqs.sh
cat $DIR/parallel_scripts/split_fastqs.sh | parallel -j 4
ls $DIR/split_fastqs/*/* | rev | cut -d "/" -f 1 | sort -u | rev > $DIR/seqs_list.txt


mkdir -p $DIR/pandaseq/assembled/../stats
while read SEQS;
	do echo "$PANDASEQ -T $CORES -o $OVERLAP -e $Q -N -F -d rbkfms -l $MINLENGTH -L $MAXLENGTH -f $DIR/split_fastqs/*_R1_*/$SEQS -r $DIR/split_fastqs/*_R2_*/$SEQS 1> $DIR/pandaseq/assembled/$SEQS.assembled.fastq 2> $DIR/pandaseq/stats/$SEQS.assembled.stats.txt.bz2";
done < $DIR/seqs_list.txt > $DIR/parallel_scripts/pandaseq.sh
cat $DIR/parallel_scripts/pandaseq.sh | parallel -j 10
find $DIR/pandaseq/assembled -type f -size 0 -exec rm {} +
ls $DIR/pandaseq/assembled/* | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1 > $DIR/seqs_list.txt

mkdir -p $DIR/quality_check/seqs_25/../chimera_removal/../final_seqs
while read SEQS;
	do echo "java -jar $RDP/RDPTools/SeqFilters.jar -Q $Q -s $DIR/pandaseq/assembled/$SEQS.assembled.fastq -o $DIR/quality_check/seqs_25 -O $SEQS; python $SCRIPTS/fastq_to_fasta.py $DIR/quality_check/seqs_25/$SEQS/NoTag/NoTag_trimmed.fastq $DIR/quality_check/seqs_25/$SEQS/$SEQS.fa";
done < $DIR/seqs_list.txt > $DIR/parallel_scripts/qc.sh
cat $DIR/parallel_scripts/qc.sh | parallel -j 20 --delay 2
cat $DIR/quality_check/seqs_25/*/*.fa > $DIR/quality_check/chimera_removal/all_combined_q25.fa

$VSEARCH --threads 20 --derep_fulllength $DIR/quality_check/chimera_removal/all_combined_q25.fa --output $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2.fa --sizeout --minuniquesize 5
$VSEARCH --threads 20 --uchime_denovo $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2.fa --chimeras $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.chimera --nonchimeras $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.good
$VSEARCH --threads 20 --uchime_ref $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.good --nonchimeras $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo_ref.good --db $CHIMERA_DB

while read SEQS;
	do echo "$VSEARCH --usearch_global $DIR/quality_check/seqs_25/$SEQS/$SEQS.fa --db $DIR/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo_ref.good --id 0.985 --matched $DIR/quality_check/final_seqs/$SEQS\_finalized.fa";
done < $DIR/seqs_list.txt > $DIR/parallel_scripts/remapping.sh
cat $DIR/parallel_scripts/remapping.sh | parallel -j 20

mkdir -p $DIR/i_file/../demultiplex/parse_index/../bins/../empty_samples/../demultiplex_finalized; split -l 1000000 $RAWDAT_DIR/*_I*.fastq* $DIR/i_file/
perl $SCRIPTS/dos2unix.pl $MAPPING_FILE > $DIR/demultiplex/tag_file.txt
python $SCRIPTS/MiSeq_rdptool_map_parser.py $DIR/demultiplex/tag_file.txt > $DIR/demultiplex/tag_file.tag

while read I;
	do mkdir $DIR/demultiplex/parse_index/$I
	echo "java -jar $RDP/RDPTools/SeqFilters.jar --seq-file $DIR/i_file/$I --tag-file $DIR/demultiplex/tag_file.tag --outdir $DIR/demultiplex/parse_index/$I"
done < $DIR/seqs_list.txt > $DIR/parallel_scripts/demultiplex.sh
cat $DIR/parallel_scripts/demultiplex.sh | parallel -j 14
awk '{print $1}' $DIR/demultiplex/tag_file.txt | tail -n +2 > $DIR/lane_list.txt

while read LANE;
	do echo "cat $DIR/demultiplex/parse_index/*/result_dir/$LANE/$LANE\_trimmed.fastq > $DIR/demultiplex/bins/$LANE\_trimmed.fastq"
done < $DIR/lane_list.txt > $DIR/parallel_scripts/cat_lanes.sh
cat $DIR/parallel_scripts/cat_lanes.sh | parallel -j 20


# here change lane_list to include only the ones you want in the analysis.

while read SEQS;
	do mkdir -p $DIR/demultiplex/parsed_fastas/$SEQS
	echo "python $SCRIPTS/bin_reads.py $DIR/quality_check/final_seqs/$SEQS\_finalized.fa $DIR/demultiplex/bins $DIR/demultiplex/parsed_fastas/$SEQS"
done < $DIR/seqs_list.txt > $DIR/parallel_scripts/bin_reads.sh
cat $DIR/parallel_scripts/bin_reads.sh | parallel -j 20

while read I;
	do echo "cat $DIR/demultiplex/parsed_fastas/*/$I*.fasta > $DIR/demultiplex/demultiplex_finalized/$I.fasta"
done < $DIR/lane_list.txt > $DIR/parallel_scripts/cat_bins.sh
cat $DIR/parallel_scripts/cat_bins.sh | parallel -j 20
find $DIR/demultiplex/demultiplex_finalized -type f -size 0 -exec mv -t $DIR/demultiplex/empty_samples/ {} +


mkdir -p $DIR/cdhit_clustering/master_otus/../R/../combined_seqs
python $SCRIPTS/renaming_seq_w_short_sample_name.py $DIR/cdhit_clustering/combined_seqs/sample_filename_map.txt $DIR/cdhit_clustering/combined_seqs/sequence_name_map.txt $DIR/demultiplex/demultiplex_finalized/*.fasta > $DIR/cdhit_clustering/combined_seqs/all_sequences.fa
cd-hit-est -i $DIR/cdhit_clustering/combined_seqs/all_sequences.fa -o $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta -c 0.95 -M 200000 -T 20
python $SCRIPTS/cdhit_clstr_to_otu.py $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta.clstr > $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt
Rscript $SCRIPTS/convert_otu_table_long_to_wide_format.R $DIR/cdhit_clustering/master_otus/cdhit_otu_table_long.txt $DIR/cdhit_clustering/R/cdhit_otu_table_wide.txt

java -Xmx24g -jar $RDP/RDPTools/classifier.jar classify -c 0.5 -f filterbyconf -o $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt -h $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf_hierarchy.txt $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta

python $SCRIPTS/rep_seq_to_otu_mapping.py $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit.fasta.clstr > $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit_rep_seq_to_cluster.map
Rscript $SCRIPTS/renaming_taxa_rep_seq_to_otus.R $DIR/cdhit_clustering/master_otus/cdhit_otu_taxa_filterbyconf.txt $DIR/cdhit_clustering/combined_seqs/combined_seqs_cdhit_rep_seq_to_cluster.map $DIR/cdhit_clustering/R/cdhit_taxa_table_w_repseq.txt


python $SCRIPTS/cdhit_clstr_to_otu.py otus_repsetOUT.fasta > test.txt
