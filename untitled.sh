module load parallel
module load Python/2.7.14
pip install --user Biopython

DIRECTORY=/mnt/research/germs/Schuyler/Projects/Xuewei
RAWDAT_DIRECTORY=$DIRECTORY/original
MAPPING_FILE=$RAWDAT_DIRECTORY/map.txt

SCRIPTS=/mnt/research/germs/Schuyler/code
RDP=/mnt/research/rdp/public
PANDASEQ=$RDP/RDP_misc_tools/pandaseq/pandaseq
CHIMERA_DB=/mnt/scratch/smithsch/DataBase_GreenGene/current_Bacteria_unaligned.fa
VSEARCH=/mnt/research/germs/Schuyler/Applications/bin/vsearch
CORES=2

## assemble paired-ends. The below parameters work well with bacterial 16S. 
OVERLAP=10 ## minimal number of overlapped bases required for pair-end assembling. Not so critical if you set the length parameters (see below)
MINLENGTH=200 #16s: "250" ## minimal length of the assembled sequence
MAXLENGTH=280 #16s: "280" ## maximum length of the assembled sequence
Q=25 ## minimal read quality score.


mkdir -p $DIRECTORY/parallel_scripts/../i_file/../demultiplex/parse_index/../bins/../empty_samples/../demultiplex_finalized; split -l 1000000 $RAWDAT_DIRECTORY/*_I*.fastq* $DIRECTORY/i_file/
ls $DIRECTORY/i_file/* | rev | cut -d "/" -f 1 | sort -u | rev > $DIRECTORY/seqs_list.txt

mkdir -p $DIRECTORY/demultiplex/parse_index/../bins/../empty_samples/../demultiplex_finalized
perl $SCRIPTS/dos2unix.pl $MAPPING_FILE > $DIRECTORY/demultiplex/tag_file.txt
python $SCRIPTS/MiSeq_rdptool_map_parser.py $DIRECTORY/demultiplex/tag_file.txt > $DIRECTORY/demultiplex/tag_file.tag
### tag_file input must be (barcode \t sample_ID) below line swaps column position if mapping file in opposite order. 
# awk '{ print $2 " " $1}' $MAPPING_FILE > $DIRECTORY/demultiplex/tag_file.tag


while read I;
	do mkdir $DIRECTORY/demultiplex/parse_index/$I
	echo "java -jar $RDP/RDPTools/SeqFilters.jar --seq-file $DIRECTORY/i_file/$I --tag-file $DIRECTORY/demultiplex/tag_file.tag --outdir $DIRECTORY/demultiplex/parse_index/$I"
done < $DIRECTORY/seqs_list.txt > $DIRECTORY/parallel_scripts/demultiplex.sh
cat $DIRECTORY/parallel_scripts/demultiplex.sh | parallel -j 20

awk '{print $2}' $DIRECTORY/demultiplex/tag_file.tag | tail -n +2 > seqs_list.txt
while read LANE;
	do echo "cat $DIRECTORY/demultiplex/parse_index/*/result_dir/$LANE/$LANE\_trimmed.fastq > $DIRECTORY/demultiplex/bins/$LANE\_trimmed.fastq"
done < $DIRECTORY/seqs_list.txt > $DIRECTORY/parallel_scripts/cat_lanes.sh
cat $DIRECTORY/parallel_scripts/cat_lanes.sh | parallel -j 20
ls $DIRECTORY/demultiplex/bins/* | rev | cut -d "/" -f 1 | sort -u | rev | cut -d "_" -f 1 > $DIRECTORY/seqs_list.txt


mkdir -p $DIRECTORY/pandaseq/assembled/../stats
$PANDASEQ -T 20 -o $OVERLAP -e $Q -N -F -d rbkfms -l $MINLENGTH -L $MAXLENGTH -f $RAWDAT_DIRECTORY/*_R1_* -r $RAWDAT_DIRECTORY/*_R2_* 1 > $DIRECTORY/pandaseq/assembled/sequences.assembled.fastq 2> $DIRECTORY/pandaseq/stats/sequences.assembled.stats.txt.bz2
find $DIRECTORY/pandaseq/assembled -type f -size 0 -exec rm {} +
mkdir -p $DIRECTORY/quality_check/seqs_25/../chimera_removal/../final_seqs
java -jar $RDP/RDPTools/SeqFilters.jar -Q $Q -s $DIRECTORY/pandaseq/assembled/sequences.assembled.fastq -o $DIRECTORY/quality_check/seqs_25 -O sequences
python $SCRIPTS/fastq_to_fasta.py $DIRECTORY/quality_check/seqs_25/sequences/NoTag/NoTag_trimmed.fastq $DIRECTORY/quality_check/seqs_25/sequences/sequences.fa


while read SEQS;
	do mkdir -p $DIRECTORY/demultiplex/demultiplex_finalized/$SEQS
	echo "python $SCRIPTS/bin_reads.py $DIRECTORY/quality_check/seqs_25/sequences/sequences.fa"
done < $DIRECTORY/seqs_list.txt > $DIRECTORY/parallel_scripts/bin_reads.sh
cat $DIRECTORY/parallel_scripts/bin_reads.sh | parallel -j 20
mv $DIRECTORY/demultiplex/bins/* $DIRECTORY/demultiplex/demultiplex_finalized/$SEQS
find $DIRECTORY/demultiplex/demultiplex_finalized -type f -size 0 -exec mv -t $DIRECTORY/demultiplex/empty_samples/ {} +
ls $DIRECTORY/demultiplex/demultiplex_finalized/* | rev | cut -d "/" -f 1 | sort -u | rev | cut -d "_" -f 1 > $DIRECTORY/seqs_list.txt


cat $DIRECTORY/demultiplex/demultiplex_finalized/* > $DIRECTORY/quality_check/chimera_removal/all_combined_q25.fa

$VSEARCH --threads 20 --derep_fulllength $DIRECTORY/quality_check/chimera_removal/all_combined_q25.fa --output $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2.fa --sizeout --minuniquesize 2
$VSEARCH --threads 20 --uchime_denovo $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2.fa --chimeras $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.chimera --nonchimeras $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.good
$VSEARCH --threads 20 --uchime_ref $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo.good --nonchimeras $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo_ref.good --db $CHIMERA_DB
$VSEARCH --threads 20 --derep_fulllength $DIRECTORY/quality_check/chimera_removal/all_combined_q25_unique_sort_min2_denovo_ref.good --relabel "U_" --output $DIRECTORY/quality_check/chimera_removal/relabeled_denovo_ref.good


mkdir -p $DIRECTORY/cdhit_clustering/master_otus/../R/../combined_seqs/../renamed_seqs
$(CDHIT)/cd-hit-est -i $DIRECTORY/quality_check/chimera_removal/relabeled_denovo_ref.good -o $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99 -c 0.99 -M 200000 -T 20

python $SCRIPTS/renaming_seq_w_short_sample_name.py "S_" $DIRECTORY/cdhit_clustering/renamed_seqs/sample_filename_map.txt $DIRECTORY/cdhit_clustering/renamed_seqs/sequence_name_map.txt $DIRECTORY/demultiplex/demultiplex_finalized/*.fasta > $DIRECTORY/cdhit_clustering/renamed_seqs/all_renamed_sequences.fa

$(CDHIT)/cd-hit-est-2d -i $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99 -i2 $DIRECTORY/cdhit_clustering/renamed_seqs/all_renamed_sequences.fa -o $DIRECTORY/cdhit_clustering/master_otus/otu_mapping_cdhit99 -c 0.99 -M 200000 -T 20

python $SCRIPTS/cdhit_otu_mapping.py $DIRECTORY/cdhit_clustering/master_otus/otu_mapping_cdhit99.clstr > $DIRECTORY/cdhit_clustering/master_otus/cdhit_otu_table_long.txt
Rscript $SCRIPTS/convert_otu_table_long_to_wide_format.R $DIRECTORY/cdhit_clustering/master_otus/cdhit_otu_table_long.txt $DIRECTORY/cdhit_clustering/R/cdhit_otu_table_wide.txt
python $SCRIPTS/rep_seq_to_otu_mapping.py $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99.clstr > $DIRECTORY/cdhit_clustering/master_otus/rep_seq_to_cluster.map 

java -Xmx24g -jar $RDP/classifier.jar classify -c 0.5 -f filterbyconf -o $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99_taxa_filterbyconf.txt -h $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99_taxa_filterbyconf_hierarchy.txt $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99
Rscript $SCRIPTS/renaming_taxa_rep_seq_to_otus.R $DIRECTORY/cdhit_clustering/master_otus/relabeled_denovo_ref_good_cdhit_99_taxa_filterbyconf.txt $DIRECTORY/cdhit_clustering/master_otus/rep_seq_to_cluster.map $DIRECTORY/cdhit_clustering/R/cdhit_taxa_table_w_repseq.txt

