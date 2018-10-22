DIR= /mnt/research/germs/STRIPS_16S_data
ORI= /mnt/research/germs/STRIPS_16S_data/original
SUBPROJECT= soilcosm

CODE= /mnt/home/yangfan1/repos/amplicon_pipelines/code
RDP= /mnt/research/rdp/public
CHIMERA_DB= /mnt/home/yangfan1/databases/current_Bacteria_unaligned.fa

VSEARCH= /mnt/research/germs/softwares/vsearch-2.4.3-linux-x86_64/bin/vsearch

FN_DELIM= _
FN_REV_INDEX= 3
FN_END= _001.fastq.gz

## assemble paired-ends. The below parameters work well with bacterial 16S. 
OVERLAP= 10## minimal number of overlapped bases required for pair-end assembling. Not so critical if you set the length parameters (see below)
MINLEN= 250#16s: "250" ## minimal length of the assembled sequence
MAXLEN= 280#16s: "280" ## maximum length of the assembled sequence
Q= 25## minimal read quality score.

all: sample_names.txt 1_rdp_pandaseq assemble 2_quality_check quality chimera_prep fq_to_fa combine_derep_sort chimera_denovo chimera_ref remapping 3_cdhit_clustering renaming clustering otu_table taxa_table

assemblage: sample_names.txt 1_rdp_pandaseq assemble
quality_check: 2_quality_check quality chimera_prep fq_to_fa combine_derep_sort chimera_denovo chimera_ref remapping
cdhit_clustering: 3_cdhit_clustering renaming clustering otu_table taxa_table

.PHONY: clean

# 1. create a list of sample names. find the sample name shared between R1 and R2.
sample_names.txt: $(ORI)
	cd $(ORI) && ls *_R*.fastq* > ../raw_seq_list.txt && rev ../raw_seq_list.txt | cut -d $(FN_DELIM) -f "$(FN_REV_INDEX)"- | rev | sort | uniq > ../sample_names.txt

# 2. assemble pair-ended readsi using RDP pandaseq.
1_rdp_pandaseq:
	mkdir 1_rdp_pandaseq 1_rdp_pandaseq/assembled 1_rdp_pandaseq/stats

assemble:
	while read i; do \
		echo "$(RDP)/RDP_misc_tools/pandaseq/pandaseq -N -o $(OVERLAP) -e $(Q) -F -d rbfkms -l $(MINLEN) -L $(MAXLEN) -f $(ORI)/"$$i""$(FN_DELIM)"R1"$(FN_END)" -r  $(ORI)/"$$i""$(FN_DELIM)"R2"$(FN_END)" 1> 1_rdp_pandaseq/assembled/"$$i"_"$(MINLEN)"-"$(MAXLEN)".fastq 2> 1_rdp_pandaseq/stats/"$$i"_assembled_"$(MINLEN)"-"$(MAXLEN)"_stats.txt" ; \
	done < $(DIR)/sample_names.txt > assembling_script.sh; \
	true
	cat assembling_script.sh | parallel -j 4	

# 3. check sequence quality.
2_quality_check:
	mkdir 2_quality_check 2_quality_check/temp 2_quality_check/fastq_q25 2_quality_check/fasta_q25 2_quality_check/chimera_removal 2_quality_check/final_good_seqs

quality:
	cd $(DIR)/1_rdp_pandaseq/assembled && \
	for i in *.fastq; do \
		echo "java -jar $(RDP)/RDPTools/SeqFilters.jar -Q $(Q) -s $$i -o $(DIR)/2_quality_check/temp/ -O $$i.q25"; \
	done > $(DIR)/quality_script.sh; \
	true && \
	cat $(DIR)/quality_script.sh | parallel -j 4

chimera_prep:
	cd $(DIR)/2_quality_check/temp && \
	for i in *.q25; do \
		echo "mv $$i/NoTag/NoTag_trimmed.fastq ../fastq_q25/$${i//.fastq.q25/.q25}.fq"; \
	done > $(DIR)/chimera_prep.sh; \
	true && \
	cat $(DIR)/chimera_prep.sh | parallel -j 4 --delay 2
	rm -r $(DIR)/2_quality_check/temp

fq_to_fa:
	cd $(DIR)/2_quality_check/fastq_q25 && \
	for i in *.fq; do \
		echo "python $(CODE)/fastq_to_fasta.py $$i ../fasta_q25/$${i//.fq/}.fa"; \
	done > $(DIR)/fq_to_fa.sh; \
	true && \
	cat $(DIR)/fq_to_fa.sh | parallel -j 4

combine_derep_sort:
	cd $(DIR)/2_quality_check/fasta_q25 && \
	cat *.fa >> $(DIR)/2_quality_check/chimera_removal/all_combined_q25.fa
	cd $(DIR)/2_quality_check/chimera_removal && \
	$(VSEARCH) --derep_fulllength all_combined_q25.fa --output all_combined_q25_unique_sort_min2.fa --sizeout --minuniquesize 2

chimera_denovo:
	cd $(DIR)/2_quality_check/chimera_removal && \
	$(VSEARCH) --uchime_denovo all_combined_q25_unique_sort_min2.fa --chimeras all_combined_q25_unique_sort_min2_denovo.chimera --nonchimeras all_combined_q25_unique_sort_min2_denovo.good

## submitted: 256G mem, used 2 hours for 420773 unique sequences.
chimera_ref:
	cd $(DIR)/2_quality_check/chimera_removal && \
	$(VSEARCH) --uchime_ref all_combined_q25_unique_sort_min2_denovo.good --nonchimeras all_combined_q25_unique_sort_min2_denovo_ref.good --db $(CHIMERA_DB)

3_cdhit_clustering:
        mkdir $(DIR)/3_cdhit_clustering $(DIR)/3_cdhit_clustering/renamed_seqs $(DIR)/3_cdhit_clustering/master_otus $(DIR)/3_cdhit_clustering/R

master:
        cd $(DIR)/3_cdhit_clustering/master_otus && \
        $(VSEARCH) --derep_fulllength $(DIR)/2_quality_check/denovo_ref/combined_denovo_ref.good --relabel "U_" --output relabeled_denovo_ref.good && \
        $(CDHIT)/cd-hit-est -i relabeled_denovo_ref.good -o relabeled_denovo_ref_good_cdhit_97 -c 0.97 -M 200000 -T 16
                
renaming:
        cd $(DIR)/2_quality_check/quality_trim && \
        python $(CODE)/renaming_seq_w_short_sample_name.py "S_" $(DIR)/3_cdhit_clustering/renamed_seqs/sample_filename_map.txt $(DIR)/3_cdhit_clustering/renamed_seqs/sequence_name_map.txt *.fa > $(DIR)/3_cdhit_clustering/renamed_seqs/all_renamed_sequences.fa

otu_mapping:
        cd $(DIR)/3_cdhit_clustering/master_otus && \
        $(CDHIT)/cd-hit-est-2d -i relabeled_denovo_ref_good_cdhit_97 -i2 ../renamed_seqs/all_renamed_sequences.fa -o otu_mapping_cdhit97 -c 0.97 -M 200000 -T 16

otu_table:
        cd $(DIR)/3_cdhit_clustering/master_otus && \
        python $(CODE)/cdhit_otu_mapping.py otu_mapping_cdhit97.clstr > cdhit_otu_table_long.txt && \
        Rscript $(CODE)/convert_otu_table_long_to_wide_format.R cdhit_otu_table_long.txt ../R/cdhit_otu_table_wide.txt && \
        python $(CODE)/rep_seq_to_otu_mapping.py relabeled_denovo_ref_good_cdhit_97.clstr > rep_seq_to_cluster.map 

taxa_table:
        cd $(DIR)/3_cdhit_clustering/master_otus && \
        java -Xmx24g -jar $(RDP)/classifier.jar classify -c 0.5 -f filterbyconf -o relabeled_denovo_ref_good_cdhit_97_taxa_filterbyconf.txt -h relabeled_denovo_ref_good_cdhit_97_taxa_filterbyconf_hierarchy.txt relabeled_denovo_ref_good_cdhit_97 && \
        Rscript $(CODE)/renaming_taxa_rep_seq_to_otus.R relabeled_denovo_ref_good_cdhit_97_taxa_filterbyconf.txt rep_seq_to_cluster.map ../R/cdhit_taxa_table_w_repseq.txt
        

clean:
        cd $(DIR) && rm *.sh *.txt && rm -r 1_rdp_pandaseq 2_quality_check 3_cdhit_clustering