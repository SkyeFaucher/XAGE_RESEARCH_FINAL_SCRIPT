# Using Bedtools at the linux command line to determine if the promoter regions of the differentially expressed genes align with ChipSeq data of FOXP4

# Establishing working directory and observing the files.
cd /Users/Skye/Desktop/New_Xage/Second_Analysis
ls

# Intersecting the human genome promoter file to the FOXP4 ChipSeq data to determine the regions of the human genome where FOXP4 is binding and saving the result in a BED file for use in further analysis.
bedtools intersect -a human_genome_promoter_region.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed > foxp4_binding_sites_in_human_geneome.bed


# Making the differentially expressed genes promoter region files compatible with the FOXP4 binding sites file.
awk '{print "chr"$0}' xage2_72_upregulated_promoter_region.bed > xage2_72_upregulated_promoter_final.bed
awk '{print "chr"$0}' xage2_72_downregulated_promoter_region.bed > xage2_72_downregulated_promoter_final.bed
awk '{print "chr"$0}' xage2_120_upregulated_promoter_region.bed > xage2_120_upregulated_promoter_final.bed
awk '{print "chr"$0}' xage2_120_downregulated_promoter_region.bed > xage2_120_downregulated_promoter_final.bed


# DE GENES BOUND BY FOXP4
# Intersecting the FOXP4 binding sites and the differentially expressed genes to get the differentially expressed genes that are bound to FOXP4.
bedtools intersect -wa -a xage2_72_upregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_72_upregulated_foxp4_bound_transcripts.txt 
bedtools intersect -wa -a xage2_72_downregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_72_downregulated_foxp4_bound_transcripts.txt
bedtools intersect -wa -a xage2_120_upregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_120_upregulated_foxp4_bound_transcripts.txt
bedtools intersect -wa -a xage2_120_downregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_120_downregulated_foxp4_bound_transcripts.txt


# DE GENES NOT BOUND BY FOXP4
# Intersecting the FOXP4 binding sites and the differentially expressed genes to get the differentially expressed genes that are not bound to FOXP4.
bedtools intersect -v -a xage2_72_upregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_72_upregulated_foxp4_unbound_transcripts.txt
bedtools intersect -v -a xage2_72_downregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_72_downregulated_foxp4_unbound_transcripts.txt
bedtools intersect -v -a xage2_120_upregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_120_upregulated_foxp4_unbound_transcripts.txt
bedtools intersect -v -a xage2_120_downregulated_promoter_final.bed -b foxp4_binding_sites_in_human_geneome.bed | sort | uniq > xage2_120_downregulated_foxp4_unbound_transcripts.txt


# NON-DE GENES NOT BOUND BY FOXP4
# Intersecting the human genome promoter file to the FOXP4 ChipSeq data to determine the regions of the human genome where FOXP4 is not binding.
bedtools intersect -v -a human_genome_promoter_region.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed > foxp4_nonbinding_sites_in_human_genome.bed

# Intersecting the resulting BED file with the differentially expressed genes to determine the genes that are not differentially expressed and not bound by FOXP4.
bedtools intersect -v -a foxp4_nonbinding_sites_in_human_genome.bed -b xage2_72_upregulated_promoter_final.bed | sort | uniq > xage2_72_upregulated_nonbinding_nonDE_transcripts.txt
bedtools intersect -v -a foxp4_nonbinding_sites_in_human_genome.bed -b xage2_72_downregulated_promoter_final.bed | sort | uniq > xage2_72_downregulated_nonbinding_nonDE_transcripts.txt
bedtools intersect -v -a foxp4_nonbinding_sites_in_human_genome.bed -b xage2_120_upregulated_promoter_final.bed | sort | uniq > xage2_120_upregulated_nonbinding_nonDE_transcripts.txt
bedtools intersect -v -a foxp4_nonbinding_sites_in_human_genome.bed -b xage2_120_downregulated_promoter_final.bed | sort | uniq > xage2_120_downregulated_nonbinding_nonDE_transcripts.txt


# NONDE GENES BOUND BY FOXP4
# Intersecting the human genome promoter file to the differentially expressed genes files to determine the genes that are not differentially expressed.
bedtools intersect -v -a human_genome_promoter_region.bed -b xage2_72_upregulated_promoter_final.bed > xage2_72_upregulated_nonDE.bed 
bedtools intersect -v -a human_genome_promoter_region.bed -b xage2_72_downregulated_promoter_final.bed > xage2_72_downregulated_nonDE.bed
bedtools intersect -v -a human_genome_promoter_region.bed -b xage2_120_upregulated_promoter_final.bed > xage2_120_upregulated_nonDE.bed
bedtools intersect -v -a human_genome_promoter_region.bed -b xage2_120_downregulated_promoter_final.bed > xage2_120_downregulated_nonDE.bed

# Intersecting the resulting files with the FOXP4 ChipSeq data to determine the genes that are not differentially expressed that are bound by FOXP4.
bedtools intersect -wa -a xage2_72_upregulated_nonDE.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed | sort | uniq > xage2_72_upregulated_nonDE_foxp4_bound_transcripts.txt
bedtools intersect -wa -a xage2_72_downregulated_nonDE.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed | sort | uniq > xage2_72_downregulated_nonDE_foxp4_bound_transcripts.txt
bedtools intersect -wa -a xage2_120_upregulated_nonDE.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed | sort | uniq > xage2_120_upregulated_nonDE_foxp4_bound_transcripts.txt
bedtools intersect -wa -a xage2_120_downregulated_nonDE.bed -b GSM5214071_ENCFF433ACM_IDR_ranked_peaks_GRCh38.bed | sort | uniq > xage2_120_downregulated_nonDE_foxp4_bound_transcripts.txt

