# Analysis of 16S rRNA gene Illumina sequencing data

All analyses were performed in USEARCH v. 10.0.240 (64-bit) with the UPARSE OTU picking method. Subsequent analyses were performed in QIIME v. 1.8. 

## Merge Paired End Reads
```
#decompress the reads
gunzip *.gz

mkdir mergedfastq

./usearch64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 300
```

## Dereplicate sequences
```
./usearch64 -fastx_uniques mergedfastq/merged.fq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout
```

## Remove Singeltons
```
./usearch64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster Sequences
```
./usearch64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Reference-based OTU picking against the 123 version of Silva
```
./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db /mnt/home/kearnspa/silva_db/silva123.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```

## Sort by size and then de novo OTU picking on sequences that failed to hit Silva
```
./usearch64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

./usearch64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom
```


# Switch to QIIME

## Assign taxonomy to Silva v.123 and UCLUST
```
assign_taxonomy.py -i full_rep_set.fna -o taxonomy -r /mnt/home/kearnspa/silva_db/silva123.fasta -t /mnt/home/kearnspa/silva_db/silva123_taxonomy.txt
```

## Add taxonomy to OTU table
```
biom add-metadata -i OTU_jsn.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```

## Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt.biom -n o__Streptophyta,o__Chlorophyta,f__mitochondria,Unassigned
```

## Align sequences to Silva with PyNast
```
align_seqs.py -i full_rep_set.fna -o alignment -t /mnt/home/kearnspa/silva_db/silva123.fasta
```

## Filter excess gaps from alignment
```
filter_alignment.py -i alignment/full_rep_set_aligned.fasta -o alignment/filtered_alignment
```

## Make phylogeny with fasttree
```
make_phylogeny.py -i alignment/filtered_alignment/full_rep_set_aligned_pfiltered.fasta -o rep_set.tre
```

## Summarize the OTU table and rarefy OTU table to lowest sequencing depth
```
biom summarize-table -i otu_table_tax_filt.biom -o otu_table_summary.txt

single_rarefaction.py -d 16618 -o single_rare.biom -i otu_table_tax_filt2.biom
```

## Calculate alpha and beta diversity
```
beta_diversity.py -m bray_curtis,unweighted_unifrac,weighted_unifrac -i single_rare.biom -o beta_div -t rep_set.tre
principal_coordinates.py -i beta_div -o coords

alpha_diversity.py -m PD_whole_tree,shannon -i single_rare.biom -o alpha -t rep_set.tre
```

## Summarize taxonomic data
```
summarize_taxa.py -i otu_table_tax_filt.biom -o taxa_sum
```

## Adonis
```
#make site-specific OTU tables by filtering out other site.

filter_samples_from_otu_table.py -i single_rare.biom -o coastal_table.biom -m mimulus_map.txt -e -s 'Site_planted:Coastal'
filter_samples_from_otu_table.py -i single_rare.biom -o inland_table.biom -m mimulus_map.txt -e -s 'Site_planted:Inland'


#global adonis (Anderson 2001, Austral Ecol)
compare_categories.py -n 10000 -i single_rare.biom -o total_genotype_adonis --method=adonis -m mimulus_map.txt -c 'Genotype'
#F=46.213, d.f.=3,37, p<0.0001, R2=0.33

compare_categories.py -n 10000 -i single_rare.biom -o total_site_adonis --method=adonis -m mimulus_map.txt -c 'Site_planted'
#F=111.34, d.f.=1,39, p<0.0001, R2=0.72

#site/genotype adonis (Anderson 2001, Austral Ecol)
compare_categories.py -n 10000 -i coastal_table.biom -m mimulus_map.txt -e -s 'Genotype' -o coastal_genotype_adonis
#F=32.653, d.f.=3,17, p<0.0001, R2=0.67

compare_categories.py -n 10000 -i inland_table.biom -m mimulus_map.txt -e -s 'Genotype' -o inland_adonis
#F=36.21, d.f.=3,17, p<0.0001, R2=0.59

#GT Origin within a site
compare_categories.py -n 10000 -i coastal_table.biom -m mimulus_map.txt -e -s 'Origin' -o coastal_genotype_adonis_origin
#F=41.24, d.f.=1,19, p<0.0001, R2=0.45

compare_categories.py -n 10000 -i inland_table.biom -m mimulus_map.txt -e -s 'Origin' -o inland_adonis_origin
#F=24.312, d.f.=1,19, p<0.0001, R2=0.51
```

## Compare community similarity between sites and GT within/between sites
```
#between sites weighted unifrac
make_distance_boxplots.py -m mimulus_alphaphied_map.txt -f 'Site_planted' -o site_boxplots --save_raw_data -d beta_div/weighted_unifrac_otu_table_tax_filt.txt

#GT between sites
make_distance_boxplots.py -m mimulus_alphaphied_map.txt -f 'Site_GT' -o site_gt_boxplots --save_raw_data -d beta_div/weighted_unifrac_otu_table_tax_filt.txt
```

## Kruskal-wallis test for determining the differences between GT at a site
```
group_significance.py -i coastal_otu_table.biom -m mimulus_alphaphied_map.txt -c Genotype -o coastal_gt_kruskal
group_significance.py -i inland_otu_table.biom -m mimulus_alphaphied_map.txt -c Genotype -o inland_gt_kruskal
```

## Comparing diversity in rare v. abundant taxa
```
#1%=270 reads/sample
#OTU table with top taxa (n=800 OTUs)
filter_otus_from_otu_table.py -i single_rare.biom -o top_taxa_rare.biom -n 270 

#convert to .txt file
 biom convert -i top_taxa_rare.biom -o top_taxa_rare.txt -b
 
 #grab OTU IDs from this file
 awk -F"\t" '{if ($1) print $1}' top_taxa_rare.txt > top_taxa.txt

#use nano to remove first 2 lines 

#OTU table with only rare taxa
 filter_otus_from_otu_table.py -i single_rare.biom -o rare_taxa_rare.biom -e top_taxa.txt
 
 #calculate PD/Shannon/OTUs obs
 alpha_diversity.py -i rare_taxa_rare.biom -o rare_diversity.txt -t rep_set.tre -m shannon,observed_species,PD_whole_tree
 alpha_diversity.py -i top_taxa_rare.biom -o top_diversity.txt -t rep_set.tre -m shannon,observed_species,PD_whole_tree
 
 #add alpha diversity to mapping file
  add_alpha_to_mapping_file.py -i rare_diversity.txt -m mimulus_map.txt -o rare_taxa_map.txt
add_alpha_to_mapping_file.py -i top_diversity.txt -m mimulus_map.txt -o top_taxa_map.txt
```
 