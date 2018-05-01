# Prepare reads for Oligotyping

This is a workflow for processing sequences from the uparse pipeline.

## Make core OTUs OTU table
```
#make OTU table of only core OTUs
filter_otus_from_otu_table.py -i otu_table_tax_filt.biom -o mim_core_otu_table.biom --negate_ids_to_exclude -e mim_core.txt

#convert to .txt file
biom convert -i mim_core_otu_table.biom -o mim_core_otu_table.txt -b --table-type='OTU table' --header-key=taxonomy

#get rep set for core OTUs
filter_fasta.py -i full_rep_set.fna -o mimulus_core_seqs.fna -b mim_core_otu_table.biom
```

## Extract names of OTUs
```
#grab OTU IDs from fasta file
grep ">" mimulus_core_seqs.fna > mim_core.txt
#remove '>', leaving only name of OTU
sed -e "s/>//" mim_core.txt > mim_core2.txt
```

## Convert fastq reads to fasta
```
#use FASTX (v 0.0.14) and convert fastq of all reads to fasta format
fastq_to_fasta -i merged.fq -o merged.fa
```

## Get OTU map for 'core OTUs' from full OTU map
```
#grab OTUs of interest from large OTU map
grep -f mim_core2.txt OTU_map.uc > core_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' core_otu_map.txt > core_otu_map2.txt

#remove OTU
awk '{print $1}' core_otu_map2.txt > core_otu_map3.txt
```

## Filter non-'core' reads
```
#filter reads not matching OTUs of interest
module load QIIME/1.8.0

filter_fasta.py -s core_otu_map3.txt -f merged.fa -o core_otus.fa

#count number of sequences
grep -c '>' core_otus.fa
263852

#make reads oligotype friendly (change "." to "_")
sed -e 's/\./_/g' core_otus.fa > core_otus_oligo.fa
```

## Align reads to Silva
```
#align sequences (submitted as cluster job b/c ~200k sequences)
align_seqs.py -i core_otus_oligo.fa -o aligned_core_otus_oligo -t silva123.fasta

#filter the alignment
filter_alignment.py -i core_otus_oligo_aligned.fasta -o filtered_alignment
```

# Oligotype analysis of each taxonomic group

Split the  alignment by taxonomic groups (Chloroflexi, Cyanobacteria, Verrucomicrobia, Acidobacteria, Actinobacteria, Betaproteobacteria, Gammaproteobacteria, Deltaproteobacteria, Bacteroidetes, Firmicutes, Alphaproteobacteria, Archaea).


## Process Chloroflexi Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/chloroflexi.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/chloroflexi_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/chloroflexi_otu_map.txt > ./mergedfastq/seqs_by_otu/chloroflexi_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/chloroflexi_otu_map2.txt > ./mergedfastq/seqs_by_otu/chloroflexi_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/chloroflexi_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/chloroflexi_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/chloroflexi_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/chloroflexi_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/chloroflexi_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/chloroflexi_aligned_seqs_final.fna
```

## Process Cyanobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/cyano.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/cyano_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/cyano_otu_map.txt > ./mergedfastq/seqs_by_otu/cyano_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/cyano_otu_map2.txt > ./mergedfastq/seqs_by_otu/cyano_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/cyano_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/cyano_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/cyano_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/cyano_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/cyano_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/cyano_aligned_seqs_final.fna
```

## Process Verrucomicrobia Sequences
```
#grab sequences of interest form large OTU map (n=9 OTUs)
grep '142010\|512563\|533650\|564262\|901183\|1107011\|1107128\|OTU_dn_94' ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/verruco_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/verruco_otu_map.txt > ./mergedfastq/seqs_by_otu/verruco_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/verruco_otu_map2.txt > ./mergedfastq/seqs_by_otu/verruco_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/verruco_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/verruco_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/verruco_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/verruco_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/verruco_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/verruco_aligned_seqs_final.fna
```

## Process Acidobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/acidobacteria.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/acido_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/acido_otu_map.txt > ./mergedfastq/seqs_by_otu/acido_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/acido_otu_map2.txt > ./mergedfastq/seqs_by_otu/acido_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/acido_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/acido_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/acido_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/acido_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/acido_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/acido_aligned_seqs_final.fna
```

## Process Actinobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/actinobacteria.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/actino_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/actino_otu_map.txt > ./mergedfastq/seqs_by_otu/actino_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/actino_otu_map2.txt > ./mergedfastq/seqs_by_otu/actino_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/actino_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/actino_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/actino_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/actino_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/actino_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/actino_aligned_seqs_final.fna
```

## Process Alphaproteobacteria Sequences
```
#grab sequences of interest form large OTU map
grep '303643\|1105814\|331835\|562311\|308836' ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/alpha_otu_map.txt


#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/alpha_otu_map.txt > ./mergedfastq/seqs_by_otu/alpha_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/alpha_otu_map2.txt > ./mergedfastq/seqs_by_otu/alpha_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/alpha_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/alpha_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/alpha_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/alpha_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/alpha_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/alpha_aligned_seqs_final.fna
```

## Process Betaproteobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/beta.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/beta_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/beta_otu_map.txt > ./mergedfastq/seqs_by_otu/beta_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/beta_otu_map2.txt > ./mergedfastq/seqs_by_otu/beta_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/beta_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/beta_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/beta_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/beta_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/beta_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/beta_aligned_seqs_final.fna
```

## Process Gammaproteobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/gamma.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/gamma_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/gamma_otu_map.txt > ./mergedfastq/seqs_by_otu/gamma_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/gamma_otu_map2.txt > ./mergedfastq/seqs_by_otu/gamma_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/gamma_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/gamma_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/gamma_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/gamma_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/gamma_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/gamma_aligned_seqs_final.fna
```

## Process Deltaproteobacteria Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/delta.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/delta_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/delta_otu_map.txt > ./mergedfastq/seqs_by_otu/delta_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/delta_otu_map2.txt > ./mergedfastq/seqs_by_otu/delta_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/delta_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/delta_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/delta_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/delta_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/delta_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/delta_aligned_seqs_final.fna
```

## Process Bacteroidetes Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/bacteroidetes.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map.txt > ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map2.txt > ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/bacteroidetes_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/bacteroidetes_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/bacteroidetes_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/bacteroidetes_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/bacteroidetes_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/bacteroidetes_aligned_seqs_final.fna
```

## Process Firmicutes Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/bacillus.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/bacillus_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/bacillus_otu_map.txt > ./mergedfastq/seqs_by_otu/bacillus_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' ./mergedfastq/seqs_by_otu/bacillus_otu_map2.txt > ./mergedfastq/seqs_by_otu/bacillus_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  1 -f ./mergedfastq/seqs_by_otu/bacillus_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/bacillus_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' ./mergedfastq/seqs_by_otu/bacillus_aligned_seqs.fna > ./mergedfastq/seqs_by_otu/bacillus_aligned_seqs2.fna

sed -e 's/1..253//g' ./mergedfastq/seqs_by_otu/bacillus_aligned_seqs2.fna > ./mergedfastq/seqs_by_otu/bacillus_aligned_seqs_final.fna
```


## Process Archaeal Sequences
```
#grab sequences of interest form large OTU map
grep -f ./mergedfastq/seqs_by_otu/archaea_oligo_otu.txt ./mergedfastq/OTU_map.uc > ./mergedfastq/seqs_by_otu/archaea_otu_map.txt

#remove excess columns keep OTU to read match
awk '{print $9,$10}' ./mergedfastq/seqs_by_otu/archaea_otu_map.txt > archaea_otu_map2.txt

#remove OTUs, leaving only sequences of interest
awk '{print $1}' core_otu_map2.txt > core_otu_map3.txt

#filter aligned fasta file, leaving only archaeal sequencces to oligotype
grep -A  2 -f ./mergedfastq/seqs_by_otu/archaea_otu_map3.txt aligned_core_otus_oligo/filtered_alignment/core_otus_oligo_aligned_pfiltered.fasta > ./mergedfastq/seqs_by_otu/archaea_aligned_seqs.fna


#purge some text from the header 
sed -e '/1..253/s/ *//' archaea_aligned_seqs.fna > archaea_aligned_seqs2.fna
sed -e 's/1..253//g' archaea_aligned_seqs2.fna > archaea_aligned_seqs_final.fna
```
