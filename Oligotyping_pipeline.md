# Run oligotyping analysis (Oligotyping v. 2.1)

##  Archaea 
```
#pad sequences that differ in length
o-pad-with-gaps archaea_aligned_seqs_final.fna -o arach_padded.fna
#run entropy
entropy-analysis arach_padded.fna --quickmkdir  -o arch_quick_oligo
  #see 3 peaks >0.2 Entropy
#do oligotyping
oligotype Desktop/arach_padded.fna Desktop/arach_padded.fna-ENTROPY -C 3 -o Desktop/arch_oligotyping
```
There's three oligotypes, but one is completly dominant relative to the other...no patterns for the archaea. :(


## Bacillus 
```
#pad sequences that differ in length
o-pad-with-gaps bacillus_aligned_seqs_final.fna -o bacillus_padded.fna

#run entropy
entropy-analysis bacillus_padded.fna --quick
	#see no meaingful peaks > 0.2 entropy

#do oligotyping

oligotype bacillus_padded.fna bacillus_padded.fna-ENTROPY -c 3 -o bacillus_oligotyping

```

## Bacteroidetes 

```
#pad sequences that differ in length
o-pad-with-gaps bacteroidetes_aligned_seqs_final.fna -o bacteroidetes_padded.fna

#run entropy
entropy-analysis bacteroidetes_padded.fna --quick
	#see 12 peaks >0.2 Entropy

#do oligotyping
oligotype bacteroidetes_padded.fna bacteroidetes_padded.fna-ENTROPY -c 12 -o bacteroidetes_oligotyping
```

## Cyanobacteria

```
#pad sequences that differ in length
o-pad-with-gaps cyano_aligned_seqs_final.fna -o cyano_padded.fna

#run entropy
entropy-analysis cyano_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/cyano_padded.fna Desktop/cyano_padded.fna-ENTROPY -C 3 -o Desktop/cyano_oligotyping
```

## Chloroflexi
```
#pad sequences that differ in length
o-pad-with-gaps chloroflexi_aligned_seqs_final.fna -o chloroflexi_padded.fna

#run entropy
entropy-analysis chloroflexi_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/chloroflexi_padded.fna Desktop/chloroflexi_padded.fna-ENTROPY -C 3 -o Desktop/chloroflexi_oligotyping
```

## Verrucomicrobia
```
#pad sequences that differ in length
o-pad-with-gaps verruco_aligned_seqs_final.fna -o verruco_padded.fna

#run entropy
entropy-analysis verruco_padded.fna --quick
	#see 100 peaks >0.2 Entropy

#do oligotyping

oligotype verruco_padded.fna verruco_padded.fna-ENTROPY -c 100 -o verruco_oligotyping -A 10
```

## Acidobacteria
```
#pad sequences that differ in length
o-pad-with-gaps acido_aligned_seqs_final.fna -o acido_padded.fna

#run entropy
entropy-analysis acido_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/acido_padded.fna Desktop/acido_padded.fna-ENTROPY -C 3 -o Desktop/acido_oligotyping
```

## Actinobacteria
```
#pad sequences that differ in length
o-pad-with-gaps actino_aligned_seqs_final.fna -o actino_padded.fna

#run entropy
entropy-analysis actino_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/actino_padded.fna Desktop/actino_padded.fna-ENTROPY -C 3 -o Desktop/actino_oligotyping
```

## Betaproteobacteria
```
#pad sequences that differ in length
o-pad-with-gaps beta_aligned_seqs_final.fna -o beta_padded.fna

#run entropy
entropy-analysis beta_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/beta_padded.fna Desktop/beta_padded.fna-ENTROPY -C 3 -o Desktop/beta_oligotyping
```

## Alphaproteobacteria
```
#pad sequences that differ in length
o-pad-with-gaps alpha_aligned_seqs_final.fna -o alpha_padded.fna

#run entropy
entropy-analysis alpha_padded.fna --quick
	#see 54 peaks >0.2 Entropy

#do oligotyping

oligotype alpha_padded.fna alpha_padded.fna-ENTROPY -c 54 -o alpha_oligotyping -A 10
```

## Gammaproteobacteria
```
#pad sequences that differ in length
o-pad-with-gaps gamma_aligned_seqs_final.fna -o gamma_padded.fna

#run entropy
entropy-analysis gamma_padded.fna --quick
	#see 59! peaks >0.2 Entropy

#do oligotyping

oligotype gamma_padded.fna gamma_padded.fna-ENTROPY -c 59 -o gamma_oligotyping -A 10
```

## Deltaproteobacteria
```
#pad sequences that differ in length
o-pad-with-gaps delta_aligned_seqs_final.fna -o delta_padded.fna

#run entropy
entropy-analysis delta_padded.fna --quick
	#see 3 peaks >0.2 Entropy

#do oligotyping

oligotype Desktop/delta_padded.fna Desktop/delta_padded.fna-ENTROPY -C 3 -o Desktop/delta_oligotyping
```
