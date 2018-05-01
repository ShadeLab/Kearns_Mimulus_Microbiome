# BLASTn
### version 2.2.27+

```
#blast against the 'nt' database, exclusing environmental sequences
blastn -query mimulus_core_seqs.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out mimulus_oligotypes_blast -negative_gilist sequence.gi -db nt
```
