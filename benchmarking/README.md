## Benchmarking mybwa to Exisitng BWA index and BWA mem 

## Timing 
To get the timing for mybwa: 
```
time python mybwa.py <genome_fasta> <reads_fastq>
```
To get the timing for existing bwa: 
```
time bwa index <genome_fasta>
time bwa mem  <genome_fasta> <reads_fastq>
```
We have to run to commands for the exisitng bwa as mybwa combines the bwt and aligning of the reads in one command, while the exisiting bwa does not. Hence, we would have to add the resulting time from the 2 commands for the exisitng bwa to compare it to mybwa.

### Example Comparison: 
```
time python mybwa.py hg19chr1.fa reads50.fq
```
output: 
```
real    0m6.707s
user    0m6.617s
sys   0m0.071s
```
```
time bwa index hg19chr1.fa
time bwa mem hg19chr1.fa reads50.fq
```
output: 
```
real    0m0.025s
user    0m0.010s
sys     0m0.008s

real    0m0.009s
user    0m0.004s
sys     0m0.004s
```
## Comparison of Aligned Reads
mybwa command: 
```
python mybwa.py <genome_fasta> <reads_fastq>
```
bwa command: 
```
time bwa index <genome_fasta>
time bwa mem  <genome_fasta> <reads_fastq> > output.sam
```

### Example Benchmarking
```
python mybwa.py hg19chr1.fa reads50.fq

output: 
0       CTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       245
1       CAGTGCACATATACTTGTTTCCCAACCTATTCTCAACTAAAGCCGATTGA
2       CGGTGTGGCGTGGTGCCCGTTACCCGGCTGCAGGTTCACAGAAATCTCAC
3       CATAATTATACCGGCCGTCACAGCGTCGTAATTCCATAATAATAACCCGC
4       GTTGTGGCAGGAGGTGCCGCATCTCCAACAAGGTCGAAGTCGCAAAAGAC
5       TGCCTTCGGTCGAGGGTGGGGGGACCCACTAAAGTGTCGAGTAGCCACTA

```
```
bwa index hg19chr1.fa
bwa mem hg19chr1.fa reads50.fq

output: 
M::process] read 6 sequences (300 bp)...
[M::mem_process_seqs] Processed 6 reads in 0.002 CPU sec, 0.001 real sec
0       0       hg19_chr1_fragment      246     50      50M     *       0       0       CTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC CTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC      NM:i:0  MD:Z:50 AS:i:50 XS:i:36
1       4       *       0       0       *       *       0       0       CAGTGCACATATACTTGTTTCCCAACCTATTCTCAACTAAAGCCGATTGACTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       AS:i:0  XS:i:0
2       4       *       0       0       *       *       0       0       CGGTGTGGCGTGGTGCCCGTTACCCGGCTGCAGGTTCACAGAAATCTCACCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       AS:i:0  XS:i:0
3       4       *       0       0       *       *       0       0       CATAATTATACCGGCCGTCACAGCGTCGTAATTCCATAATAATAACCCGCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       AS:i:0  XS:i:0
4       4       *       0       0       *       *       0       0       GTTGTGGCAGGAGGTGCCGCATCTCCAACAAGGTCGAAGTCGCAAAAGACCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       AS:i:0  XS:i:0
5       4       *       0       0       *       *       0       0       TGCCTTCGGTCGAGGGTGGGGGGACCCACTAAAGTGTCGAGTAGCCACTACTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       AS:i:0  XS:i:0

```
As we can see, mybwa aligns the reads to position 245 and bwa aligns it to 246. This is because the implementation of mybwa uses 0 indexing while the exisiting bwa uses one indexing. 






