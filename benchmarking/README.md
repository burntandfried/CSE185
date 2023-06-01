## Benchmarking mybwa to exisitng bwa 

## Timing 
To get the timing for mybwa: 
```
time python mybwa.py <genome_fasta> <reads_fastq>
```
To get the timing for existing bwa: 
```
time bwa index <genome_fasta>
time bwa mem  <genome_fasta> <reads_fastq> > output.sam
```
We have to run to commands for the exisitng bwa as mybwa combines the bwt and aligning of the reads in one command, while the exisiting bwa does not. Hence, we would have to add the resulting time from the 2 commands for the exisitng bwa to compare it to mybwa.

### Example comparison: 
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
time bwa mem hg19chr1.fa reads50.fq > output.sam
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
## Comparison of aligned reads
mybwa command: 
```
python mybwa.py <genome_fasta> <reads_fastq>
```
bwa command: 
```
time bwa index <genome_fasta>
time bwa mem  <genome_fasta> <reads_fastq> > output.sam
```

### example 
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
bwa mem hg19chr1.fa reads50.fq > output.sam

output in sam file: 

```






