# CSE185-BWA-project
This is the project for CSE185, done by Jade Chng, Ananya Prasad and Anna-Sophia Dinov. It implements a Burrows-Wheeler Transfrom to index the reference genome and then performs Burrows Wheeler Alignment to align reads to the genome. 

# Install Instructions

# Basic Usage Instructions 
The basic usage of mybwa is: 

```
python mybwa.py <genome_fasta> <reads_fastq>
```

you can run mybwa on the provded fasta and fastq files located in the example_files folder.

for example: 

```
python mybwa.py hg19chr1.fa reads50.fq
```

This should produce to standard output: 

```
0       CTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACC       245
1       CAGTGCACATATACTTGTTTCCCAACCTATTCTCAACTAAAGCCGATTGA
2       CGGTGTGGCGTGGTGCCCGTTACCCGGCTGCAGGTTCACAGAAATCTCAC
3       CATAATTATACCGGCCGTCACAGCGTCGTAATTCCATAATAATAACCCGC
4       GTTGTGGCAGGAGGTGCCGCATCTCCAACAAGGTCGAAGTCGCAAAAGAC
5       TGCCTTCGGTCGAGGGTGGGGGGACCCACTAAAGTGTCGAGTAGCCACTA

```
The first column is the read name, the second is the read and the columns after would be the position/index (0-based) that the match was found. 



# File format 
Genome file: A fasta file with the first line being a header line followed by a single string that represents the database to be queried on the next line.\
Query file: A fastq file where lines with the query sequence are a single string (no new line or space characters)

# Implementation Details
We optimized our Burrow's Wheeler alignment by...

# Contributors 

# Testing 

* bonus bandges

