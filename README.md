# CSE185-BWA-project
This is the project for CSE185, done by Jade Chng, Ananya Prasad and Anna-Sophia Dinov. It implements a Burrows-Wheeler Transfrom to index the reference genome and then performs Burrows Wheeler Alignment to align reads to the genome. 

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
The standard algorithm for BWT alignment uses high and low pointers to search the BWT for reads and uses the suffix array of the genome to determine the positions where any matches occurred. To optimize this exact read matching, we implemented the alignment using a Count dictionary, which for each possible symbol (A, C, G, and T), stores an array whose value at position i is the sum of the number of times that symbol occurs in the first i positions of the BWT. This stores information that is used to update the high and low pointers more efficiently, which lowers the amount of time needed to search for each query, and therefore speeds up the runtime of our alignment algorithm.

# Contributors 

# Testing 

* bonus bandges

