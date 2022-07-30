# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
    A quality score around 30 is a good cutoff for index and biological reads because it's a standard cutoff often seen in papers and studies published in the past several years. This is also based on the histograms, where the lowest average was just above 30. 

    3. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)
    /usr/bin/time -v zcat ($index1 or $index2) | sed -n '2~4p'| grep -c 'N'
    index1: 3976613
    index2: 3328051 
    
## Part 2
1. Define the problem

We are given 4 large fasta files that have some assortment of index pairs and reads. We want to be able to sort the reads by index pairs for data analysis. The index pairs correspond to different experiments. The index1 and index2 files are connected in that the second index is supposed to be the reverse complement of the first. 

2. Describe output

The output of the python script should be 52 different files. There will be 48 files for different index pairs and the two corresponding reads, 2 files for unknown indexes, and 2 files for unmatched indexes. All of these will be very large fastq files. 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
