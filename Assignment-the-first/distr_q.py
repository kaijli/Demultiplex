
from statistics import mean
import numpy as np
import bioinfo as b
import argparse, gzip

def get_args():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-r", "--reads", help="read length of each fq entry", type=int, required=False)

	return parser.parse_args()

args = get_args()
reads = args.reads
# reads = 101

read1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
index1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
index2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
read2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

file_dict = {read1:101, index1:8, index2:8, read2:101}

# def rev_comp(seq: str, RNAflag: bool = False):
#     """
#     Takes a sequence and returns its reverse compliment.
#     Not case sensitive.
#     """
#     reverse_complement = ""
#     seq = seq.upper()
#     seq = seq.strip()
#     if b.validate_base_seq(seq, RNAflag):
#         if not RNAflag:
#             reverse_complement = "".join(b.DNAcomplement.get(base, base) for base in reversed(seq))
#         else:
#             reverse_complement = "".join(b.RNAcomplement.get(base, base) for base in reversed(seq))
#     return reverse_complement

# print(rev_comp("ACTGGTCACA"))
# print(rev_comp("asfkjhdj"))
# print(rev_comp("ACUGGUCACA", True))
# print(rev_comp("actgg"))

def mean_qscores(file, reads):
    '''
    takes a zipped fq file and returns the mean qscore for each 
    read position in a numpy array. 
    '''
    q_totals = np.zeros(reads, dtype=np.float64)
    # with gzip.open(file,"r") as fh: 
    with open(file,"r") as fh: 
        file_line = 0
        for line in fh:
            file_line +=1
            # line = line.decode()
            line = line.strip()
            if file_line%4==0:
                lin_arr = np.array(list(line))
                print(lin_arr)
                # q_scores += b.convert_phred(line)
                print(q_totals)
                for index in range(len(lin_arr)):
                    q_totals[index] += b.convert_phred(lin_arr[index])
    # print(q_totals)
    # print(q_totals/reads)
    return q_totals/reads


# mean_qscores("/Users/kaitlynli/Documents/UFO/bioinformatics/Bi621/PS/ps4-kaijli/test.fastq")

for key,value in file_dict.items():
    arr_name = str(key)+"mean"
    arr_name = mean_qscores(key, value)
    print(f"{arr_name=}")