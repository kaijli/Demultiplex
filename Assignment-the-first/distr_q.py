
import numpy as np
import bioinfo as b
import matplotlib.pyplot as plt
import argparse, gzip


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-r", "--reads", help="read length of each fq entry", type=int, required=False)
    parser.add_argument("-f", "--filename", help="file for processing", type=str, required=False)
    parser.add_argument("-i", "--iteration", help="run number for naming plots", type=int, required=False)
    parser.add_argument("-d", "--directory", help="output directory for plots", type=str, required=False, default=".")
    return parser.parse_args()

args = get_args()
filename = args.filename
reads = args.reads
it = args.iteration
dir = args.directory 
# reads = 101

# read1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
# index1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# index2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# read2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

# file_dict = {read1:101, index1:8, index2:8, read2:101}

def mean_qscores(file, reads):
    '''
    takes a zipped fq file and returns the mean qscore for each 
    read position in a numpy array. 
    '''
    q_totals = np.zeros(reads, dtype=np.float64)
    with gzip.open(file,"r") as fh:  # comment out for non zipped test files.
    # with open(file,"r") as fh: 
        file_line = 0
        for line in fh:
            file_line +=1
            line = line.decode()   # comment out for non zipped test files.
            line = line.strip()
            if file_line%4==0:
                lin_arr = np.array(list(line))
                for index in range(len(lin_arr)):
                    q_totals[index] += b.convert_phred(lin_arr[index])
        q_avgs = q_totals/(file_line/4)
    return q_avgs

def plot_dists(arr):
    '''
    takes an array and plots a distribution plot
    '''
    fig, ax = plt.subplots()
    ax.bar(range(len(arr)), arr)
    plt.xlabel('Read Position')
    plt.ylabel('Mean Quality Score')
    plt.title(f"Mean Quality Scores of Illumina Read Base Pairs\n{filename}")
    plt.savefig(f"{dir}/hist{it}.png")

means = mean_qscores(filename, reads)
plot_dists(means)



