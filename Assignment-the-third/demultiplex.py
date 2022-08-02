import numpy as np
import bioinfo as b
import matplotlib.pyplot as plt
import itertools as it
import argparse, gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-1", "--read1", help="the first file for processing. this is read1 file", type=str, required=False)
    parser.add_argument("-2", "--index1", help="the second file for processing. this is index1 file", type=str, required=False)
    parser.add_argument("-3", "--index2", help="the third file for processing. this is index2 file.", type=str, required=False)
    parser.add_argument("-4", "--read2", help="the fourth file for processing. this is read2 file.", type=str, required=False)
    parser.add_argument("-t", "--trialnum", help="trial number for running python script", type=int, required=False)
    return parser.parse_args()

args = get_args()
read1 = args.read1
index1 = args.index1
index2 = args.index2
read2 = args.read2
trial = args.trialnum

r1_rec = np.array(4)
i1_rec = np.array(4)
i2_rec = np.array(4)
r2_rec = np.array(4)

with (
    open(read1, "r") as r1,
    open(index1, "r") as i1,
    open(index2, "r") as i2,
    open(read2, "r") as r2
    ):
    








