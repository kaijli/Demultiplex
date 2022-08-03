import numpy as np
import bioinfo as b
import matplotlib.pyplot as plt
import itertools as it
import argparse, gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-k", "--knownindexes", help="file with known indexes in format sample ID, index columns", type=str, required=False)
    parser.add_argument("-1", "--read1", help="the first file for processing. this is read1 file", type=str, required=False)
    parser.add_argument("-2", "--index1", help="the second file for processing. this is index1 file", type=str, required=False)
    parser.add_argument("-3", "--index2", help="the third file for processing. this is index2 file.", type=str, required=False)
    parser.add_argument("-4", "--read2", help="the fourth file for processing. this is read2 file.", type=str, required=False)
    parser.add_argument("-t", "--trialnum", help="trial number for running python script", type=int, required=False)
    parser.add_argument("-c", "--cutoff", help="cutoff for quality score", type=int, required=False)
    return parser.parse_args()

args = get_args()
known = args.knownindexes
read1 = args.read1
index1 = args.index1
index2 = args.index2
read2 = args.read2
trial = args.trialnum
cutoff = args.cutoff

index1 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/i1_unit_test.fq"
index2 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/i2_unit_test.fq"
read1 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/r1_unit_test.fq"
read2 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/r2_unit_test.fq"

'''
define variables
'''
r1_rec = np.array(4)
i1_rec = np.array(4)
i2_rec = np.array(4)
r2_rec = np.array(4)
index_IDs = dict()
indexes = set()
hopped, unknown = 0,0
matched = dict()
# read_files = {read1:"r1", read2:"r2", index1:"i1", index2:"i2"}
write_files = dict()

# collect indexes into 
known = "indexes.txt"
with open (known, "r") as fh:
    for line in fh:
        line = line.strip()
        (id, ind)=line.split()
        index_IDs[ind] = id
        matched[ind] = 0
# print(sampleIDs)
for key in index_IDs:
    indexes.add(key)
print(indexes)



'''
open 4 files in parallel
need to open all files in parallel
'''
# dict_variable = {key:value for (key,value) in dictionary.items()}
write_files = {index:(open(f"fq_out/r1_{sampleID}_{index}.fq", "w"), open(f"fq_out/r2_{sampleID}_{index}.fq", "w")) for (index, sampleID) in index_IDs.items()}
write_files["hopped"] = (open(f"fq_out/r1_hopped.fq", "w"), open(f"fq_out/r2_hopped.fq", "w"))
write_files["unknown"] = (open(f"fq_out/r1_unknown.fq", "w"), open(f"fq_out/r2_unknown.fq", "w"))

def write_fqs(category:str, r1_rec: np.ndarray, r2_rec: np.ndarray):
    '''
    this function takes category of index, hopped, or unknown and writes to opened fq files
    according to which read and which category with amended header line
    '''
    write_files[category][0].write(f"{r1_rec[0]} {r1_rec[1]}:{b.rev_comp(r2_rec[1])}\n")
    write_files[category][1].write(f"{r2_rec[0]} {r1_rec[1]}:{b.rev_comp(r2_rec[1])}\n")
    for line in r1_rec[1:]:
        write_files[category][0].write(f"{line}\n")
    for line in r2_rec[1:]:
        write_files[category][1].write(f"{line}\n")

def collect_record(file):
    # def next_n_lines(file_opened, N):
    return np.array([x.strip() for x in it.islice(file, 4)])

with (
    open(read1, "r") as r1,
    open(index1, "r") as i1,
    open(index2, "r") as i2,
    open(read2, "r") as r2
    ):
    while True:
        '''
        collect 4 lines at a time into an array
        '''
        hopped, unknown = 0,0
        matched = dict()
        r1_rec = collect_record(r1)
        r2_rec = collect_record(r2)
        i1_rec = collect_record(i1)
        i2_rec = collect_record(i2)
        print(i2_rec)
        if np.size(r1_rec)==0:
            break
            
        '''
        begin processing indexes
        '''
        i2_rev = b.rev_comp(i2_rec[1])
        # are indexes known?
        if not b.check_index(i1_rec[1], indexes) or not b.check_index(i2_rec[1], indexes):
            unknown += 1
            write_fqs("unknown", r1_rec, r2_rec)
            
        # are q score above cutoff
        elif not b.check_qscores(i1_rec[3], cutoff) or not b.check_qscores(i2_rec[3], cutoff):
            unknown += 1
            write_fqs("unknown", r1_rec, r2_rec)
            
        elif i1_rec[1] != i2_rev:
            hopped += 1
            write_fqs("hopped", r1_rec, r2_rec)
            
        elif i1_rec[1] == i2_rev:
            matched[i1_rec[1]] = 1
            write_fqs(i1_rec[1], r1_rec, r2_rec)
            
'''
close all files. 
'''
for file_tuple in write_files.values():
    file_tuple[0].close()
    file_tuple[1].close()