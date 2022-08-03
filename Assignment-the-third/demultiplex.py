from typing import Any
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
    parser.add_argument("-c", "--cutoff", help="cutoff for quality score", type=int, required=False, default = "30")
    parser.add_argument("-s", "--cutoffstyle", help="cutoff style, avg or indv q scores", type=str, required=False, default = "avg")
    return parser.parse_args()

args = get_args()
known = args.knownindexes
read1 = args.read1
index1 = args.index1
index2 = args.index2
read2 = args.read2
trial = args.trialnum
cutoff = args.cutoff
style = args.cutoffstyle

'''
test files
'''
# index1 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/i1_unit_test.fq"
# index2 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/i2_unit_test.fq"
# read1 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/r1_unit_test.fq"
# read2 = "/projects/bgmp/kli8/bioinformatics/Bi622/Demultiplex/Assignment-the-first/unit_tests/r2_unit_test.fq"
# known = "indexes.txt"

'''
define and set variables
'''
r1_rec = np.array(4)
i1_rec = np.array(4)
i2_rec = np.array(4)
r2_rec = np.array(4)
index_IDs = dict()
indexes = set()
hopped_indexes, hopped_pair, matched_indexes = dict(),dict(),dict()
write_files = dict()
unknown, hopped, matched = 0,0,0

'''
collect indexes from txt file
'''
with open (known, "r") as fh:
    for line in fh:
        line = line.strip()
        (id, ind)=line.split()
        index_IDs[ind] = id
        matched_indexes[ind],hopped_indexes[ind]  = 0,0
# print(sampleIDs)
for key in index_IDs:
    indexes.add(key)
print(indexes)


'''
functions to use while demultiplexing
'''
# dict_variable = {key:value for (key,value) in dictionary.items()}
def write_fqs(category:str, r1_rec: np.ndarray, r2_rec: np.ndarray):
    '''
    this function takes category of index, hopped, or unknown and writes to opened fq files
    according to which read and which category with amended header line
    '''
    write_files[category][0].write(f"{r1_rec[0]} {i1_rec[1]}:{b.rev_comp(i2_rec[1])}\n")
    write_files[category][1].write(f"{r2_rec[0]} {i1_rec[1]}:{b.rev_comp(i2_rec[1])}\n")
    for line in r1_rec[1:]:
        write_files[category][0].write(f"{line}\n")
    for line in r2_rec[1:]:
        write_files[category][1].write(f"{line}\n")

def collect_record(file):
    '''
    Collects 4 lines at a time from fq files to retain current record
    '''
    return np.array([x.strip() for x in it.islice(file, 4)])

def hopped_counter(indexpair: tuple=("unknown", "unknown")):
    '''
    counts the different hopped pairings
    '''
    if indexpair not in hopped_indexes.keys():
        hopped_pair[indexpair] = 1
    else:
        hopped_pair[indexpair] += 1

'''
open all writing files
'''
write_files = {index:(open(f"fq_out_{trial}/r1_{sampleID}_{index}.fq", "w"), open(f"fq_out_{trial}/r2_{sampleID}_{index}.fq", "w")) for (index, sampleID) in index_IDs.items()}
write_files["hopped"] = (open(f"fq_out_{trial}/r1_hopped.fq", "w"), open(f"fq_out_{trial}/r2_hopped.fq", "w"))
write_files["unknown"] = (open(f"fq_out_{trial}/r1_unknown.fq", "w"), open(f"fq_out_{trial}/r2_unknown.fq", "w"))


'''
open 4 input files for reading in parallel
'''
with ( #for non test files use gzip.open
    gzip.open(read1, "rt") as r1,
    gzip.open(index1, "rt") as i1,
    gzip.open(index2, "rt") as i2,
    gzip.open(read2, "rt") as r2
    # open(read1, "rt") as r1,
    # open(index1, "rt") as i1,
    # open(index2, "rt") as i2,
    # open(read2, "rt") as r2
    ):
    # print(f"{unknown=}\n{hopped=}\n{matched=}")
    print("Files opened")
    while True:
        '''
        collect 4 lines at a time into an array
        '''
        r1_rec = collect_record(r1)
        r2_rec = collect_record(r2)
        i1_rec = collect_record(i1)
        i2_rec = collect_record(i2)
        # print(i2_rec)
        if np.size(r1_rec)==0:
            # print(f"{unknown=}\n{hopped=}\n{matched=}")
            break
            
        '''
        begin processing indexes
        first checks if either index is unknown, 
        then if either index is a low quality read,
        then if the index pair doesn't match,
        then if the index pair matches. 
        '''
        i2_rev = b.rev_comp(i2_rec[1])
        index_pair = (i1_rec[1], i2_rev)
        # b.check_indexes(i1_rec[1], i2_rev, indexes)
        # are indexes known?
        if not b.check_indexes(index_pair, indexes):
            # print(i1_rec[1] in indexes)
            unknown += 1
            write_fqs("unknown", r1_rec, r2_rec)
            # print("index fail")
        # are q score above cutoff
        elif not b.check_qscores(i1_rec[3],i2_rec[3], stat=style):
            unknown += 1
            write_fqs("unknown", r1_rec, r2_rec)
            # print("qscore fail")
        elif i1_rec[1] != i2_rev:
            hopped_counter(index_pair)
            hopped_indexes[i1_rec[1]] +=1
            hopped_indexes[i2_rev] +=1
            hopped +=1
            write_fqs("hopped", r1_rec, r2_rec)
            # print("match fail")
        elif i1_rec[1] == i2_rev:
            matched_indexes[i1_rec[1]] +=1
            matched += 1
            write_fqs(i1_rec[1], r1_rec, r2_rec)
            # print("matched")

'''
close all fastq files. 
'''
for file_tuple in write_files.values():
    file_tuple[0].close()
    file_tuple[1].close()

print("Files Closed")

'''
calculations / manipulations for summary information
'''
proportions = dict()
total = hopped + matched + unknown
proportions["hopped"] = hopped/total
proportions["matched"] = matched/total
proportions["unknown"] = unknown/total
sorted_matched = dict(sorted(matched_indexes.items(), key=lambda item: item[1], reverse = True))

'''
function for summary statistics
'''
def write_dict_tsv(data_dict):
    with open(f"fq_out_{trial}/{data_dict}_counts.tsv", "w") as fh:
        for index, count in data_dict.items():
            fh.write(f"{index}\t{count}\n")

'''
write out summary statistics
'''
with open(f"fq_out_{trial}/run_summary.txt", "w") as fh:
    fh.write(f"Demultiplexing Run Number {trial} Summary\n")
    fh.write(f"Read 1 File:\t{read1}\n")
    fh.write(f"Index 1 File:\t{index1}\n")
    fh.write(f"Index 2 File:\t{index2}\n")
    fh.write(f"Read 2 File:\t{read2}\n")
    fh.write(f"Index File:\t{known}\n")
    fh.write(f"QScore Cutoff:\t{cutoff}\n")
    fh.write(f"QScore Cutoff Style:\t{style}\n")
    fh.write(f"Number of Unknown Records:\t{unknown}\n")
    fh.write(f"Number of Hopped Records:\t{hopped}\n")
    fh.write(f"Number of Matched Records:\t{matched}\n")

write_dict_tsv(sorted_matched)
write_dict_tsv(hopped_pair)


'''
plotting various histograms
'''
print("plotting")
# proportion of matched, hopped, and unknown records. 
fig, ax = plt.subplots()
ax.bar(list(proportions.keys()), list(proportions.values()))
plt.xlabel("Index Category")
plt.ylabel("Proportion")
plt.title("Proportion of Demux Categories")
plt.savefig(f"fq_out_{trial}/proportion_hist.png", bbox_inches = "tight")

# the frequency of matched indexes
fig, ax = plt.subplots()
x = list(sorted_matched.keys())
y = list(sorted_matched.values())
ax.bar(x, y)
plt.xticks(rotation = 35, ha = "right", fontsize=9)
plt.xlabel("Index Pairs")
plt.ylabel("Record Count")
plt.title("Frequency of Matched Indexes")
plt.savefig(f"fq_out_{trial}/matched_hist.png", bbox_inches = "tight")

# the frequency of hopped indexes
fig, ax = plt.subplots()
ax.bar(list(hopped_indexes.keys()), list(hopped_indexes.values()))
plt.xticks(rotation = 35, ha = "right", fontsize=9)
plt.xlabel("Hopped Index")
plt.ylabel("Record Count")
plt.title("Frequency of Index Hops")
plt.savefig(f"fq_out_{trial}/hopped_hist.png", bbox_inches = "tight")

# the frequency of hopped index pairings
x =  [str(i) for i in list(hopped_pair.keys())]
fig, ax = plt.subplots()
ax.bar(x, list(hopped_pair.values()))
plt.xticks(rotation = 35, ha = "right", fontsize=9)
plt.xlabel("Hopped Index Pair")
plt.ylabel("Record Count")
plt.title("Frequency of Hopped Index Pairs")
plt.savefig(f"fq_out_{trial}/hopped_pairs_hist.png", bbox_inches = "tight")

print("Finished!")