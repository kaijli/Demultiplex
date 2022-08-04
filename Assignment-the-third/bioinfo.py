#!/user/bin/env python

# Author: Kaitlyn Li <kaitlyn.j.li@gmail.com>
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module
"""This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
"""
__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

import numpy as np
import itertools as it

DNA_bases = set("ACGTNacgtn") #each individual base as capital or lowercase
RNA_bases = set("ACGUNacgun")
DNAcomplement = {"A": "T", "C": "G", "G": "C", "T": "A", "N":"N"} 
RNAcomplement = {"A": "U", "C": "G", "G": "C", "U": "A", "N":"N"} 

print("in my Bioinfo module, __name__ is:", __name__)

def init_arr(size: int) -> np.ndarray:
    """
    initiates a numpy array based off requested size.
    """
    arr = np.zeros(int(size), dtype=np.float64)
    return arr
        
def convert_phred(letter: str, coding: int = 33) -> int:
    """Converts a single character into a phred score"""
    value = ord(letter)-coding
    return value

def qual_score(phred_score: str) -> float:
    """This function takes unmodified phred_score string as parameter and calculates the average quality
    score of the whole phred string"""
    total = 0
    for index in range(len(phred_score)):
        total+=convert_phred(phred_score[index])
    return total/len(phred_score)

def qs_arr(phred_score: str) -> np.ndarray:
    """This function takes unmodified phred_score string as parameter and returns an array
    containing each score position translated to its quality score"""
    arr = np.zeros(int(len(phred_score)), dtype=np.float64)
    for index in range(len(phred_score)):
        arr[index] = convert_phred(phred_score[index])
    return arr

def validate_base_seq(seq: str,RNAflag: bool = False) -> bool:
    """This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive."""
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str) -> float:
    gc_bases = seq.count("C") + seq.count("G")
    return gc_bases/len(seq)

def oneline_fasta(filein: str, fileout: str):
    """
    Thie function takes a fasta file where the sequences have line breaks in them.
    It takes two filenames, and makes a new files with the name input second. 
    """
    collect_seqs = {}
    with open(filein, "r") as fh:
        header = ""
        for line in fh:
            line = line.strip()
            if ">" in line:
                header = line
            else:
                if header not in collect_seqs.keys():
                    collect_seqs[header] = line
                else:
                    collect_seqs[header] += line
    with open(fileout, "w") as fh:
        for header, seq in collect_seqs.items():
            fh.write(f"{header}\n{seq}\n")

def rev_comp(seq: str, RNAflag: bool = False) -> str:
    """
    Takes a sequence and returns its reverse compliment.
    Not case sensitive.
    """
    reverse_complement = ""
    seq = seq.upper()
    seq = seq.strip()
    if validate_base_seq(seq, RNAflag):
        if not RNAflag:
            reverse_complement = "".join(DNAcomplement.get(base, base) for base in reversed(seq))
        else:
            reverse_complement = "".join(RNAcomplement.get(base, base) for base in reversed(seq))
    return reverse_complement

def check_indexes(index_pair: tuple, index_set: set) -> bool:
    """
    takes index1 string and reverse of index2 string
    checks both are in a set of indexes
    """
    known = False
    if index_pair[0] in index_set and index_pair[1] in index_set:
        known = True
    # print(known)
    return known

def check_qscores(phred_score1: str, phred_score2: str, cutoff: int = 30, stat: str = "avg") -> bool:
    """
    takes two phred strings for set of indexes
    returns whether both are above cutoff
    stat indicates whether it"s each qscore or average of string
    indv = checks if any q score in phred string is below cutoff
    avg = checks if average q score of phred string is below cutoff
    cutoff is inclusive
    """
    good = False
    if stat == "indv":
        arr1 = qs_arr(phred_score1)
        arr2 = qs_arr(phred_score2)
        new_arr1 = np.array(it.takewhile(lambda x: x >=cutoff, arr1))
        new_arr2 = np.array(it.takewhile(lambda x: x >=cutoff, arr2))
        if np.array_equal(arr1, new_arr1) and np.array_equal(arr2, new_arr2):
            good = True
    elif stat == "avg":
        if qual_score(phred_score1) >= cutoff and qual_score(phred_score2) >= cutoff:
            good = True
    return good


import os
# __name__ is name associated with all the code. behind the scenes python thing
if __name__ == "__main__":
    # write tests for functions above
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    assert convert_phred("E") == 36, "Does not work on phred letter"
    assert qual_score("EEEEEEEE") == 36, "Does not work on phred string"
    assert np.array_equal(qs_arr("EEEE"), np.array([36, 36, 36, 36]))==True, "Does not work on phred string"
    assert gc_content("CTCACTTGGA") == .5, "Does not work on nucleotide seq"
    assert rev_comp("ACTGGTCACA") == "TGTGACCAGT", "Reverse complement does not work on DNA"
    assert rev_comp("ACUGGUCACA", True) == "UGUGACCAGU", "Reverse complement does not work on RNA"

    print("Passed DNA and RNA tests")
    print("Passed phred tests")
# tests for 1 line fa files
    with open("tester.txt", "w") as fh:
        fh.write(">lakjdfkjhskjfdhlk\n")
        fh.write("alsdhflkuhsdjfhlsdjf\n")
        fh.write("alsdhflkuhsdjfhlsdjf\n")
    oneline_fasta("tester.txt", "one_lined.txt")
    lines = 0
    with open("one_lined.txt", "r") as fh:
        for line in fh:
            lines += 1
    assert lines == 2, "Does not work on removing \n from seqs"

    if os.path.exists("one_lined.txt"):
        os.remove("one_lined.txt")
        print("removed one_lined.txt")
    else:
        print("The file does not exist") 
    if os.path.exists("tester.txt"):
        os.remove("tester.txt")
        print("removed tester.txt")
    else:
        print("The file does not exist")