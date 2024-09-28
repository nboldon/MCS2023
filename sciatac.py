#!/usr/bin/env python
import argparse
import subprocess
import sys
import os
import gzip
import io
import itertools
import time
import json
import collections
import re
from multiprocessing import Pool
import multiprocessing
import time
import glob

global bc1_dict, bc2_dict, bc3_list, cellranger_fastq1, cellranger_fastq2, cellranger_fastq3, parallelChunks, parallelSamples, CHUNK_SIZE
global bc2_correction_map, bc2_correction_map
global bc1_list, bc1ToIndex, bc2_list, bc2ToIndex, bc1_pat, whiteBarcodes, MEseqLen
global sampleDict, sampleList
global outputDir, input1, input2 

bc1_dict={71:['TGAGAGTTAT', 'GTCGCTAGCA', 'CGATCGTCGC', 'GCACTTGGCG', 'CACTTCGGTG', 'CTTACACTGA', 'GAGAGCTCTA', 'ATTCATCTTG', 'TTACACAACA', 'AACTTGATGC', 'ACGATCGCGA', 'CGTTGAACTT'], 72:['TGGAGCAATT', 'AGGCTTGTAG', 'GCATGATAAT', 'ACTATTAGCA', 'CAGGCACAAG', 'TATCTACGCG', 'TCGTCGTATG', 'CGCTACGCGA', 'GTAGAACGCG', 'TACGAGAAGA', 'GCTCTGTTAG', 'CTCACTGCAA']}

bc2_dict={51:['AAAGTCTGAC', 'GTTAATATTA', 'CAAGTCACCA', 'ACGAGCCATT', 'TCTAGAGTTA', 'TGCGACCTTC', 'TGTCATCCAT', 'TGATTACCTT'], 52:['AGCCGTGGAT', 'GCCTGACATA', 'GACGTGGTCA', 'CACGCTAATG', 'TGATGTGAAG', 'CCTAACTATA', 'TGTTTCAGAT', 'TTGACACTTT'], 53:['GAAGGACGAT', 'AGAGCGTACG', 'AGTTGGTGAC', 'GCACAGTGTA', 'TTTGCCGTCT', 'AGTCTATGTT', 'TACCATTATT', 'CTCTTGGTTC'], 54:['AGCGTAACCG', 'ACAAAGCAAC', 'GATAGAAGGC', 'ACACAAACTT', 'ACATACCGCT', 'CTTCGTGTCA', 'GCGCATTACG', 'TGACTTGTTA'], 55:['GATTGCTGCT', 'TTTCTTGAAC', 'ATTGAATTTC', 'TGAAGAAACT', 'TCTTATGAGA', 'ATAGCCTCTA', 'GCCTATATGC', 'CTAAGATCAA'], 56:['AGATTATTGA', 'CGAGTAAGGT', 'ACGTTACGAA', 'TACCACACAA', 'CGACTTTGAG', 'TGCTTGGACT', 'CCAGTCCTAT', 'TTGCAATAAA'], 57:['CATCACGGAT', 'GATGACCAAC', 'GTATAAGCTC', 'TGCGGAGAGA', 'GCGATAGGTT', 'TGTACCGCTC', 'CCTTAGGTCG', 'TCAAAGATGT']}


#bc2_dict={51:['AAAGTCTG', 'GTTAATAT', 'CAAGTCAC', 'ACGAGCCA', 'TCTAGAGT', 'TGCGACCT', 'TGTCATCC', 'TGATTACC'], 52:['AGCCGTGG', 'GCCTGACA', 'GACGTGGT', 'CACGCTAA', 'TGATGTGA', 'CCTAACTA', 'TGTTTCAG', 'TTGACACT'], 53:['GAAGGACG', 'AGAGCGTA', 'AGTTGGTG', 'GCACAGTG', 'TTTGCCGT', 'AGTCTATG', 'TACCATTA', 'CTCTTGGT'], 54:['AGCGTAAC', 'ACAAAGCA', 'GATAGAAG', 'ACACAAAC', 'ACATACCG', 'CTTCGTGT', 'GCGCATTA', 'TGACTTGT'], 55:['GATTGCTG', 'TTTCTTGA', 'ATTGAATT', 'TGAAGAAA', 'ACATTTAA', 'ATAGCCTC', 'GCCTATAT', 'CTAAGATC'], 56:['AGATTATT', 'CGAGTAAG', 'ACGTTACG', 'TACCACAC', 'CGACTTTG', 'TGCTTGGA', 'CCAGTCCT', 'TTGCAATA']}

bc3_list = ['ACAGCCGTCA', 'GCAGGCAAAC', 'GCTCCTTGCA', 'ACCAGCGTTG', 'CCAACAAGGC', 'CGGCACATCT', 'GTGGTCTCCA', 'GCCACACTGT', 'ACATATGCAG', 'GTCCAGAGTT', 'AGCCTTCGCA', 'AGTGTTGGAC', 'GTGTCCGCAT', 'GTGCGCACTG', 'GTCCGTGTTC', 'GCAGTAGATA', 'CCGAACCGAT', 'GCTTACCACG', 'AGCCTAGTGT', 'CGAGTTTCGG', 'TAAGCCAACT', 'GCGGACTTAG', 'TGGCTGGACG', 'CACGCGTAAC', 'CAACTTCGTG', 'CGGAGCCTGA', 'GTTGTGGCTC', 'TAGCAATCGG', 'GACTTTATCA', 'GTCGACTGGC', 'ATGTCAATAA', 'CTCGGAGCCA', 'ATATTCGACG', 'TAGGAGACTA', 'CCGTGCGGTA', 'TGGAAATGCT', 'GTTGCGAAAT', 'TAACGGTGGA', 'CCTTAGGCCT', 'GCTATAACAC', 'TCGGATCTTC', 'TGATCAGGAG', 'CCACTTAAGT', 'CCATCGTTTC', 'AGTTTGTCTA', 'TGCTCTCACG', 'GTCGTGATGG', 'TTACCGATGC', 'TACGGAATGT', 'TCGTGACCGC', 'GTGTATCGGT', 'CACATCGCAT', 'TGCGTGGCGA', 'CCGCTAAGAA', 'CTACTACACA', 'TTTCCTGTCG', 'TCTTGGCTTG', 'TTAAACCGGA', 'GTCTGAAGAG', 'TATATGAGCA', 'GTGCTAGTAG', 'TACGATCAAA', 'TCAATCCTCG', 'TGCGGTTGTG', 'TGTTGCTCGG', 'CCTGTGTGTG', 'CGTCTCACAG', 'TACGTTTCCT', 'CTTGAGAGGA', 'CCGTTGACGG', 'CTTTGATGGC', 'TAGTTCTCAC', 'TCGCTCGCGT', 'TGTATAGTTG', 'TTCTGCCGTC', 'TTGTGTATGG', 'TTCCTCAGCG', 'CTCTGCACGT', 'TCCTGCTTCT', 'CTTTATCCTG', 'TTGTACGTCA', 'TCCTTAGCCG', 'TTAGTTGGCA', 'CTGGTTACTT', 'TTCGCCACAA', 'CTTACTATTA', 'TTCCAAGGAA', 'TTCTCTGCGC', 'TTTGGTCATT', 'CTGATTTGAG', 'TCCGTACGAT', 'CCGTTTGTAT', 'TCGGTGTAAA', 'TTTAGCGCTA', 'TTCTAATCTT', 'TCCATCTGTA']

#bc3_list=['ACAGCCGT', 'GCAGGCAA', 'GCTCCTTG', 'ACCAGCGT', 'CCAACAAG', 'CGGCACAT', 'GTGGTCTC', 'GCCACACT', 'ACATATGC', 'GTCCAGAG', 'AGCCTTCG', 'AGTGTTGG', 'GTGTCCGC', 'GTGCGCAC', 'GTCCGTGT', 'GCAGTAGA', 'CCGAACCG', 'GCTTACCA', 'AGCCTAGT', 'CGAGTTTC', 'TAAGCCAA', 'GCGGACTT', 'TGGCTGGA', 'CACGCGTA', 'CAACTTCG', 'CGGAGCCT', 'GTTGTGGC', 'TAGCAATC', 'GACTTTAT', 'GTCGACTG', 'ATGTCAAT', 'CTCGGAGC', 'ATATTCGA', 'TAGGAGAC', 'CCGTGCGG', 'TGGAAATG', 'GTTGCGAA', 'TAACGGTG', 'CCTTAGGC', 'GCTATAAC', 'TCGGATCT', 'TGATCAGG', 'CCACTTAA', 'CCATCGTT', 'AGTTTGTC', 'TGCTCTCA', 'GTCGTGAT', 'TTACCGAT', 'TACGGAAT', 'TCGTGACC', 'GTGTATCG', 'CACATCGC', 'TGCGTGGC', 'CCGCTAAG', 'CTACTACA', 'TTTCCTGT', 'TCTTGGCT', 'TTAAACCG', 'GTCTGAAG', 'TATATGAG', 'GTGCTAGT', 'TACGATCA', 'TCAATCCT', 'TGCGGTTG', 'TGTTGCTC', 'CCTGTGTG', 'CGTCTCAC', 'TACGTTTC', 'CTTGAGAG', 'CCGTTGAC', 'CTTTGATG', 'TAGTTCTC', 'TCGCTCGC', 'TGTATAGT', 'TTCTGCCG', 'TTGTGTAT', 'TTCCTCAG', 'CTCTGCAC', 'TCCTGCTT', 'CTTTATCC', 'TTGTACGT', 'TCCTTAGC', 'TTAGTTGG', 'CTGGTTAC', 'TTCGCCAC', 'CTTACTAT', 'TTCCAAGG', 'TTCTCTGC', 'TTTGGTCA', 'CTGATTTG', 'TCCGTACG', 'CCGTTTGT', 'TCGGTGTA', 'TTTAGCGC', 'TTCTAATC', 'TCCATCTG']

MEseqLen=19
parallelSamples =20
parallelChunks = 8
CHUNK_SIZE = 100000
#cell ranger fastq file
cellranger_fastq1="_S1_L001_R1_001.fastq.gz"
cellranger_fastq2="_S1_L001_R2_001.fastq.gz"
cellranger_fastq3="_S1_L001_R3_001.fastq.gz"


def init_processes():
    global sample2fileHandle,sampleList 
    sample2fileHandle = {}
    pid=multiprocessing.current_process().pid
    for sample in sampleList:
        file1 = f"{outputDir}/{sample}/{pid}_R1.fq"
        file2 = f"{outputDir}/{sample}/{pid}_R2.fq"
        file3 = f"{outputDir}/{sample}/{pid}_R3.fq"
        sample2fileHandle[sample] = []
        sample2fileHandle[sample].append(open(file1, "wt"))
        sample2fileHandle[sample].append(open(file2, "wt"))
        sample2fileHandle[sample].append(open(file3, "wt"))
        
        
        
def get_chunks():
    global CHUNK_SIZE, input1, input2
    # If you have N processors,
    # then we need memory to hold 2 * (N - 1) chunks (one processor
    # is reserved for the main process).
    # The size of a chunk is CHUNK_SIZE * average-line-length.
    # If the average line length were 100, then a chunk would require
    # approximately 1_000_000 bytes of memory.
    # So if you had, for example, a 16MB machine with 8 processors,
    # you would have more
    # than enough memory for this CHUNK_SIZE.
    
    if1 = gzip.open(input1, "rb")
    if2 = gzip.open(input2, "rb")
    chunk = []
    linecount =0 
    while True:
        r1_name = if1.readline()
        if not r1_name: 
            break
        chunk.append(r1_name)
        chunk.append(if1.readline())
        t=if1.readline()
        chunk.append(if1.readline())
        
        t=if2.readline()
        chunk.append(if2.readline())
        t=if2.readline()
        chunk.append(if2.readline())
        linecount+=1
        if linecount == CHUNK_SIZE:
            yield chunk
            linecount=0
            chunk=[]
    if chunk:
        yield chunk

def process_chunk(chunk):
    global bc2_correction_map, bc3_correction_map, sampleDict, sampleList
    global sample2fileHandle
    global MEseqLen
    pid=multiprocessing.current_process().pid
    
    chunksize = len(chunk)
    
    total=0
    bc1err =0
    bc2err =0
    bc3err =0
    kept=0
        
    for i in range(0, chunksize, 5):
        total+=1
        r1_name = chunk[i].decode('utf-8').strip()
        r1_seq= chunk[i+1].decode('utf-8').strip()
        r1_qual= chunk[i+2].decode('utf-8').strip()
        r2_seq= chunk[i+3].decode('utf-8').strip()
        r2_qual= chunk[i+4].decode('utf-8').strip()      

        
        # Get barcodes and correct
        bc1, bc2, bc3, umi, R2_seqstr = get_barcode_seqs2(r1_name, r2_seq)     
        bc2 = correct_barcode(bc2, bc2_correction_map)
        bc3 = correct_barcode(bc3, bc3_correction_map)
        if not bc1:
            bc1err +=1
        if not bc2:
            bc2err +=1
        if not bc3:
            bc3err +=1
        if bc1 and bc2 and bc3:
            kept+=1
            bcindexStr = f"{bc1ToIndex[bc1]}-{bc2ToIndex[bc2]}"
            
            if bcindexStr in sampleDict:
                sample = sampleDict[bcindexStr]
                r2_name=re.sub(" 1", " 2", r1_name)
                r3_name=re.sub(" 1", " 3", r1_name)
                
                r2_seq = R2_seqstr[MEseqLen:]
                r2_qual = r2_qual[-(len(r2_seq)):]
                
                bcstr = bc1+bc2+bc3
                (sample2fileHandle[sample][0]).write(f"{r1_name}\n{r1_seq}\n+\n{r1_qual}\n")
                (sample2fileHandle[sample][1]).write(f"{r2_name}\t{bcstr}\n")
                (sample2fileHandle[sample][2]).write(f"{r3_name}\n{r2_seq}\n+\n{r2_qual}\n")
                #WH.write(f"{r1_name}\t{bcindexStr}\t{bc1+bc2+bc3}\t{r1_seq}\t{r1_qual}\t{R2_seqstr}\t{r2_qual}\n")
    for sample in sampleList:
        sample2fileHandle[sample][0].flush()
        sample2fileHandle[sample][1].flush()
        sample2fileHandle[sample][2].flush()
    return [total, bc1err, bc2err, bc3err, kept]

def closefh(x):
    global sample2fileHandle, sampleList
    pid=multiprocessing.current_process().pid
    print (f"Close filehandles for process {pid}")
    for sample in sampleList:
        sample2fileHandle[sample][0].close()
        sample2fileHandle[sample][1].close()
        sample2fileHandle[sample][2].close()
    time.sleep(1)
        
def swapbc(sample):
    global outputDir, cellranger_fastq1, cellranger_fastq2, cellranger_fastq3, whiteBarcodes
    print (f"Processing {sample}")
    processDir = f"{outputDir}/{sample}"
    R1files = glob.glob(f"{processDir}/*R1.fq")
    R2files = [re.sub(r'R1.fq$', 'R2.fq', i) for i in R1files]
    R3files = [re.sub(r'R1.fq$', 'R3.fq', i) for i in R1files]
    
    outfileR1 = f"{processDir}/{sample}{cellranger_fastq1}"
    outfileR2 = f"{processDir}/{sample}{cellranger_fastq2}"
    outfileR3 = f"{processDir}/{sample}{cellranger_fastq3}"
    
    FQ2H = gzip.open(outfileR2, "wt")
    bc2crbc={}
    whiteBarcodeIndex=0
    readcount = 0
    for r2f in R2files:
        with open (r2f, "rt") as R2H:
            for line in R2H:
                readcount +=1
                name, cellbc = line.strip().split("\t")
                crbc = ""
                if cellbc in bc2crbc:
                    crbc = bc2crbc[cellbc]
                else:
                    crbc = whiteBarcodes[whiteBarcodeIndex]
                    whiteBarcodeIndex+=1
                    bc2crbc[cellbc] = crbc
                
                FQ2H.write(f"{name}\n{crbc}\n+\n{qstr}\n")
    
    FQ2H.close()
    os.system(f"rm {processDir}/*_R2.fq")
    os.system(f"cat {' '.join(R1files)} | gzip -c > {outfileR1}")
    os.system(f"rm {processDir}/*_R1.fq")
    os.system(f"cat {' '.join(R3files)} | gzip -c > {outfileR3}")
    os.system(f"rm {processDir}/*_R3.fq")
    return [sample, readcount, whiteBarcodeIndex]
            
def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of whitelists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for mismatch_whitelist in mismatch_map:
        corrected = mismatch_whitelist.get(barcode, None)

        if corrected:
            return corrected

    return None


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    Yield:
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode whitelist.
    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = [None] * (edit_distance + 1)

    mismatch_to_whitelist_map[0] = {k: k for k in whitelist}

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(list(whitelist))

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                del mismatch_to_whitelist_map[mismatch_count][mismatch]

    return mismatch_to_whitelist_map


def reverse_complement(x):
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    xrev = x[::-1]
    xrevcomp = ''.join([complements[z] for z in xrev])
    return xrevcomp



def get_barcode_seqs2(r1_name, r2_seq):
    global bc1_pat
    """
    Extract the correct sequences from the R1 name.
    """
    matches = re.search("(\\w+)\\+(\\w+)$", r1_name)
    if not matches:
        return "", "", "", "", ""
    bc3=matches[1]
    bc2=matches[2]
    
    matches = re.match(bc1_pat, r2_seq)
    umi=""
    bc1=""
    seqstr = ""
    if matches:
        umi=matches[1]
        bc1=matches[2]
        seqstr=matches[3]
    return bc1, bc2, bc3, umi, seqstr

            
def pre_process_barcode(editDist):
    global bc2_correction_map, bc3_correction_map, bc1_dict, bc2_dict, bc3_list
    global bc1_list, bc1ToIndex, bc2_list, bc2ToIndex
    global bc1_pat  
    
    bc1_list=[]
    bc1ToIndex = {}
    for myindex in bc1_dict:
        mytmplist = bc1_dict[myindex]
        for bc in mytmplist:
            bc1ToIndex[bc] = myindex
            bc1_list.append(bc)

    bc2_list=[]
    bc2ToIndex = {}
    for myindex in bc2_dict:
        mytmplist = bc2_dict[myindex]
        for bc in mytmplist:
            bc2ToIndex[bc] = myindex
            bc2_list.append(bc)

    bc1_pat = "|".join(bc1_list)
    bc1_pat = "(\\w{4,9})("+bc1_pat + ")(\\w+)"
    bc2_correction_map = construct_mismatch_to_whitelist_map(set(bc2_list), editDist)
    bc3_correction_map = construct_mismatch_to_whitelist_map(set(bc3_list), editDist)


if __name__ == '__main__':
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description='A program to convert custom single cell atac data to Cellranger compatible data files')
    parser.add_argument('-1', '--input1', required=True, help='Input R1.fastq.gz file. Must be a .gz file')
    parser.add_argument('-2', '--input2', required=True, help='Input R2.fastq.gz file. Must be a .gz file')
    parser.add_argument('-s', '--samplesheet', required=True, help='Samplesheet describing the layout of the samples. It is a tab-delimited text file with three columns: sampleName, bc1, bc2')
    parser.add_argument('-o', '--outputdir', required=True, help='Output directory name')
    parser.add_argument('-f', '--overwrite', action=argparse.BooleanOptionalAction, help='(Optional) Overwrite existing output directory')
    parser.add_argument('-e', '--editdistance', required=False, default=1, help='(Optional) Edit distance for barcode matching, default is 1')
    parser.add_argument('-w', '--whitelist', required=False, default=scriptDir+"/cratac_curated.txt.gz", help='(Optional) By default, the script will use a file named cratac_curated.txt.gz located in the same directory as the script. This file was curated from the Cellranger whiteList file, based on the original whitelist file cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz. See README appendix how the file is curated')

    args = parser.parse_args()

    input1 = args.input1
    input2 = args.input2
    editDist = args.editdistance
    whiteList = args.whitelist
    sampleFile = args.samplesheet
    outputDir = args.outputdir
    logfile=f"{outputDir}/run_stats"
    
    if isinstance(editDist, str):
        if editDist.isnumeric():
            editDist=int(editDist)
        else:
            print("-e --editdistance must be an integer value")
            sys.exit()

    pre_process_barcode(editDist)
    

    
    #validate input argument

    if not (input1.endswith(".gz") and input2.endswith(".gz")):
        print("Input fastq files must be .gz files")
        sys.exit()
    if not (whiteList.endswith(".gz")):
        print("The white list file must be a .gz file")
        sys.exit()
    if not os.path.isfile(input1):
        print(f"The input file {input1} does not exist")
        sys.exit()
    if not os.path.isfile(input2):
        print(f"The input file {input1} does not exist")
        sys.exit()
    if not os.path.isfile(whiteList):
        print(f"The whitelist barcode file {whiteList} does not exist. Specify the file path with -w option")
        sys.exit()

    if os.path.exists(outputDir):
        if args.overwrite:
            print(f"The output directory {outputDir} exists. As -f is specified, it will be overwritten.")
            os.system(f"rm -fr {outputDir}")
        else:
            print(f"The output directory {outputDir} exists. Specify the -f to overwrite.")
            sys.exit()

    os.makedirs(outputDir)
    
    #process samplesheet
    sampleDict = {}
    sampleList = []
    sampleCount=0
    with open (sampleFile, "r") as IN:
        for line in IN:
            line= line.strip()
            if len(line)>0:
                sampleCount+=1
                sample, indexBC1, indexBC2 = line.split("\t")
                sample = re.sub("\s", "", sample)
                indexBC1 = re.sub("\s", "", indexBC1)
                indexBC2 = re.sub("\s", "", indexBC2)
                sampleDict[f"{indexBC1}-{indexBC2}"] = sample
                if sample not in sampleList:
                    sampleList.append(sample)
                    os.mkdir(f"{outputDir}/{sample}")
    IN.close()
    if sampleCount==0:
        print(f"There is no sample in samplesheet {sampleFile}. Exit now!")
        sys.exit()
    else:
        print (f"There are {sampleCount} samples in sample sheet.\nStart stage 1: demultiplexing samples")
    sum_total =0
    sum_bc1err=0
    sum_bc2err=0
    sum_bc3err=0
    sum_kept=0
    with Pool(parallelChunks, initializer=init_processes,) as pool:
        for t in pool.imap_unordered(process_chunk,  get_chunks()):
            sum_total +=t[0]
            sum_bc1err+=t[1]
            sum_bc2err+=t[2]
            sum_bc3err+=t[3]
            sum_kept+=t[4]
            print(f"Processed reads: {sum_total}")
            print(t)
            
        pool.map(closefh, range(0, parallelChunks))
        pool.close()
        pool.join()
    
    print("Stage 1 de-multiplexing finished")
    
    
    print("Start stage 2: replace cell barcode with cellranger whitelist")
    whiteBarcodes = []
    
    with gzip.open(whiteList, "rb") as IN:
        for line in IN:
            whiteBarcodes.append(line.decode('utf-8').strip())
    IN.close()
    qstr="I"*(len(whiteBarcodes[0]))
    sample2count = {}
    with Pool(parallelSamples) as pool:
        for t in pool.imap_unordered(swapbc, sampleList):
            print(t)
            sample2count[t[0]] = f"{t[1]}\t{t[2]}"
    pool.close()
    pool.join()
    
    print("Stage 2 change_to_cellranger_barcode finished")
    LOGW=open(logfile, "wt")
    LOGW.write(f"TotalReads:\t{sum_total}\n")
    LOGW.write(f"BC1err:\t{sum_bc1err}\n")
    LOGW.write(f"BC2err:\t{sum_bc2err}\n")
    LOGW.write(f"BC3err:\t{sum_bc3err}\n")
    LOGW.write(f"UsedReads:\t{sum_kept}\n")
    LOGW.write(f"Sample\tReads\tCells\n")
    for sample in sampleList:
        LOGW.write(f"{sample}\t{sample2count[sample]}\n")
        
    LOGW.close()



