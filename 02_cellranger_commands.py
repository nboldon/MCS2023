#!/usr/bin/python

import argparse
import os
import sys

parser=argparse.ArgumentParser(description="A program to create parallel commands to run cellranger")
parser.add_argument('-i', '--input', required=True, help='Input directory')
parser.add_argument('-o', '--output', required=True, help='Output commands file')
parser.add_argument('-r', '--ref', required=True, help='reference')

args = parser.parse_args()

inputdir = args.input
outputcommands = args.output
reference = args.ref

#validate
if not os.path.exists(inputdir):
    print(f"Error: input directory {inputdir} does not exist! ")
    sys.exit()

    
if not os.path.exists(reference):
    print(f"Error: reference {reference} does not exist! ")
    sys.exit()

inputdir = os.path.abspath(inputdir)
reference = os.path.abspath(reference)

samples = os.listdir(inputdir)


WH = open (outputcommands, "wt")
for sample in samples:
    WH.write(f"cd {inputdir}; cellranger-atac count --id={sample}_cr --reference={reference} --fastqs={inputdir}/{sample} --sample={sample} --localcores=32 --localmem=160\n")
WH.close()
