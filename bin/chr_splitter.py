#!/usr/bin/env python
from Bio import SeqIO
import argparse

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--chr", type=str, required=True)
    parser.add_argument("-r", "--ref", type=str, required=True)
    parser.add_argument("-p", "--pie", type=int, required=True)
    args = parser.parse_args()
    return args
arg = GetArguments()

chr = arg.chr
ref = SeqIO.parse(arg.ref,"fasta")
piece = arg.pie

for rec in ref:
    id = rec.id
    if id == chr:
        nparts = len(rec.seq)/piece
        if nparts < 1:
            line = chr
            print(line)
        else:
            for i in range(1,int(nparts) + 1):
                if i == 1:
                    start = 0
                    end = piece * i
                elif i == int(nparts):
                    start = piece * i + 1
                    end = len(rec.seq)
                else:
                    start = piece * (i-1) + 1
                    end = piece * i
                line = chr+":"+str(start)+"-"+str(end)
                print(line+"\n")
