from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import os
import gzip
import regex
from pickle import dump,load
import subprocess,os,itertools,gzip
from glob import glob
from itertools import combinations, zip_longest
import multiprocessing as mp
from math import ceil
import timeit
import concurrent.futures
import numpy as np
from itertools import groupby
from pickle import load


minProteinLength = 10
maxProteinLength = 150

putativePeptideDictPath = '/data/Groups/Theme2/putative_CDS.pkl'
putativePeptideDict = load(open(putativePeptideDictPath,'rb'))
codingGroupsNoLeader = "/home/u1474301/ripp_allVall/codingGroups_noleader.fa"
codingGroupsFile = "/home/u1474301/ripp_allVall/codingGroups.fa"

putativePeptides = set()
for peptides in putativePeptideDict.values():
    putativePeptides.update(peptides)

codingGroups = glob('/data/Groups/Theme2/ripp_codingGroups_corrected/*.fa')
print(len(codingGroups))
for idx,clusterFile in enumerate(codingGroups):
    if idx % 1000 == 0:
        print(idx/len(codingGroups) * 100)
    no_leader = []
    cluster = os.path.splitext(os.path.split(clusterFile)[1])[0]
    seqs = [x for x in SeqIO.parse(clusterFile,'fasta')]
    for seq in seqs:
        new_seq_id = '{}|{}'.format(cluster,seq.id)
        new_seq_name = '{}|{}'.format(cluster, seq.name)
        if seq.id not in putativePeptides:
            no_leader.append(seq)
        seq.id = new_seq_id
        seq.name  = new_seq_name

    with open(codingGroupsFile,'a') as outfile:
        SeqIO.write(seqs,outfile,'fasta')
    with open(codingGroupsNoLeader,'a') as outfile:
        SeqIO.write(no_leader,outfile,'fasta')


