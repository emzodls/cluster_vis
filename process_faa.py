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


def runFastaParser(path):
    filter_orfs = []
    baseName = os.path.split(os.path.splitext(path)[0])[1] + '_filtered.fa'
    for entry in SeqIO.parse(path, "fasta"):
        if len(entry.seq) >= minProteinLength and len(entry.seq) <= maxProteinLength:
            filter_orfs.append(entry)
    return baseName,filter_orfs

def getContigGroups(path):
    groups = []
    for x in SeqIO.parse(path, 'fasta'):
        idParse = x.id.split('|')
        sequence = str(x.seq)
        if not '*' in sequence[:-1]:
            idParse.append(sequence)
            groups.append(idParse)
    species = groups[0][0]
    codingGroups = [dict((x[3],x[-1]) for x in list(g)) for k,g in groupby(groups,lambda x: x[2])]
    codingGroupPeptides = []
    peptidesAdded = set()
    if species in putativePeptideDict:
        for peptide in putativePeptideDict[species]:
            for group in codingGroups:
                if peptide in group and peptide not in peptidesAdded:
                    codingGroupPeptides.append(group)
                    peptidesAdded.update(group.keys())
                    break

    return species,codingGroupPeptides

if __name__ == "__main__":
    taskList =  glob('/data/Groups/Theme2/rep_genomes_cds/*fna')
    #taskList = glob('/Volumes/lab_data/Ripps/sample_genomes_faa/samples/*fna')
    # for task in taskList:
    #     outname, filter_orfs =  runFastaParser(task)
    #     with open(os.path.join('/Volumes/lab_data/sample_ripp', outname), 'w') as handle:
    #         SeqIO.write([orf[3] for orf in filter_orfs], handle, 'fasta')
    #alreadyDone  =  glob('/data/Groups/Theme2/rep_genomes_cds_filtered/*_filtered.fa')
    #orfs = set(os.path.split(x)[1] for x in alreadyDone)
    with concurrent.futures.ProcessPoolExecutor(max_workers=35) as executor:
        futures = [executor.submit(getContigGroups,path)
                   for path in taskList]
        for future in concurrent.futures.as_completed(futures):
            species, codingGroupPeptides = future.result()
            print(species)
            for idx,codingGroup in enumerate(codingGroupPeptides):
                #with open(os.path.join('/Volumes/lab_data/sample_ripp_fix', outname), 'w') as handle:
                if len(codingGroup) > 2:
                    with open(os.path.join('/data/Groups/Theme2/ripp_codingGroups_corrected/',species) + '_{}.fa'.format(idx), 'w') as handle:
                        for k,v in codingGroup.items():
                            handle.write('>{}\n{}\n'.format(k,v))


