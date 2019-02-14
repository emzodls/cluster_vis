import concurrent.futures
from collections import Counter,defaultdict
from pickle import dump,load
from glob import glob
import os
from Bio import SeqIO
# os.chdir('/Volumes/lab_data/Ripps/ripp_codingGroups')
# id_pfam_file = 'id_distance_pfam.csv'
# cand_pfam_file = 'distance_pfam.csv'
# asm2acc = load(open('acc2asm.pkl','rb'))
# codingSeqsDict = load(open('codingGroup2cds.pkl','rb'))
# putativePeptides = load(open('putativePeptides.pkl','rb'))
#
# candidate_clusters = set()
#
# for line in open(cand_pfam_file):
#     if not line.startswith('#'):
#         line_parse = line.split(',')
#         cluster,group = line_parse[0],line_parse[1]
#         families = load(open('ripp_clustering/{}.groups.pkl'.format(cluster),'rb'))
#         for coding_group in families[int(group)]:
#             asm = '_'.join(coding_group.split('_')[:-1])
#             if asm in asm2acc:
#                 candidate_clusters.add(coding_group)
#
# coding_seqs_to_update = defaultdict(dict)
# putativePeptideDict = dict()
# print(len(candidate_clusters))
# for candidate_cluster in candidate_clusters:
#     asm = '_'.join(candidate_cluster.split('_')[:-1])
#     acc = asm2acc[asm]
#     coding_seqs_to_update[acc][candidate_cluster] = codingSeqsDict[candidate_cluster]
#     putativePeptideDict[candidate_cluster] = codingSeqsDict[candidate_cluster] & putativePeptides

# dump(coding_seqs_to_update,open('coding_seqs_to_update.pkl','wb'))
# dump(putativePeptideDict,open('putativePeptideDict.pkl','wb'))

def get_coding_seqs(pair):
    coding_id_map = dict()
    file_header,coding_groups = pair
    coding_ids = set()
    for x in coding_groups.values():
        coding_ids.update(x)
    fasta_file_search = glob('/data/Groups/Theme2/rep_genomes_cds/{}*.fna'.format(file_header))
    #print(file_header,fasta_file_search)
    if len(fasta_file_search) == 1:
        fasta_file = fasta_file_search[0]
        # print(fasta_file)
        for entry in SeqIO.parse(fasta_file,'fasta'):
            entry_id_parse = entry.id.split('|')
            if entry_id_parse[3] in coding_ids:
                coding_id_map[entry_id_parse[3]] = entry.id
    return coding_id_map

if __name__ == "__main__":
    coding_seqs_to_update = load(open('coding_seqs_to_update.pkl','rb'))
    coding_map = dict()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(get_coding_seqs,path)
                   for path in coding_seqs_to_update.items()]
        for future in concurrent.futures.as_completed(futures):
            coding_id_map_part = future.result()
            coding_map.update(coding_id_map_part)
    dump(coding_map,open('coding_id_map.pkl','wb'))