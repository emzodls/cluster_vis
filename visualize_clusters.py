import json
import numpy as np
from pickle import load
from calculate_cluster_distances import calculateDistBitScore
codingGroup2cds = load(open('../codingGroup2cds.pkl','rb'))


def generate_cluster_dict(cluster,proteinMap = None):
    ## assumes protein is in standard clusterTools format, otherwise a map needs to be provided
    cluster_dict = dict()
    if proteinMap:
        mapped_cluster = [proteinMap[i] for i in cluster]
        back_map = dict(zip(mapped_cluster,cluster))
        ordered_cluster = sorted(mapped_cluster,
                                 key = lambda x: int(x.split('|')[1].split('-')[0]))
    else:
        back_map = None
        ordered_cluster = sorted(cluster,
                                 key = lambda x: int(x.split('|')[1].split('-')[0]))
    first_prot = ordered_cluster[0]
    last_prot = ordered_cluster[-1]
    cluster_dict['start'] = int(first_prot.split('|')[1].split('-')[0])
    cluster_dict['end'] = int(last_prot.split('|')[1].split('-')[0])
    cluster_dict['offset'] = cluster_dict['start'] - 1
    cluster_dict['size'] = cluster_dict['end'] - cluster_dict['offset']

    orfs = []
    for protein in ordered_cluster:
        protein_dict = dict()
        protein_parse = protein.split('|')
        if back_map:
            protein_dict['id'] = back_map[protein]
        else:
            protein_dict['id'] = protein_parse[-1]
        location = protein_parse[1].split('-')
        protein_dict['start'] =  int(location[0]) - cluster_dict['offset']
        protein_dict['end'] = int(location[1]) - cluster_dict['offset']
        if protein_parse[2] == '+':
            protein_dict['strand'] = 1
        else:
            protein_dict['strand'] = -1

        protein_dict['color'] = 'rgb(211,211,211)'
        orfs.append(protein_dict)
        cluster_dict['orfs'] = orfs
    return cluster_dict


def generate_json(clusterDict,hitDict,proteinMap = None):
    ## Find reference matrix
    json = []
    cluster_order = sorted(list(clusterDict.keys()))
    dist_matrix = np.zeros((len(cluster_order), len(cluster_order)), dtype=np.float32)
    pair_dict = dict()
    for idxI, clusterI in enumerate(cluster_order):
        for idxJ, clusterJ in enumerate(cluster_order):
            order = tuple(sorted([clusterI, clusterJ]))
            if clusterI != clusterJ and order not in pair_dict:
                dist_matrix[idxI,idxJ], pair_dict[order] = calculateDistBitScore([hitDict[i] for i in clusterDict[clusterI]],
                                                                      [hitDict[i] for i in clusterDict[clusterJ]])
    dist_matrix = dist_matrix + dist_matrix.T - np.diag(dist_matrix.diagonal())
    sorted_indices = np.argsort(dist_matrix.mean(1))
    ref_cluster = cluster_order[sorted_indices[0]]

    if proteinMap:
        ref_cluster_dict = generate_cluster_dict(clusterDict[ref_cluster],proteinMap = proteinMap)
    else:
        ref_cluster_dict = generate_cluster_dict(clusterDict[ref_cluster])
    ref_cluster_dict['id'] = ref_cluster

    json.append(ref_cluster_dict)

    for idx in sorted_indices[1:]:
        cluster_id = cluster_order[idx]
        if proteinMap:
            cluster_dict = generate_cluster_dict(clusterDict[cluster_id],proteinMap = proteinMap)
        else:
            cluster_dict = generate_cluster_dict(clusterDict[cluster_id])
        cluster_dict['id'] = cluster_id

        json.append(cluster_dict)

    scale = max(x['size'] for x in json)

    return scale, json

