import json,math,colorsys
import numpy as np
from collections import defaultdict
from pickle import load
from calculate_cluster_distances import calculateDistBitScore
codingGroup2cds = load(open('../codingGroup2cds.pkl','rb'))
from fractions import Fraction
from itertools import chain,count

'''
color generation from:
http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors
'''
def zenos_dichotomy():
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in count():
        yield Fraction(1,2**k)

def getfracs():
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield 0
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5) # can be used for the v in hsv to map linear values 0..1 to something that looks equidistant

def genhsv(h):
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range

genrgb = lambda x: colorsys.hsv_to_rgb(*x)

flatten = chain.from_iterable

def _get_colors_Janus(num_colors):
    fracGen = getfracs()
    fracs = [next(fracGen) for i in range(int((num_colors+1)/2))]
    rgbs = list(map(genrgb,flatten(list(map(genhsv,fracs)))))
    return rgbs[:num_colors]


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
        cluster_dict['dist'] = '{:4f}'.format(dist_matrix.mean(1)[idx])
        json.append(cluster_dict)

    scale = max(x['size'] for x in json)

    protein_dict = dict()
    for cluster_dict in json:
        protein_dict.update({x['id']:(x,cluster_dict['id']) for x in cluster_dict['orfs']})

    # Start with pairs from the reference cluster and make your way down
    homolog_groups = defaultdict(set)
    proteins_in_groups = set()

    for ctr, idxI in enumerate(sorted_indices):
        clusterI = cluster_order[idxI]
        for idxJ in sorted_indices[ctr+1:]:
            clusterJ = cluster_order[idxJ]
            order = tuple(sorted([clusterI, clusterJ]))
            pairs = pair_dict[order]

    for ctr, idxI in enumerate(sorted_indices):
        clusterI = cluster_order[idxI]
        for idxJ in sorted_indices[ctr+1:]:
            clusterJ = cluster_order[idxJ]
            order = tuple(sorted([clusterI, clusterJ]))
            pairs = pair_dict[order]
            for (dist,proteinA, proteinB) in pairs:
                if proteinA not in proteins_in_groups or proteinB not in proteins_in_groups:
                    proteinA_dict,clusterA = protein_dict[proteinA]
                    if clusterA == clusterI:
                        proteinI = proteinA
                        proteinJ = proteinB
                    elif clusterA == clusterJ:
                        proteinI = proteinB
                        proteinJ = proteinA

                    homolog_groups[proteinI].add((proteinI,0))
                    homolog_groups[proteinI].add((proteinJ,dist))
                    proteins_in_groups.add(proteinI)
                    proteins_in_groups.add(proteinJ)
    num_colors = len(homolog_groups)

    colors = _get_colors_Janus(num_colors)
    hsv_grey = (0.0, 0.0, Fraction(117, 255))
    for idx,ref in enumerate(homolog_groups.keys()):
        color = colors[idx]
        print(list(map(lambda x: int(x * 255), color)))
        for homolog,dist in homolog_groups[ref]:
            color_hsv = colorsys.rgb_to_hsv(*color)
            scaled_h = (hsv_grey[0]-color_hsv[0]) * dist + color_hsv[0]
            scaled_s = (hsv_grey[1]-color_hsv[1]) * dist + color_hsv[1]
            scaled_v = (hsv_grey[2]-color_hsv[2]) * dist + color_hsv[2]
            scaled_color = list(map(lambda x: int(x*255),colorsys.hsv_to_rgb(scaled_h,scaled_s,scaled_v)))
            protein_dict[homolog][0]['color'] = 'rgb({},{},{})'.format(*scaled_color)
            if homolog != ref:
                protein_dict[homolog][0]['pair'] = '{}'.format(ref)
                protein_dict[homolog][0]['dist'] = '{:.4f}'.format(dist)

    return scale, json

