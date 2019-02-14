from math import log
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import AffinityPropagation
import os
from itertools import combinations
import concurrent.futures
from collections import Counter,defaultdict
from pickle import dump
from glob import glob

class hitDict():
    def __init__(self,species,id):
        self.species = species
        self.id = id
        self.score_dict = dict()
        self.maxscore = 1.
        self.hits = set()

    def add_hit(self,hit_id,score):
        if score > self.maxscore:
            self.maxscore = score
        if hit_id in self.hits:
            if self.score_dict[hit_id] < score:
                self.score_dict[hit_id] = score
        else:
            self.score_dict[hit_id] = score
            self.hits.add(hit_id)
    def get(self,key,default=None):
        return self.score_dict.get(key,default)
    def __repr__(self):
        return repr(self.score_dict)

def parseBLAST(path,swapQuery=False,evalCutoff=10,scoreCutoff=0):
    proteins = dict()
    if swapQuery:
        queryIdx = 1
        hitIdx = 0
    else:
        queryIdx = 0
        hitIdx = 1
    with open(path) as blast_handle:
        for line in blast_handle:
            try:
                line_parse = line.split('\t')
                species_id,query = line_parse[queryIdx].split('|')
                protein = proteins.setdefault(query,hitDict(species_id,query))
                hit_id = line_parse[hitIdx].split('|')[1]
                score  = float(line_parse[11])
                eval = float(line_parse[10])
                if eval <= evalCutoff and score >= scoreCutoff:
                    protein.add_hit(hit_id,score)
            except (ValueError,IndexError):
                pass
    return proteins

def calculate_distance(protein, target,linearDist = True):
    if protein.id == target.id:
        return 0.
    else:
        ## handle non-reciprocal hits
        if target.id in protein.hits and \
                        protein.id not in target.hits:
            if not linearDist:
                normScore = -log(protein.score_dict.get(target.id) / protein.maxscore)
                maxDist = -log(1 / protein.maxscore)
                if maxDist == 0:
                    return 1
                else:
                    return normScore / maxDist
            else:
                distance = 1 - (protein.score_dict.get(target.id) / protein.maxscore)
                return distance

        elif target.id not in protein.hits and \
                        protein.id in target.hits:
            if not linearDist:
                normScore = -log(target.score_dict.get(protein.id) / target.maxscore)
                maxDist = -log(1 / target.maxscore)
                if maxDist == 0:
                    return 1
                else:
                    return normScore / maxDist
            else:
                distance = 1 - (target.score_dict.get(protein.id) / target.maxscore)
                return distance

        else:
            if len(protein.hits) == 0  or len(target.hits) == 0:
                proteinAbundanceFactor = 1
                targetAbundanceFactor = 1
            else:
                proteinAbundanceFactor = len(protein.hits & target.hits) / float(len(protein.hits))
                targetAbundanceFactor =  len(protein.hits & target.hits)  / float(len(target.hits))

            # if both the protein abundance factor and the target abundance factor are 0 this means that there
            # aren't any common blast hits between the proteins so the distance is 1
            if proteinAbundanceFactor == 0 and targetAbundanceFactor == 0:
                return 1.
            elif not linearDist:
                normScore = protein.score_dict.get(target.id, 1) / protein.maxscore
                normReciprocalScore = target.score_dict.get(protein.id, 1) /target.maxscore
                maxDist = -log((proteinAbundanceFactor / protein.maxscore +
                                targetAbundanceFactor / target.maxscore) / 2)
                distance = -log(((normScore * proteinAbundanceFactor) + (normReciprocalScore * targetAbundanceFactor)) / 2)
                if maxDist == 0:
                    return 1
                else:
                    return distance / maxDist
            else:
                normScore = protein.score_dict.get(target.id, 0) / protein.maxscore
                normReciprocalScore = target.score_dict.get(protein.id, 0) / target.maxscore
                distance = (proteinAbundanceFactor * (1 - normScore) \
                           + targetAbundanceFactor * (1 - normReciprocalScore)) \
                           / (proteinAbundanceFactor+targetAbundanceFactor)
                return distance


def calculateDistBitScore(cluster1,cluster2,linearDist = True):
    '''
    another way to measure distance that will pair up the proteins using matching
    then using the bitscores of the pairs of proteins added up will calculate the distance
    in a similar way to how protein distance is calculated:

    linear: 1 - (sum bit score of paired matches/max bit score of cluster)

    scaling will work this way:

    If there are reciprocal scores average the distance. Otherwise, if hits are on smaller
    cluster, (small clus/large clus)*dist + (1 - small clus/large clus)

    '''
    clus1Size = len(cluster1)
    clus2Size = len(cluster2)

    scoreMatrix =  np.ndarray((clus1Size,clus2Size))

    # populate the score matrix if there are any proteins that are "close together"
    for i,proteinI in enumerate(cluster1):
        for j,proteinJ in enumerate(cluster2):
            scoreMatrix[i,j] = calculate_distance(proteinI,proteinJ,linearDist=linearDist)

    # get the pairings
    pairings = [(x,y) for x,y in zip(*linear_sum_assignment(scoreMatrix)) if (x<clus1Size) and (y<clus2Size)]
    pairScores = [(scoreMatrix[(x,y)],cluster1[x],cluster2[y]) for x,y in pairings]
    pairs = [(x,y,z) for x,y,z in pairScores if x < 1.]
    pairs.sort(key=lambda x:x[0])
    # figure out which cluster has the hits and calculate coverage

    clus1maxBitScore = sum(protein.maxscore for protein in cluster1) + 0.00000001
    clus2maxBitScore = sum(protein.maxscore for protein in cluster2) + 0.00000001
    # precision errors
    clus1cumBitScore = sum(proteinI.score_dict.get(proteinJ.id,0) for dist,proteinI,proteinJ in pairs)
    clus2cumBitScore = sum(proteinJ.score_dict.get(proteinI.id,0) for dist,proteinI,proteinJ in pairs)
    assert clus1cumBitScore <= clus1maxBitScore and clus2cumBitScore <= clus2maxBitScore
    return 0.5*(1 - clus1cumBitScore/clus1maxBitScore) + 0.5*(1-clus2cumBitScore/clus2maxBitScore),[(x,y.id,z.id) for x,y,z in pairs]

def group_clusters(proteins,cutoff=0.7):
    species_dict = {}
    for protein in proteins.values():
        species_prots = species_dict.setdefault(protein.species, [])
        species_prots.append(protein)
    distance_dict = dict()
    cluster_order = sorted(species_dict.keys())
    simMatrix = np.zeros((len(cluster_order), len(cluster_order)), dtype=np.float32)
    for idxI, clusterI in enumerate(cluster_order):
        for idxJ, clusterJ in enumerate(cluster_order):
            order = tuple(sorted([clusterI, clusterJ]))
            if clusterI == clusterJ:
                simMatrix[idxI, idxJ] = 1
            else:
                if clusterI != clusterJ and order not in distance_dict:
                    distance_dict[order] = calculateDistBitScore(species_dict[clusterI], species_dict[clusterJ])[0]
                if distance_dict[order] <= cutoff:
                    simMatrix[idxI, idxJ] = 1 - distance_dict[order]
                else:
                    simMatrix[idxI, idxJ] = 0
    af = AffinityPropagation(damping=0.9, max_iter=1000, convergence_iter=200, affinity="precomputed").fit(simMatrix)
    labelsSub = af.labels_

    return distance_dict,sorted(list(zip(cluster_order, labelsSub)), key=lambda x: x[1])

def summarize_clusters(distance_dict,classification):
    groups = defaultdict(set)
    ave_dist = {}

    for cluster_id, group in classification:
        groups[group].add(cluster_id)
    counts = {k:len(v) for k,v in groups.items()}
    for group in groups.keys():
        distance = 0
        ctr = 0
        for pair in combinations(groups[group],2):
            order = tuple(sorted(pair))
            distance += distance_dict[order]
            ctr += 1
        if ctr == 0:
            ave_dist[group] = 0
        else:
            ave_dist[group] = distance/ctr
    return groups,counts,ave_dist


def batch_runner(path):
    try:
        base_name = os.path.split(path)[1].split('.')[0]
        proteins = parseBLAST(path)
        distance_dict, clusters = group_clusters(proteins)
        result = summarize_clusters(distance_dict, clusters)
        return base_name,result
    except AssertionError:
        return base_name, None

if __name__ == "__main__":
    taskList = glob('/data/Groups/Theme2/id_ripp_families/*tsv')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(batch_runner,path)
                   for path in taskList]
        for future in concurrent.futures.as_completed(futures):
            outname, result = future.result()
            print(outname)
            if result:
                groups, counts, ave_dist = result
                with open(os.path.join('/data/Groups/Theme2/id_ripp_clustering/','{}.groups.pkl'.format(outname)),'wb') as outfile:
                    dump(groups,outfile)
                with open(os.path.join('/data/Groups/Theme2/id_ripp_clustering/', '{}.group_summary.tsv'.format(outname)),
                          'w') as outfile:
                    for group in groups.keys():
                        outfile.write('{},{},{}\n'.format(group,counts[group],ave_dist[group]))
                with open('/data/Groups/Theme2/id_ripp_clustering/summary.csv','a') as outfile:
                    outfile.write('{},{}\n'.format(outname, len(groups.keys())))
            else:
                with open('/data/Groups/Theme2/id_ripp_clustering/errors.csv', 'a') as outfile:
                    outfile.write('{}\n'.format(outname))