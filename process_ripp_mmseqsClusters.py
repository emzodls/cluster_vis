from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import gzip
from pickle import dump,load
import subprocess,os,itertools,gzip
import igraph
from glob import glob
from itertools import combinations, zip_longest
import multiprocessing as mp
from math import ceil
import timeit
import concurrent.futures
import numpy as np
from itertools import groupby
from pickle import load

import os
from collections import defaultdict,Counter



os.chdir('/Users/emzodls/Documents')


clusters = defaultdict(set)
links = defaultdict(set)


for line in open('codingGroups_clus_cluster.tsv'):
    line = line.strip()
    tagA,tagB = line.split('\t')
    
    clusterA,geneA = tagA.split('|')
    clusterB, geneB = tagB.split('|')
    speciesA = geneA.split('_CDS_')[0]
    speciesB = geneB.split('_CDS_')[0]
    
    clusters[clusterA].add(geneA)
    clusters[clusterB].add(geneB)
    if geneA != geneB and speciesA != speciesB:
        links[geneA].add(tagB)
        links[geneB].add(tagA)

clusters = {k:v for k,v in clusters.items() if len(v) > 1}

connections = defaultdict(Counter)

clusterOrder = sorted(clusters.keys())

for cluster in clusterOrder:
    for gene in clusters[cluster]:
        for link in links[gene]:
            link_cluster, link_gene = link.split('|')

            connections[cluster][link_cluster] += 1
# for clusterID,genes in clusters.items():
#     for gene in genes:
#         for link in links[gene]:
#             link_cluster, link_gene = link.split('|')
#             connections[clusterID][link_cluster] += 1


dump(connections,open('connections.pkl','wb'))
g = igraph.Graph()
print('Found {} different assemblies, converting to nodes'.format(len(clusters)))
for cluster in clusterOrder:
    g.add_vertex(cluster)

edgesAndWeights = []

for idx,clusterA in enumerate(clusterOrder):
    print(idx+1,len(clusters),clusterA)
    for clusterB in clusterOrder:
        if clusterA != clusterB and connections[clusterA][clusterB] > 0 and connections[clusterB][clusterA] > 0:
            calculatedSimilarity = (connections[clusterA][clusterB]/len(clusters[clusterA]) +
                                    connections[clusterB][clusterA]/len(clusters[clusterB]))/2
            edgesAndWeights.append(((clusterA,clusterB),calculatedSimilarity))



# for edgenum,(species1,species2) in enumerate(simDict.keys()):
#     print('Added {} edge of {}'.format(edgenum,totalEdges))
#     g.add_edge(species1,species2, weight=simDict[(species1,species2)])
g.add_edges([x for x,y in edgesAndWeights])
g.es["weight"] = [y for x,y in edgesAndWeights]
print(edgesAndWeights)
print('Added {} edges'.format(len(edgesAndWeights)))
## get maximal cliques
#g.maximal_cliques(file=open('/Users/emzodls/Documents/cliques.txt','w'))

### Group with label propagation
subgraphs = g.community_multilevel()
membership = subgraphs.membership
with open('/Users/emzodls/Dropbox/Lab/Warwick/RiPP_nnets/node_spec_filt_multilevel.csv','w') as outfile:
    for vertex,membership in zip(g.vs,membership):
        outfile.write('{},{}\n'.format(vertex["name"],membership))
# ### Group with walkthrough
# dendrogram = g.community_walktrap()
# clusters_walktrap = dendrogram.as_clustering()
# membership_walktrap = clusters_walktrap.membership
# print('Finding neighborhoods with walktrap')
# with open('/Users/emzodls/Documents/node_membership_walktrap.csv','w') as outfile:
#     for vertex,membership in zip(g.vs,membership_walktrap):
#         outfile.write('{},{}\n'.format(vertex["name"],membership))
