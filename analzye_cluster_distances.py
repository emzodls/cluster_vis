import os
from glob import glob

os.chdir('/Volumes/lab_data/Ripps/ripp_codingGroups/id_ripp_clustering')

ave_spec_files = glob('*.avedist.csv')
ave_coding_files = glob('*.group_summary.tsv')

ave_coding_dist = {}
ave_spec_dist = {}

for ave_spec_file in ave_spec_files:
    cluster = ave_spec_file.split('.')[0]
    for line in open(ave_spec_file):
        line = line.strip()
        group,comps,ave_dist = line.split(',')
        ave_spec_dist[(cluster,group)] = (ave_dist,comps)

for ave_coding_file in ave_coding_files:
    cluster = ave_coding_file.split('.')[0]
    for line in open(ave_coding_file):
        line = line.strip()
        group,size,ave_dist = line.split(',')
        ave_coding_dist[(cluster,group)] = (ave_dist,size)

with open('distance_summaries.csv','w') as outfile:
    outfile.write('Cluster,Group,Size,Ave Cluster Dist, Ave Species Dist, Comps\n')
    for (cluster,group),(ave_cluster_dist,size) in ave_coding_dist.items():
        if int(size) >= 5 and (cluster,group) in ave_spec_dist:
            spec_dist,comps = ave_spec_dist[(cluster,group)]
            outfile.write('{},{},{},{},{},{}\n'.format(cluster,group,size,ave_cluster_dist,spec_dist,comps))

