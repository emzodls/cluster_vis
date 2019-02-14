import os
from collections import defaultdict
from pickle import load

def get_pfam_doms(cluster,group,pfam_path='id_ripp_families/',groups_path='id_ripp_clustering/'):
    coding_groups_dict = load(open(os.path.join(groups_path,'{}.groups.pkl'.format(cluster)),'rb'))
    coding_groups = coding_groups_dict[int(group)]
    domtbl_file_antismash = os.path.join(pfam_path,'{}.antismash.domtbl'.format(cluster))
    domtbl_file_prism = os.path.join(pfam_path, '{}.prism.domtbl'.format(cluster))
    pfam_dict = defaultdict(set)
    # print(coding_groups)
    for line in open(domtbl_file_antismash):
        line = line.strip()
        if not line.startswith('#'):
            line_parse = line.split()
            coding_group,cds = line_parse[0].split('|')
            pfam_hit =  line_parse[3]
            if coding_group in coding_groups:
                pfam_dict[coding_group].add(pfam_hit)
    for line in open(domtbl_file_prism):
        line = line.strip()
        if not line.startswith('#'):
            line_parse = line.split()
            coding_group,cds = line_parse[0].split('|')
            pfam_hit =  line_parse[3]
            if coding_group in coding_groups:
                pfam_dict[coding_group].add(pfam_hit)
    return pfam_dict


os.chdir('/Volumes/lab_data/Ripps/ripp_codingGroups/')

with open('id_distance_ripp_hmms.csv','w') as outfile:
    for line in open('id_ripp_clustering/distance_summaries.csv'):
        if line.startswith('#'):
            outfile.write(line)
        else:
            line = line.strip()
            line_parse = line.split(',')
            pfam_dict = get_pfam_doms(line_parse[0],line_parse[1])
            if pfam_dict:
                outfile.write('{},{}\n'.format(line,set.intersection(*pfam_dict.values())))
            else:
                print(line_parse[0],line_parse[1])

# pfam_dict = get_pfam_doms('170133','3')
# print(pfam_dict)
# print(set.intersection(*pfam_dict.values()))