import os
from collections import defaultdict
from pickle import dump
from Bio import SeqIO

os.chdir('/Users/emzodls/Dropbox/Lab/Warwick/RiPP_nnets')
cluster_assignments = defaultdict(set)
cluster2group = dict()
for line in open('node_spec_filt_multilevel.csv'):
    line = line.strip()
    lineParse = line.split(',')
    cluster_assignments[lineParse[1]].add(lineParse[0])
    cluster2group[lineParse[0]] = lineParse[1]

counts = [(k,len(v)) for k,v in cluster_assignments.items()]
counts =sorted(counts,key=lambda x:x[1],reverse=True)

os.chdir('/Volumes/lab_data/Ripps/ripp_codingGroups')

cluster_ripp_hits = defaultdict(set)

for line in open('ripp_doms.domtbl'):
    line = line.strip()
    if line.startswith('#'):
        cluster = line[1:]
    else:
        lineParse = line.split()
        gene = lineParse[0]
        ripp_hit = lineParse[3]
        evalue = float(lineParse[12])
        if evalue <= 1e-5:
            cluster_ripp_hits[cluster].add((gene,ripp_hit))

cluster_doms = {k:set(x[1] for x in v) for k,v in cluster_ripp_hits.items()}

ripp_classes = defaultdict(list)

score_cutoff_dict = dict()

for line in open('/Users/emzodls/antismash/antismash/generic_modules/hmm_detection/hmmdetails.txt'):
    line = line.strip()
    lineParse = line.split('\t')
    score_cutoff_dict[lineParse[0]] = float(lineParse[-2])

antismash_ripps = defaultdict(set)

#parse ripp_domtbl based on antismash cutoff scores
for line in open('ripp_domains_antismash.domtbl'):
    line = line.strip()
    if line.startswith('#'):
        cluster = line[1:].split('.antismash')[0]
    else:
        lineParse = line.split()
        gene = lineParse[0]
        ripp_hit = lineParse[3]
        score = float(lineParse[7])
        if score >= score_cutoff_dict[ripp_hit]:
            antismash_ripps[cluster].add((gene,ripp_hit))

    # line = line.strip()
    # if not line.startswith('#'):
    #     lineParse = line.split()
    #     name = lineParse[0]
    #     cluster,gene = name.split('|')
    #     ripp_hmm = lineParse[3]
    #     score = float(lineParse[7])
    #     if score >= score_cutoff_dict[ripp_hmm]:
    #         antismash_ripps[cluster].add((gene,ripp_hmm))

antismash_cluster_doms = {k:set(x[1] for x in v) for k,v in antismash_ripps.items()}

antismash_ripp_classes = defaultdict(list)

for cluster,doms in antismash_cluster_doms.items():
    if ('LANC_like' in doms and 'Flavoprotein' in doms) or ('LANC_like' in doms and'Trp_halogenase' in doms) or ('LANC_like' in doms and 'p450' in doms) or \
        ('LANC_like' in doms and 'Pkinase' in doms) or ('LANC_like' in doms and 'DUF4135' in doms) or ('LANC_like' in doms and 'Lant_dehydr_N' in doms) or \
        ('LANC_like' in doms and 'Lant_dehydr_C' in doms) or ('LANC_like' in doms and 'adh_short' in doms) or ('LANC_like' in doms and 'adh_short_C2' in doms) or \
        'TIGR03731' in doms or 'Antimicr18' in doms or 'Gallidermin' in doms or 'L_biotic_A' in doms or 'leader_d' in doms or \
        'leader_eh' in doms or 'leader_abc' in doms or 'mature_d' in doms or 'mature_ab' in doms or 'mature_a' in doms or 'mature_b' in doms or \
        'mature_ha' in doms or 'mature_h_beta' in doms or 'lacticin_l' in doms or 'lacticin_mat' in doms or 'LD_lanti_pre' in doms or \
        'strep_PEQAXS' in doms:

        antismash_ripp_classes['lanti'].append(cluster)

    elif 'strepbact' in doms or 'Antimicrobial14' in doms or 'Bacteriocin_IId' in doms or 'BacteriocIIc_cy' in doms or 'Bacteriocin_II' in doms or \
        'Bacteriocin_IIi' in doms or 'Lactococcin' in doms or 'Antimicrobial17' in doms or 'Lactococcin_972' in doms or 'Bacteriocin_IIc' in doms or \
        'LcnG-beta' in doms or 'Cloacin' in doms or 'Linocin_M18' in doms or 'TIGR03603' in doms or 'TIGR03604' in doms or 'TIGR03605' in doms or \
        'TIGR03651' in doms or 'TIGR03678' in doms or 'TIGR03693' in doms or 'TIGR03798' in doms or 'TIGR03882' in doms or 'TIGR03601' in doms or \
        'TIGR03602' in doms or 'TIGR03795' in doms or 'TIGR03793' in doms or 'TIGR03975' in doms or 'DUF692' in doms or 'TIGR01193' in doms:

        antismash_ripp_classes['bacteriocin'].append(cluster)

    elif ('Lant_dehydr_C' in doms and 'YcaO' in doms) or ('Lant_dehydr_N' in doms and 'YcaO' in doms) or 'thiostrepton' in doms or \
        ('YcaO' in doms and 'PF06968' in doms and 'thio_amide' in doms) or ('YcaO' in doms and 'PF04055' in doms and 'thio_amide' in doms) or \
        ('YcaO' in doms and 'PF07366' in doms and 'thio_amide' in doms) or ('YcaO' in doms and 'PF07366' in doms) or \
        ('YcaO' in doms and 'PF06968' in doms) or ('YcaO' in doms and 'PF04055' in doms):

        antismash_ripp_classes['thio'].append(cluster)

    elif 'cypemycin' in doms and 'cypI' in doms:
        antismash_ripp_classes['lina'].append(cluster)

    elif 'cyanobactin_synth' in doms:
        antismash_ripp_classes['cyano'].append(cluster)

    elif 'glycocin' in doms or 'sublancin' in doms:
        antismash_ripp_classes['glycocin'].append(cluster)

    elif 'goadsporin_like' in doms:
        antismash_ripp_classes['lap'].append(cluster)

    elif 'PF13471' in doms and 'PF00733' in doms:
        antismash_ripp_classes['lasso'].append(cluster)

    elif 'subtilosin' in doms or 'thuricin' in doms or 'TIGR04404' in doms or 'TIGR03973' in doms:
        antismash_ripp_classes['sacti'].append(cluster)

    elif 'botH' in doms:
        antismash_ripp_classes['bottro'].append(cluster)

    elif 'Subtilosin_A' in doms or 'skfc' in doms:
        antismash_ripp_classes['htt'].append(cluster)

    elif 'micJ25' in doms or 'mcjC' in doms:
        antismash_ripp_classes['micro'].append(cluster)

    elif 'mvd' in doms or 'mvnA' in doms:
        antismash_ripp_classes['mvd'].append(cluster)

    elif 'PoyD' in doms:
        antismash_ripp_classes['prot'].append(cluster)

for cluster,doms in cluster_doms.items():
    if 'AgrB' in doms and 'AgrD' in doms:
        ripp_classes['aip'].append(cluster)
    elif 'DUF95' in doms and 'HTT_precursor' in doms:
        ripp_classes['htt'].append(cluster)
    elif 'BotA' in doms and 'BotC' in doms:
        ripp_classes['bot'].append(cluster)
    elif 'ComQ' in doms and 'ComX' in doms:
        ripp_classes['comX'].append(cluster)
    elif 'PatA' in doms and 'PatE' in doms and ('PatG' in doms or 'PatG_ox' in doms):
        ripp_classes['cyano'].append(cluster)
    elif 'SunA' in doms and 'SunS' in doms:
        ripp_classes['gly'].append(cluster)
    elif 'lanB' in doms and 'lanC' in doms:
        ripp_classes['lanI'].append(cluster)
    elif 'lanM' in doms:
        ripp_classes['lanII'].append(cluster)
    elif 'lanKC' in doms:
        ripp_classes['lanIII'].append(cluster)
    elif ('CypA' in doms or 'LegA' in doms) and ('CypH' in doms or 'LegH' in doms)  and 'CypL' in doms :
        ripp_classes['lina'].append(cluster)
    elif (('McbC' in doms or 'McbD' in doms) and 'McbB' in doms) or 'GodG' in doms:
        ripp_classes['lap'].append(cluster)
    elif 'MdnA' in doms and ('MdnB' in doms or 'MdnC' in doms):
        ripp_classes['micro'].append(cluster)
    elif 'ProcA' in doms:
        ripp_classes['prochl'].append(cluster)
    elif 'PoyA' in doms and 'PoyD' in doms:
        ripp_classes['prot'].append(cluster)
    elif 'SboA' in doms and 'AlbA' in doms:
        ripp_classes['sacti'].append(cluster)
    elif 'StrA' in doms and 'StrB' in doms and 'StrC' in doms:
        ripp_classes['sacti'].append(cluster)
    elif 'LazA' in doms and 'LazB' in doms and 'LazC' in doms:
        ripp_classes['thio'].append(cluster)
    elif 'TfxA' in doms and 'TfxB' in doms and 'TfxB' in doms:
        ripp_classes['trifo'].append(cluster)
    elif 'TvaA' in doms and 'TvaH' in doms:
        ripp_classes['thioviri'].append(cluster)
    elif 'YmA' in doms and 'YmF' in doms:
        ripp_classes['ym'].append(cluster)


antismash_ripp_cts = {k:len(v) for k,v in antismash_ripp_classes.items()}
ripp_classes_counts = {k:len(v) for k,v in ripp_classes.items()}

family_classes = defaultdict(set)
no_family = list()
for key,clusters in antismash_ripp_classes.items():
    for cluster in antismash_ripp_classes[key]:
        if cluster in cluster2group:
            family = cluster2group[cluster]
            family_classes[family].add('antismash_{}'.format(key))
        else:
            no_family.append(cluster)

for key,clusters in ripp_classes.items():
    for cluster in ripp_classes[key]:
        if cluster in cluster2group:
            family = cluster2group[cluster]
            family_classes[family].add('prism_{}'.format(key))
        else:
            no_family.append(cluster)
# families to consider
families_to_consider = []
clusters_to_consider = set()

identified_clusters = []

for family in counts:
    if len(family_classes[family[0]]) == 0 and family[1] >= 100:
        families_to_consider.append(family)
    elif len(family_classes[family[0]]) > 0 and family[1] >= 100:
        identified_clusters.append(family)

print(len(identified_clusters))
for family,count in identified_clusters:
    print(family,family_classes[family])

dump(identified_clusters,open('identified_families.pkl','wb'))
# print(len(families_to_consider))
#
for family in identified_clusters:
    for cluster in cluster_assignments[family[0]]:
        clusters_to_consider.add(cluster)

print(len(clusters_to_consider))

coding_group_seqs = SeqIO.parse('codingGroups.fa','fasta')

for entry in coding_group_seqs:
    cluster = entry.id.split('|')[0]
    if cluster in clusters_to_consider:
        with open('identified_families/{}.fa'.format(cluster2group[cluster]),'a') as handle:
            SeqIO.write(entry,handle,'fasta')
