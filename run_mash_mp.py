import subprocess,os,glob,itertools
import os,sys
from glob import glob
from itertools import combinations, zip_longest
import multiprocessing as mp
from math import ceil
import timeit
import concurrent.futures
from pickle import load
#from multiprocessing.pool import ThreadPool

def runIQtree(path):
    infile,outfile = path
    cmd = ["iqtree", "-quiet", "-nt", '1', "-s", infile, "-m", 'GTR', "-pre", outfile,'-bb','1000']
    subprocess.call(cmd)

def runMash(genome1,genome2):
    output = subprocess.check_output(['mash dist -s 10000 {} {}'.format(genome1, genome2)],
                shell=True)
    output = output.decode('utf-8')
    return output  # just to simulate some calculation

def runMashChunks(chunk):
    results = []
    for pair in chunk:
        if pair:
            genome1,genome2 = pair
            results.append(runMash(genome1,genome2))
    return(results)

# def feed(queue, taskList):
#     for pair in taskList:
#             queue.put(pair)
#
# def calc(queueIn, queueOut):
#     while True:
#         try:
#             pair = queueIn.get(block = False)
#             #print("dealing with {}".format(pair))
#             result = runMash(pair)
#             queueOut.put((pair,result))
#         except:
#             break
#
# def write(queue, fname):
#     with open(fname, "a") as fhandle:
#         while True:
#             try:
#                 par, result = queue.get(block = True)
#                 print(result)
#                 fhandle.write(result)
#             except:
#                 break
def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

if __name__ == "__main__":
#
#     mashFileSource = '/home/u1474301/mash_calc/actino_sketches'
#     outputDir = '/home/u1474301/mash_calc/'
#     mashFiles = glob('{}/*msh'.format(mashFileSource))
#     pairs = itertools.combinations(mashFiles, 2)
#     numThreads = mp.cpu_count()
#     chunkSize = 10000
#     chunks = grouper(chunkSize, pairs)
#     outputName = "distances.tsv"
    acc2asmDict = load(open('acc2asm.pkl','rb'))
    clusters = glob('/data/Groups/Theme2/id_ripp_clustering/*.pkl')
    mashFiles = glob('/data/Groups/Theme2/rep_genomes_fna/*.msh')
    genomesSketched = set(os.path.splitext(x)[0] for x in mashFiles)
    for cluster in clusters:
        groupDict = load(open(cluster,'rb'))
        aveDistDict = {}
        clusterName = os.path.split(cluster)[1].split('.')[0]
        print(clusterName)
        for groupNum,codingGroups in groupDict.items():
            genomesToCompare = []
            for codingGroup in codingGroups:
                ncbiAcc = '_'.join(codingGroup.split('_')[:-1])
                if ncbiAcc in acc2asmDict:
                    accBase = acc2asmDict[ncbiAcc].split('.')[0]
                    path = glob('/data/Groups/Theme2/rep_genomes_fna/{}*.gz'.format(accBase))
                    if path:
                        genomesToCompare.append(path[0])

            pairs = combinations(genomesToCompare,2)
            pairs = list(pairs)
            num_pairs = len(pairs)
            print(num_pairs)
            total_distance = 0
            for genome1,genome2 in pairs:
                if genome1 not in genomesSketched:
                    subprocess.check_output(['mash sketch -s 10000 {}'.format(genome1)],
                                            shell=True)
                    genomesSketched.add(genome1)


                if genome2 not in genomesSketched:
                    subprocess.check_output(['mash sketch -s 10000 {}'.format(genome2)],
                                            shell=True)
                    genomesSketched.add(genome2)

                output = runMash(genome1+'.msh',genome2+'.msh')
                total_distance += float(output.split('\t')[2])
            if num_pairs > 0:
                ave_dist = total_distance/num_pairs
                aveDistDict[groupNum] = (ave_dist,num_pairs)
        print(aveDistDict)
        with open('/data/Groups/Theme2/id_ripp_clustering/{}.avedist.csv'.format(clusterName),'w') as outfile:
            for groupNum,(ave_dist,size) in aveDistDict.items():
                outfile.write('{},{},{}\n'.format(groupNum,size,ave_dist))


    #
    #
    #
    # alignFiles = glob('/home/u1474301/uplb/mlst_genes/hahella/trim/*fasta')
    #
    # outdir = '/home/u1474301/uplb/mlst_genes/hahella/trees'
    # paths = []
    #
    # for alignment in alignFiles:
    #     baseName = os.path.splitext(os.path.split(alignment)[1])[0]
    #     outpath = outdir + baseName
    #     paths.append((alignment,outpath))
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     futures = [executor.submit(runIQtree, path) for path in paths]
    #     concurrent.futures.wait(futures)
    #     # for future in concurrent.futures.as_completed(futures):
    #     #     result = future.result()
    #     #     with open(os.path.join(outputDir, outputName), 'a') as outfile:
    #     #         outfile.writelines(result)
    #     #     print('Finished Chunk {} of {}'.format(chunkNum, len(futures)))
    #     #     chunkNum += 1
    #
    #
    # def multi():
    #     # mashFileSource = '/home/u1474301/mash_calc/actino_sketches'
    #     # outputDir = '/home/u1474301/mash_calc/'
    #     mashFileSource = '/Users/emzodls/Downloads/genome_assemblies/ncbi-genomes-2018-02-06/'
    #     outputDir = '/Users/emzodls/Downloads/genome_assemblies/'
    #     mashFiles = glob('{}/*msh'.format(mashFileSource))[:50]
    #     # mashFileNames = set(os.path.split(x)[1] for x in mashFiles)
    #     pairs = itertools.combinations(mashFiles, 2)
    #     numThreads = mp.cpu_count()
    #     fname_multi = "test_multi.tsv"
    #     fname_single = "test.tsv"
    #     chunkSize = 100
    #     chunks = grouper(chunkSize, pairs)
    #     pool = mp.Pool(numThreads,maxtasksperchild=10)
    #     with concurrent.futures.ProcessPoolExecutor() as executor:
    #         chunkNum = 1
    #         futures = [executor.submit(runMashChunks,chunk) for chunk in chunks]
    #         for future in concurrent.futures.as_completed(futures):
    #             result = future.result()
    #             with open(os.path.join(outputDir, fname_multi), 'a') as outfile:
    #                 outfile.writelines(result)
    #             print('Finished Chunk {} of {}'.format(chunkNum, len(futures)))
    #             #task.add_done_callback(write)
    #             #dist = pool.apply_async(runMashChunks,args=[chunk])
    #             #result = task.result()
    #             # with open(os.path.join(outputDir,fname_multi),'a') as outfile:
    #             #     outfile.writelines(result)
    #             chunkNum += 1
    #
    # def single():
    #     # mashFileSource = '/home/u1474301/mash_calc/actino_sketches'
    #     # outputDir = '/home/u1474301/mash_calc/'
    #     mashFileSource = '/Users/emzodls/Downloads/genome_assemblies/ncbi-genomes-2018-02-06/'
    #     outputDir = '/Users/emzodls/Downloads/genome_assemblies/'
    #     mashFiles = glob('{}/*msh'.format(mashFileSource))[:50]
    #     # mashFileNames = set(os.path.split(x)[1] for x in mashFiles)
    #     pairs = itertools.combinations(mashFiles, 2)
    #
    #     numThreads = mp.cpu_count()
    #     fname_multi = "test_multi.tsv"
    #     fname_single = "test.tsv"
    #     chunkSize = 100
    #     chunks = grouper(chunkSize, pairs)
    #     for chunk in chunks:
    #         with open(os.path.join(outputDir, fname_single), 'a') as outfile:
    #             outfile.writelines(runMashChunks(chunk))
    #             print('SpawnThread')
    # #
    #
    # t1 = timeit.Timer('multi()', "from __main__ import multi")
    #
    # t2 = timeit.Timer('single()', "from __main__ import single")
    #
    # print('Multi:',t1.timeit(2),'Single:',t2.timeit(2))
    #


# def grouper(n, iterable, padvalue=None):
#     "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
#     return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)
#
#
# def processChunk(chunk,chunkIdx,outpath):
#     outfile = open('{}/distances_{}.tsv'.format(outpath,chunkIdx),'a')
#     for idx, pair in enumerate(chunk):
#         if pair:
#             #print('number {} on chunk {}'.format(idx,chunkIdx))
#             genome1,genome2 = pair
#             subprocess.Popen(
#                 ['mash dist {} {}'.format(genome1, genome2)],
#                 stdout=outfile, shell=True)