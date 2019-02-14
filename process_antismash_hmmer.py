from glob import glob
taskList = glob('*antismash*domtbl')

with open('../ripp_domains_antismash.domtbl','w') as outfile:
	for idx,domtbl in enumerate(taskList):
		basename = domtbl.split('*.antismash.domtbl')[0]
		outfile.write('#{}\n'.format(basename))
		print('{}/{}\n'.format(idx+1,len(taskList)))
		for line in open(domtbl):
			if not line.startswith('#'):
				outfile.write(line)
