f=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1/inspector_wtdbgp17/contig_to_ref.debreak.temp','r')
a=f.read().split('\n')[:-1]
f.close()
f=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1/chr1_1.fa','r')
ref=f.read().split('>')[1:]
f.close()

refinfo={}
for  c in ref:
	seq=''
	c=c.split('\n')[:-1]
	refname=c[0].split(' ')[0]
	for cc in c[1:]:
		seq+=cc
	refinfo[refname]=seq

contig=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1/wtdbgp17/chr1_diploid.ctg.fa','r').read().split('>')[1:]
contigs=[]
for c in contig:
	c=c.split('\n')[:-1]
	length=0
	for cc in c:
		length+=len(cc)
	if length>1000000:
		contigs+=[c[0].split(' ')[0]]

f=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1/inspector_wtdbgp17/contig_to_ref.debreak.temp_filtered','w')
for c in a:
	if c.split('\t')[4] not in contigs:
		continue
	if 'D-' in c or 'INV-' in c:
		start=int(c.split('\t')[1])
		end=int(c.split('\t')[2])+start
		seq=refinfo[c.split('\t')[0]][start:end]
		if seq.count('N')>0.5*len(seq):
			continue
		f.write(c+'\n')
	if 'I-' in c:
		start=int(c.split('\t')[1])
		seq=refinfo[c.split('\t')[0]][start-100:start+100]
		if seq.count('N')>30:
			continue
		f.write(c+'\n')
f.close()

