filepath='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/'

filename=filepath+'inspector_flye/contig_to_ref.debreak.temp'
allae=open(filename,'r').read().split('\n')[:-1]

#allcontig=open(filepath+'canu/canu_ccs/chr1_hg002.contigs.fasta','r').read().split('>')[1:]
#allcontig=open(filepath+'flye/assembly.fasta','r').read().split('>')[1:]
allcontig=open(filepath+'flye/assembly.fasta','r').read().split('>')[1:]

contiglen={}
for c in allcontig:
	length=0
	ctgname=c.split('\n')[0].split(' ')[0]
	for m in c.split('\n')[1:-1]:
		length+=len(m)
	contiglen[ctgname]=length

bb=[]
for c in allae:
	ctg=c.split('\t')[4]
	if contiglen[ctg]>=1000000:
		bb+=[c]



a=open('/data/user/maggic/svstudy/data/reference/hg38.fa','r').read().split('>')[1:]
refse={}
for c in a:
	name=c.split('\n')[0].split(' ')[0]
	seq=''
	for mm in c.split('\n')[1:-1]:
		seq+=mm
	refse[name]=seq

dd=[]
for c in bb:
	if 'D-' in c:
		seq=refse[c.split('\t')[0]][int(c.split('\t')[1])-500:int(c.split('\t')[1])+500+int(c.split('\t')[2])]
		if 'NNN' not in seq:
			dd+=[c]
		continue
	if 'I-' in c:
		seq=refse[c.split('\t')[0]][int(c.split('\t')[1])-500:int(c.split('\t')[1])+500]
		if 'NNN' not in seq:
			dd+=[c]

f=open(filename+'.filtered','w')
for c in dd:
	f.write(c+'\n')
f.close()
