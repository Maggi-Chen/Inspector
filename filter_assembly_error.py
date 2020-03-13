filepath='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_wtdbgp17/'


def filter_1():
	a=open(filepath+'assembly_errors.bed-gt.refpos','r').read().split('\n')[:-1]
	expinfo=open(filepath+'deletion-merged','r').read().split('\n')[:-1]
	colinfo=open(filepath+'insertion-merged','r').read().split('\n')[:-1]

	b=[]
	for c in a:
		if 'ContigNotMapped' in c:
			continue
		# Test Function
		if int(c.split('\t')[5])<=10000:
			continue
		if 'Exp' in c:
			for d in expinfo:
				if c.split('\t')[4]==d.split('\t')[0] and c.split('\t')[5]==d.split('\t')[1] and c.split('\t')[2]==d.split('\t')[2]:
					if float(d.split('\t')[4])>=10:
						b+=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[3]+'\t'+d.split('\t')[4]+'\tExpansion']
					break
		if 'Col' in c:
			for d in colinfo:
				if c.split('\t')[4]==d.split('\t')[0] and c.split('\t')[5]==d.split('\t')[1] and c.split('\t')[2]==d.split('\t')[2]:
					if float(d.split('\t')[4])>=10 :
						b+=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[3]+'\t'+d.split('\t')[4]+'\tCollapse']
					break
	f=open(filepath+'assembly_errors.bed-gt.refpos_filtered','w')
	for c in b:
		f.write(c+'\n')
	f.close()


def filter_2():
	tools='flye'
	#filepath='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_'+tools+'/'
	filepath='/data/scratch/maggic/denovo_evaluation/hg002_ccs/inspector_'+tools+'/'


	a=open(filepath+'assembly_errors.bed-gt','r').read().split('\n')[:-1]
	ctgs=open(filepath+'valid_contig.fa','r').read().split('>')[1:]
	ctglen={}
	for c in ctgs:
		name=c.split('\n')[0]
		length=len(c.split('\n')[1])
		ctglen[name]=length

	large_ctg=[]
	for c in ctglen:
		if ctglen[c]>=1000000:
			large_ctg+=[c]
	b=[]
	for c in a:
		pos1=int(c.split('\t')[1])
		pos2=int(c.split('\t')[2])
		ctg=ctglen[c.split('\t')[0]]
		if pos1<=10000 or (ctg-pos2)<=10000:
			continue
		b+=[c]

	large=[c for c in b if c.split('\t')[0] in large_ctg]
	f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/assembly_error/ccs_'+tools+'.ae','w')
	for c in b:
		f.write(c+'\n')
	f.close()
	f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/assembly_error/ccs_'+tools+'_large.ae','w')
	for c in large:
		f.write(c+'\n')
	f.close()






filter_2()

'''

def sortae(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

filepath='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_wtdbgp17/'
#filename='contig_to_ref.debreak.temp_filtered'
#filename='assembly_errors.bed-gt.refpos_filtered'
#filename='missed_assembly_error'
filename='FP_assembly_error'
allae=open(filepath+filename,'r').read().split('\n')[:-1]
allae.sort(key=sortae)
testif=0

while testif==0:
	testif=1
	print 'repeat'
	for i in range(len(allae)-1):
		print i
		if allae[i].split('\t')[0]!=allae[i+1].split('\t')[0]:
			continue
		
		if 'Exp' in allae[i]:
		#if 'I-' in allae[i]:
			if int(allae[i].split('\t')[1])+int(allae[i].split('\t')[2])+1000>=int(allae[i+1].split('\t')[1]):
				testif=0
		
		if 'Col' in allae[i]:
		#if 'D-' in allae[i]:
			if int(allae[i].split('\t')[1])+1000>=int(allae[i+1].split('\t')[1]):
				testif=0
		
		if int(allae[i].split('\t')[1])+int(allae[i].split('\t')[2])+1000>=int(allae[i+1].split('\t')[1]):
			testif=0
		print testif
		print i
		if testif==0:
			print 'bad'
			if 'Exp' in allae[i]:
			#if 'I-' in allae[i]:
				size1=int(allae[i].split('\t')[2])
			else:
				size1=-int(allae[i].split('\t')[2])
			if 'Exp' in allae[i+1]:
			#if 'I-' in allae[i+1]:
				size2=int(allae[i+1].split('\t')[2])
			else:
				size2=-int(allae[i+1].split('\t')[2])
			size=size1+size2
			mergetype=''
			if size>50:
				mergetype='Expansion'
				#mergetype='I-merge'
			if size<-50:
				mergetype='Collapse'			
				#mergetype='D-merge'
			if mergetype!='':
				allae[i]=allae[i].split('\t')[0]+'\t'+allae[i].split('\t')[1]+'\t'+str(abs(size))+'\t'+mergetype+'\tmerged\t'+allae[i]+'\t'+allae[i+1]
				allae.remove(allae[i+1])
			else:
				allae.remove(allae[i])
				allae.remove(allae[i+1])
			print 'merge  '+str(i)
			break
f=open(filepath+filename+'_merged','w')
for c in allae:
	f.write(c+'\n')
f.close()
'''

