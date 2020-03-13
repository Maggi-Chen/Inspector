tools='canu'
#f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_'+tools+'/contig_to_ref.sam','r')
f=open('/data/scratch/maggic/denovo_evaluation/hg002_clr/inspector_'+tools+'/contig_to_ref.sam','r')

a=f.read().split('\n')[:-1]
f.close()



def calcu_len(cigar):
	num=''
	length=0
	for c in cigar:
		if  c in '1234567890':
			num+=c; continue
		if c in '=XISH':
			length+=int(num); num=''; continue
		if c in 'D':
			num=''; continue
		print c
	return length


def findpos(c,ctglength):
	temppos=[]
	refpos=int(c.split('\t')[3])
	cigar=c.split('\t')[5]
	if cigar.split('M')[0].split('=')[0].split('S')[0].split('H')[0] =='':
		ctgpos=0
	else:
		ctgpos=int(cigar.split('M')[0].split('=')[0].split('S')[0].split('H')[0])
	currpos=ctgpos
	
	num=''
	
	for m in cigar:
		if m in '1234567890':
			num+=m; continue
		if m in 'M=X':
			ctgpos+=int(num); refpos+=int(num); num=''
		if m=='I':
			ctgpos+=int(num);num=''
		if m =='D':
			refpos+=int(num); num='';continue
		if m in 'SH':
			if refpos>int(c.split('\t')[3]):
				break
			else:
				num=''
		while  ctgpos>=currpos+2000:
			if ctgpos>=currpos+4000:
				temppos+=[[currpos+2000,refpos]]
			else:
				temppos+=[[ctgpos,refpos]]
			currpos+=2000
	if c.split('\t')[1] in ['16','2064']:
		newtemp=[]
		for mm in temppos:
			newtemp+=[[ctglength-mm[0],mm[1]]]
		temppos=newtemp
	return temppos

contig_len={}

samectg=''
processedlength=0
samepos=[]
allpos=[]
lastlength=0
for c in a:
	if c[0]=='@' or c.split('\t')[1] in ['256','272','4'] :
		continue

	temppos=[]
	ctgname=c.split('\t')[0]
	if ctgname not in contig_len:
		ctglength=calcu_len(c.split('\t')[5])
		contig_len[ctgname]=ctglength
	else:
		ctglength=contig_len[ctgname]

	if ctglength < 1000000:
		continue

	temppos=findpos(c,ctglength)
	
	if ctgname==samectg:	
		samepos+=temppos
	else:
		samepos.sort()
		if samectg == 'tig00001510':
			f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/contig_linearity/clr_'+tools+'_primary.pos','w')
			for c in samepos:
				f.write(str(c[0])+'\t'+str(c[1])+'\n')
			f.close()

		newpos=[]
		for mm in samepos:
			newpos+=[[mm[0]+processedlength,mm[1]]]
		allpos+=newpos
		processedlength+=lastlength
		lastlength=ctglength
		samepos=[]
		samectg=ctgname
		samepos=temppos


if samepos!=[]:
	samepos.sort()
	newpos=[]
	for mm in samepos:
		newpos+=[[mm[0]+processedlength,mm[1]]]
	allpos+=newpos


f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/contig_linearity/clr_'+tools+'.pos','w')
for c in allpos:
	f.write(str(c[0])+'\t'+str(c[1])+'\n')
f.close()

