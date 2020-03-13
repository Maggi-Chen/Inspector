

#filepath='/data/scratch/maggic/denovo_evaluation/hg002_clr/'
#contig_to_ref_sam=filepath+'inspector_flye/contig_to_ref.sam'
#svfile=filepath+'inspector_flye/assembly_errors.bed-gt'
#contig=open(filepath+'flye/assembly.fasta','r').read().split('>')[1:]


filepath='/scratch/maggic/denovo_evaluation/hg002_ccs/racon_polish/inspector_wtdbg_racon2/'
contig_to_ref_sam=filepath+'contig_to_ref.sam'
#contig_to_ref_sam='/data/scratch/maggic/simu/assembly_errors.bed-gt'
svfile=filepath+'assembly_errors.bed-gt'
contig=open('/scratch/maggic/denovo_evaluation/hg002_ccs/racon_polish/wtdbg_racon_polish2.fasta','r').read().split('>')[1:]

#contig_to_ref_sam='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/inspector_flye/contig_to_ref.sam'
#svfile='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero//inspector_flye/assembly_errors.bed-gt'
#contig=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero//flye/assembly.fasta','r').read().split('>')[1:]



def cigarinfo(cigar):
	num=''
	matchlen=0
	info=[0,0,0]
	for c in cigar:
		if c in '1234567890':
			num+=c
		if c in 'MI=X':
			matchlen+=int(num); num=''
		if c in 'D':
			num=''
		if c in 'SH':
			if matchlen==0:
				info[0]=int(num); num=''
			else:
				info[2]=int(num); num=''
	info[1]=matchlen
	return info

def findcontigpos(cigar,start,pos,clipinfo,flag,size,svtype):
	num=''
	if int(flag) in [0,2048]:
		target=pos-clipinfo[0]
	else:
		target=clipinfo[1]-pos+clipinfo[2]
	location=int(start)
	for c in cigar:
		if c in '1234567890':
			num+=c; continue
		if c in 'M=X':
			if target<int(num):
				location+=target; return location
			else:
				target-=int(num); location+=int(num); num=''
		if c in 'D':
			if c==svtype and ((target<=500 and 0.8<=float(num)/size<=1.2 ) or (size>=1000 and 0.9<=float(num)/size<=1.1 and target<=1000)):
				print location
				return location
			location+=int(num); num='';
		if c in 'I':
			if c==svtype and ((target<=500 and 0.8<=float(num)/size<=1.2) or (size>=1500 and 0.9<=float(num)/size<=1.1 and target<=1000)):
				return location
			target-=int(num)
			num=''
		if c in 'SH':
			num=''
		if target<=0:
			return location


f=open(contig_to_ref_sam,'r')
a=f.readline()
while a[0]=='@':
	a=f.readline()

contigs=[]

while a!='':
	if a.split('\t')[1] in ['256','272']:
		a=f.readline(); continue
	contigname=a.split('\t')[0]
	ref=a.split('\t')[2]
	refstart=a.split('\t')[3]
	cigar=a.split('\t')[5]
	flag=int(a.split('\t')[1])
	clipinfo=cigarinfo(cigar)
	if a.split('\t')[1] in ['16','2064']:
		contigstart=clipinfo[2]
		contigstop=clipinfo[1]+contigstart
	else:
		contigstart=clipinfo[0]
		contigstop=clipinfo[0]+clipinfo[1]

	contigs+=[[contigname,contigstart,contigstop,ref,refstart,cigar,clipinfo,flag]]
	a=f.readline()

f.close()


contigs=[c for c in contigs if c[2]-c[1]>=20000]



f=open(filepath+'contigmapinfo','w')
for c in contigs:
	f.write(str(c)+'\n')
f.close()
quit()


large_contigs=[]
for c in contig:
	c=c.split('\n')[:-1]
	length=0
	for cc in c:
		length+=len(cc)
	if length>1000000:
		large_contigs+=[c[0].split(' ')[0]]


allsv=open(svfile,'r').read().split('\n')[:-1]
allsvpos=[]
iii=0
small_contig=0

for c in allsv:
	contigname=c.split('\t')[0]
	if contigname not in large_contigs:
		small_contig+=1
		continue
	contigpos=int(c.split('\t')[1])
	contiginfo=[]
	size=''
	svtype=''
	if 'Expansion' in c:
		size=str(int(c.split('\t')[2])-int(c.split('\t')[1]))
		svtype='I'
	if 'Collapse' in c:
		size=c.split('\t')[5].split('=')[1]
		svtype='D'
	if 'Inversion' in c:
		size=str(int(c.split('\t')[2])-int(c.split('\t')[1]))
		svtype='None'
	for d in contigs:
		if d[0]==contigname and d[1]<contigpos<d[2]:
			contiginfo=d
			break
	if contiginfo!=[]:
		refpos=findcontigpos(contiginfo[5],contiginfo[4],contigpos,contiginfo[6],contiginfo[7],float(size),svtype)
		#allsvpos+=[c+'\t'+d[3]+'\t'+str(refpos)+'\t'+size]
		allsvpos+=[contiginfo[3]+'\t'+str(refpos)+'\t'+size+'\t'+c.split('\t')[3]+'\t'+c]

	else:
		#allsvpos+=[c+'\tContigNotMapped']
		iii+=1
'''

for d in contigs:
	if d[0]=='ctg4' and d[1]<1580658<d[2]:
		contiginfo=d;break
refpos=findcontigpos(d[5],d[4],1580658,d[6],d[7])
print refpos
'''




print str(iii)+' svs are not mapped'
print str(small_contig)+' SVs are in small contigs'

f=open(svfile+'.refpos','w')
for c in allsvpos:
	f.write(c+'\n')
f.close()










