

'''
contig_to_ref_sam='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/inspector_canu/contig_to_ref.sam'
#svfile='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_flye_new/assembly_basepair_error.vcf'
svfile='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/inspector_canu/assembly_basepair_error.vcf'
#contig=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/wtdbg/chr1_diploid.ctg.fa','r').read().split('>')[1:]
contig=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/canu/canu_out/chr1_diploid.contigs.fasta','r').read().split('>')[1:]
#contig=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/flye/assembly.fasta','r').read().split('>')[1:]
'''
contig_to_ref_sam='/data/scratch/maggic/denovo_evaluation/hg002_ccs/inspector_canu/contig_to_ref.sam'
svfile='/data/scratch/maggic/denovo_evaluation/hg002_ccs/inspector_canu/assembly_basepair_error.vcf'
#contig=open('/data/scratch/maggic/denovo_evaluation/hg002_clr/wtdbg/chr1_hg002.ctg.fa','r').read().split('>')[1:]
#contig=open('/data/scratch/maggic/denovo_evaluation/hg002_ccs/flye/assembly.fasta','r').read().split('>')[1:]
contig=open('/data/scratch/maggic/denovo_evaluation/hg002_ccs/canu/canu_ccs/chr1_hg002.contigs.fasta','r').read().split('>')[1:]




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
				info[2]==int(num); num=''
	info[1]=matchlen
	return info

def findcontigpos(cigar,start,pos,clipinfo,flag):
	num=''
	if flag in [0,2048]:
		target=pos-clipinfo[0]
	else:
		target=clipinfo[1]-pos+clipinfo[2]
	location=int(start)
	for c in cigar:
		if c in '1234567890':
			num+=c
		if c in 'M=X':
			if target<int(num):
				location+=target; return location
			else:
				target-=int(num); location+=int(num); num=''
		if c in 'D':
			location+=int(num); num=''
		if c in 'I':
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
	if c[0]=='#':
		continue
	contigname=c.split('\t')[0]
	if contigname not in large_contigs:
		small_contig+=1
		continue
	contigpos=int(c.split('\t')[1])
	contiginfo=[]
	size=''
	svtype=''
	if 'TYPE=del' in c:
		size=str(len(c.split('\t')[3])-len(c.split('\t')[4].split(',')[0]))
		svtype='Expansion'
	if 'TYPE=ins' in c:
		size=str(len(c.split('\t')[4].split(',')[0])-len(c.split('\t')[3]))
		svtype="Collapse"
	if 'TYPE=snp' in c:
		size='1'; svtype="SNP"
	if 'TYPE=mnp' in c:
		size=str(len(c.split('\t')[4])); svtype="MNP"
	if svtype=='' or int(size)>10:
		continue
	for d in contigs:
		if d[0]==contigname and d[1]<contigpos<d[2]:
			contiginfo=d
			break
	if contiginfo!=[]:
		refpos=findcontigpos(d[5],d[4],contigpos,d[6],d[7])
		#allsvpos+=[c+'\t'+d[3]+'\t'+str(refpos)+'\t'+size]
		allsvpos+=[d[3]+'\t'+str(refpos)+'\t'+size+'\t'+svtype+'\t'+c]

	else:
		allsvpos+=[c+'\tContigNotMapped']
		iii+=1

print str(iii)+' svs are not mapped'
print str(small_contig)+' SVs are in small contigs'

f=open(svfile+'.refpos','w')
for c in allsvpos:
	f.write(c+'\n')
f.close()











