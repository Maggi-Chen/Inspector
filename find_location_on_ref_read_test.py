


import sys 
name=sys.argv[1]
svfile='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/test_aligenr/flye_largekmer/'+name



contig_to_ref_sam='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/inspector_flye_polished/contig_to_ref.sam'
#svfile='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/test_aligenr/test_ngmlr.wtdbgp17.read_ctg.sam'
#svfile='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/inspector_wtdbgp17/read_to_contig.sam'
contig=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/flye/scaffolds.fasta','r').read().split('>')[1:]



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

def findcontigpos(cigar,start,pos,clipinfo,flag):
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
			location+=int(num);  num='';
		if c in 'I':
			target-=int(num) 
			num=''
		if c in 'SH':
			num=''
		if target<=0:
			return location
	print 'Wrong' 
	quit()

def process_readinfo(readinfo):
	
	refstart=int(readinfo.split('\t')[3])
	cigar=readinfo.split('\t')[5]
	clip=[0,0,0]
	num=''
	readlen=0
	refend=refstart
	for  c in cigar:
		if c in '1234567890':
			num+=c; continue
		if c in 'M=X':
			refend+=int(num)
			readlen+=int(num)
			num=''; continue
		if c in 'D':
			refend+=int(num)
			num=''; continue
		if c in 'I':
			readlen+=int(num)
			num=''; continue
		if c in 'SH':
			if readlen==0:
				clip[0]=int(num)
			else:
				clip[2]=int(num)
			num=''; continue
	clip[1]=readlen
	if readinfo.split('\t')[1] in ['0','2048']:
		readstart=clip[0]+1
		readend=clip[0]+clip[1]
	else:
		readstart=clip[2]+1
		readend=clip[2]+clip[1]

	return [readstart,readend,refstart,refend]

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

'''
f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/test_aligenr/flye_minimap2.contig.info','w')
for c in contigs:
	f.write(str(c)+'\n')
f.close()
'''





large_contigs=[]
for c in contig:
	c=c.split('\n')[:-1]
	length=0
	for cc in c:
		length+=len(cc)
	if length>1000000:
		large_contigs+=[c[0].split(' ')[0]]


f=open(svfile,'r')
g=open(svfile+'_refpos','w')
readinfo=f.readline()
iii=0
small_contig=0

while readinfo!='':
	if readinfo[0]=='@' or readinfo.split('\t')[1] in ['256','272']:
		readinfo=f.readline()
		continue
	if readinfo.split('\t')[1] in ['4']:
		iii+=1
		readinfo=f.readline()
		continue

	if readinfo.split('\t')[1] in ['0','2048']:
		reverse=False
	else:
		reverse=True

	contigname=readinfo.split('\t')[2]
	if contigname not in large_contigs:
		if readinfo.split('\t')[1] in ['0','16']:
			small_contig+=1
		readinfo=f.readline()
		continue

	processedinfo=process_readinfo(readinfo)
	

	contiginfo=[]

	for d in contigs:
		if d[0]==contigname and d[1]<processedinfo[2] and processedinfo[3] <d[2]:
			contiginfo=d
			break


	if contiginfo!=[]:
		if reverse:
			refpos=findcontigpos(contiginfo[5],contiginfo[4],processedinfo[3],contiginfo[6],contiginfo[7])
		else:
			refpos=findcontigpos(contiginfo[5],contiginfo[4],processedinfo[2],contiginfo[6],contiginfo[7])

		if contiginfo[7] in [0,2048]:
			ctgrev=False
		else:
			ctgrev=True
		if reverse == ctgrev:
			g.write(readinfo.split('\t')[0]+'\t+\t'+str(processedinfo[0])+'\t'+str(processedinfo[1])+'\t'+contiginfo[3]+'\t'+str(refpos)+'\n')
		else:
			g.write(readinfo.split('\t')[0]+'\t-\t'+str(processedinfo[0])+'\t'+str(processedinfo[1])+'\t'+contiginfo[3]+'\t'+str(refpos)+'\n')


	else:
		iii+=1

	readinfo=f.readline()

'''

for d in contigs:
	if d[0]=='ctg4' and d[1]<1580658<d[2]:
		contiginfo=d;break
refpos=findcontigpos(d[5],d[4],1580658,d[6],d[7])
print refpos
'''




print str(iii)+' svs are not mapped flye'
print str(small_contig)+' SVs are in small contigs'

'''
f=open(svfile+'.refpos','w')
for c in allsvpos:
	f.write(c+'\n')
f.close()
'''












