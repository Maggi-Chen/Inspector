snp=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/makeref/sv_snp','r').read().split('\n')[:-1]
path='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/inspector_wtdbgp17/'
allerror=open(path+'assembly_basepair_error_ref','r').read().split('\n')[:-1]

b=[]
same=0
svs=[]
for c in allerror:
	if 'SNP' not in c:
		b+=[c]; continue
	ifover=0
	for d in snp:
		if d.split('\t')[0]==c.split('\t')[0] and abs(int(d.split('\t')[1])-int(c.split('\t')[1]))<=100:
			same+=1; ifover=1; svs+=[d]; break
	if ifover==0:
		b+=[c]

f=open(path+'assembly_basepair_error_ref_groundtruth','w')
for c in b:
	f.write(c+'\n')
f.close()

print 'overlap with snp set:  '+str(same)
print 'homo:  '+str(len([c for c in svs if c.split('\t')[4]=='3']))
