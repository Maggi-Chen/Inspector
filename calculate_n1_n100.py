tools='wtdbg'
#filepath='/data/scratch/maggic/denovo_evaluation/hg002_clr/inspector_'+tools+'/'
filepath='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/inspector_'+tools+'/'

f=open(filepath+'valid_contig.fa','r')
a=f.readline()
size=[]
while a!='':
	if '>' in a:
		a=f.readline()
	else:
		size+=[len(a)-1]
		a=f.readline()

size.sort(reverse=True)


totalsize=sum(size)

curr=0
n100=[]
iii=1
for c in size:
	curr+=c
	while curr>=iii*totalsize/100:
		n100+=[c]
		iii+=1

print len(n100)
f=open('/data/scratch/maggic/simulation/denovosimulation/halploid_chr1_noNread/n1_100/sim_hal_'+tools+'.txt','w')
for  i in range(100):
	f.write(str(i+1)+'\t'+str(n100[i])+'\n')
f.close()
