import os
import subprocess
def simple(contig,outpath,min_size,min_size_assemblyerror):
	f=open(contig,'r')
	allcontig=f.read().split('>')[1:]
	f.close()
	length=[]
	maxlen=0
	large_contigs=[]
	all_contigs=[]
	totallength=0
	totallength_large=0
	f=open(outpath+'valid_contig.fa','w')
	for c in allcontig:
		c=c.split('\n')[:-1]
		contig_name=c[0].split(' ')[0]
		all_contigs+=[contig_name]
		seq=''
		length1=0
		for cc in c[1:]:
			seq+=cc
			length1+=len(cc)
		length+=[length1]
		if length1>maxlen:
			maxlen=length1
			maxcontig=contig_name
		if length1>=min_size:
			f.write('>'+contig_name+'\n'+seq+'\n')
			totallength+=length1
		if length1>=min_size_assemblyerror:
			large_contigs+=[contig_name]
			totallength_large+=length1

	f.close()

	length.sort(reverse=True)
	f=open(outpath+'summary_statistics','w')
	f.write('Statics of contigs:\n')
	f.write('Number of contigs\t'+str(len(length))+'\n')
	f.write('Total length\t'+str(sum(length))+'\n')
	f.write('Longest contig\t'+str(length[0])+'\n')
	if len(length)>1:
		f.write('Second longest contig length\t'+str(length[1])+'\n')

	iii=0
	total=sum(length)/2
	for c in length:
		iii+=c
		if iii>=total:
			n50=c; break
	f.write('N50\t'+str(n50)+'\n')
	length=[c for c in length if c > 1000000]
	f.write('Number of contigs >1Mbp\t'+str(len(length))+'\n')
	f.write('Total length of contigs >1Mbp\t'+str(sum(length))+'\n')
	iii=0; total=sum(length)/2; n50=0
	for c in length:
		iii+=c
		if iii>total:
			n50=c; break
	f.write('N50 of contigs >1Mbp\t'+str(n50)+'\n\n\n')
	f.close()
	return [all_contigs,large_contigs,totallength,totallength_large,maxcontig,maxlen]

def sort_sv(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

def assembly_info(outpath,chromosomes):
	os.system("cat "+outpath+"del-merged-* > "+outpath+"deletion-merged")
	os.system("cat "+outpath+"ins-merged-* > "+outpath+"insertion-merged")
	os.system("cat "+outpath+"dup-merged-* > "+outpath+"duplication-merged")
	os.system("cat "+outpath+"inv-merged-* > "+outpath+"inversion-merged")
	#os.system("cat "+outpath+"tra-merged-* > "+outpath+"translocation-merged")
	os.system("rm "+outpath+"*-info-* "+outpath+"*-merged-*")

	f=open(outpath+'deletion-merged','r')
	alldel=f.read().split('\n')[:-1]
	f.close()
	allins=open(outpath+'insertion-merged','r').read().split('\n')[:-1]
	alldup=open(outpath+'duplication-merged','r').read().split('\n')[:-1]
	allins+=alldup
	allsv=alldel+allins
	allsv.sort(key=sort_sv)
	f=open(outpath+'assembly_errors.bed','w')
	for c in allsv:
		if 'Del' in c:
			c=c.split('\t')
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tExpansion\n')
		if 'Ins' in c:
			c=c.split('\t')
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+1)+'\t'+c[3]+'\tCollapse\tSize='+c[2]+'\n')
	allinv=open(outpath+'inversion-merged','r').read().split('\n')[:-1]
	for c in allinv:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tInversion\n')
	f.close()

	f=open(outpath+'summary_statistics','a')
	f.write('Number of assembly collapse\t'+str(len(allins))+'\n')
	f.write('Number of assembly expansion\t'+str(len(alldel))+'\n')
	f.write('Number of assembly inversion\t'+str(len(allinv))+'\n')
	alldel=[c for c in alldel if c.split('\t')[0] in chromosomes]
	allins=[c for c in allins if c.split('\t')[0] in chromosomes]
	allinv=[c for c in allinv if c.split('\t')[0] in chromosomes]

	f.write('Number of assembly collapse in large contigs\t'+str(len(allins))+'\n')
	f.write('Number of assembly expansion in large contigs\t'+str(len(alldel))+'\n')
	f.write('Number of assembly inversion in large contigs\t'+str(len(allinv))+'\n\n\n')

	f.close()
	#os.system("rm "+outpath+"*ion-merged")
	return 0

def assembly_info_ref(outpath):
	os.system("cat "+outpath+"del-merged-* > "+outpath+"deletion-merged_ref")
	os.system("cat "+outpath+"ins-merged-* > "+outpath+"insertion-merged_ref")
	os.system("cat "+outpath+"dup-merged-* > "+outpath+"duplication-merged_ref")
	os.system("cat "+outpath+"inv-merged-* > "+outpath+"inversion-merged_ref")
	os.system("rm "+outpath+"*-info-* "+outpath+"*-merged-*")

	f=open(outpath+'deletion-merged_ref','r')
	alldel=f.read().split('\n')[:-1]
	f.close()
	allins=open(outpath+'insertion-merged_ref','r').read().split('\n')[:-1]
	alldup=open(outpath+'duplication-merged_ref','r').read().split('\n')[:-1]
	allins+=alldup
	allsv=alldel+allins
	allsv.sort(key=sort_sv)
	f=open(outpath+'assembly_errors_ref.bed','w')
	for c in allsv:
		if 'Ins' in c or 'Dup' in c:
			c=c.split('\t')
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\tExpansion\n')
		if 'Del' in c:
			c=c.split('\t')
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+1)+'\tCollapse\tSize='+c[2]+'\n')
	allinv=open(outpath+'inversion-merged_ref','r').read().split('\n')[:-1]
	for c in allinv:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\tInversion\n')
	f.close()

	f=open(outpath+'summary_statistics','a')
	f.write('Assembly errors from contig to reference:\n')
	f.write('Number of assembly collapse\t'+str(len(alldel))+'\n')
	f.write('Number of assembly expansion\t'+str(len(allins))+'\n')
	f.write('Number of assembly inversion\t'+str(len(allinv))+'\n\n\n')
	f.close()
	#os.system("rm "+outpath+"*ion-merged")
	return 0


def basepair_error(outpath):
	f=open(outpath+'assembly_basepair_error.vcf','r')
	a=f.readline()
	mismatch=0
	dels=0
	ins=0
	snp=0
	mnp=0
	while a!='':
		if a[0]=='#':
			a=f.readline(); continue
		if 'TYPE=snp' in a:
			mismatch+=1
			snp+=1
		if 'TYPE=ins' in a:
			mismatch+=len(a.split('\t')[4])-len(a.split('\t')[3])
			ins+=1
		if 'TYPE=del' in a:
			mismatch+=len(a.split('\t')[3])-len(a.split('\t')[4])
			dels+=1
		if 'TYPE=mnp' in a:
			mismatch+=len(a.split('\t')[4])
			mnp+=1
		a=f.readline()
	accuracy=1-mismatch/100000.0
	f=open(outpath+'summary_statistics','a')
	f.write('Number of small collapse\t'+str(ins)+'\n')
	f.write('Number of small expansion\t'+str(dels)+'\n')
	f.write('Number of single basepair error\t'+str(snp)+'\n')
	f.write('Number of multiple basepair error\t'+str(mnp)+'\n')
	f.write('Base pair accuracy\t'+str(accuracy)+'\n\n\n')
	return 0

def basepair_error_ref(outpath,largestchr):
	f=open(outpath+'contig_to_ref.sam','r')
	a=f.readline()
	mismatch=0
	totallength=0
	ins=0; dels=0; snp=0
	svs=[]
	while a!='':
		if a[0]=='@':
			a=f.readline(); continue
		if a.split('\t')[0]!=largestchr or a.split('\t')[1] in ['256','272']:
			a=f.readline(); continue
		cigar=a.split('\t')[5]
		num=''
		length=0
		chrom=a.split('\t')[2]
		pos=int(a.split('\t')[3])
		for c in cigar:
			if c in '1234567890':
				num+=c
			if c in 'SH':
				num=''
			if c in 'M=':
				length+=int(num); num=''
			if c == 'X':
				if int(num)==1:
					svs+=[chrom+'\t'+str(pos+length)+'\t'+str(pos+length)+'\tSNP']
				else:
					svs+=[chrom+'\t'+str(pos+length)+'\t'+str(pos+length+int(num)-1)+'\tMNP\tsize='+num]
				length+=int(num); mismatch+=int(num);num=''
				snp+=1
			if c == 'I':
				if int(num)<=10:
					svs+=[chrom+'\t'+str(pos+length)+'\t'+str(pos+length+1)+'\tExpansion\tsize='+num]
					mismatch+=int(num); ins+=1
				num=''
			if c == 'D':
				if int(num)<=10:
					svs+=[chrom+'\t'+str(pos+length)+'\t'+str(pos+length+int(num))+'\tCollapse']
					mismatch+=int(num); dels+=1
				length+=int(num)
				num=''

		totallength+=length
		a=f.readline()
	accuracy=round((1-mismatch/float(totallength))*10000)/10000.0
	f=open(outpath+'summary_statistics','a')
	f.write('Base pair accuracy of longest contig from contig to reference:\n')
	f.write('Number of small assembly collapse\t'+str(dels)+'\n')
	f.write('Number of small assembly extension\t'+str(ins)+'\n')
	f.write('Number of single basepair error\t'+str(snp)+'\n')
	f.write('Base pair accuracy\t'+str(accuracy)+'\n\n\n')
	f.close()
	f=open(outpath+'assembly_basepair_error_ref','w')
	for c in svs:
		f.write(c+'\n')
	f.close()
	return 0




def get_ref_chroms(outpath):
	f=open(outpath+'contig_to_ref.sam','r')
	a=f.readline()
	chroms=[]
	length=0
	longestlen=0
	longestchr=''
	while a[0]=='@':
		if a[:3]=='@SQ':
			chroms+=[a.split('\t')[1].split(':')[1]]
			length+=int(a.split('\t')[2].split(':')[1])
			if int(a.split('\t')[2].split(':')[1])>longestlen:
				longestlen=int(a.split('\t')[2].split(':')[1])
				longestchr=a.split('\t')[1].split(':')[1]
		a=f.readline()
	a=int(subprocess.check_output("awk \'$3==0\' "+outpath+'contig_to_ref.depth | wc -l',shell=True))
	covered=length-a
	return (chroms,length,longestchr,longestlen,covered)


def check_depth_ref(outpath,ref):
	cov0=int(subprocess.check_output("awk \'$3==0\' "+outpath+'contig_to_ref.depth | wc -l',shell=True))
	cov1=int(subprocess.check_output("awk \'$3==1\' "+outpath+'contig_to_ref.depth | wc -l',shell=True))
	cov2=int(subprocess.check_output("awk \'$3==2\' "+outpath+'contig_to_ref.depth | wc -l',shell=True))
	cov3=int(subprocess.check_output("awk \'$3>2\' "+outpath+'contig_to_ref.depth | wc -l',shell=True))

	total=cov0+cov1+cov2+cov3
	refn=int(subprocess.check_output("grep -o 'N' "+ref+" |wc -l",shell=True).split(' ')[0])
	total-=refn
	cov0-=refn
	
	f=open(outpath+'summary_statistics','a')
	f.write('#BP with cov=0:   '+str(cov0)+',  '+str(cov0*100.00/total)+'\n')
	f.write('#BP with cov=1:   '+str(cov1)+',  '+str(cov1*100.00/total)+'\n')
	f.write('#BP with cov=2:   '+str(cov2)+',  '+str(cov2*100.00/total)+'\n')
	f.write('#BP with cov>2:   '+str(cov3)+',  '+str(cov3*100.00/total)+'\n')
	f.close()
	'''
	print '#BP with cov=0:   '+str(cov0)+',  '+str(cov0*100.00/total)
	print '#BP with cov=1:   '+str(cov1)+',  '+str(cov1*100.00/total)
	print '#BP with cov=2:   '+str(cov2)+',  '+str(cov2*100.00/total)
	print '#BP with cov>2:   '+str(cov3)+',  '+str(cov3*100.00/total)
	'''

if __name__ == '__main__':
	outpath='/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/ionspector_wtdbgp17/'
	ref='/data/scratch/maggic/simulation/denovosimulation/halploid_chr1/chr1.fa'
	check_depth_ref(outpath,ref)
