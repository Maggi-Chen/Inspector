import os
import subprocess
import pysam
import sys

def simple(contigfile,outpath,min_size,min_size_assemblyerror):
	if len(contigfile)==2:
		halp=True
	else:
		halp=False
	halpnum=1
	f=open(outpath+'valid_contig.fa','w')
	length=[]
	maxlen=0
	largecontiglength={}
	all_contigs=[]
	largecontigs=[]
	map_contigs=[]
	totallength=0
	totallength_large=0
	contig_length_info=[]
	for contig in contigfile:
		allcontig=open(contig,'r').read().split('>')[1:]
		for c in allcontig:
			c=c.split('\n')[:-1]
			contig_name=c[0].split(' ')[0]
			if halp:
				contig_name='HAP_'+str(halpnum)+'_'+contig_name
			all_contigs+=[contig_name]
			seq=''
			length1=0
			for cc in c[1:]:
				seq+=cc
				length1+=len(cc)
			length+=[length1]
			contig_length_info+=[contig_name+'\t'+str(length1)]
			if length1>maxlen:
				maxlen=length1
				maxcontig=contig_name
			if length1>=min_size:
				f.write('>'+contig_name+'\n'+seq+'\n')
				totallength+=length1
				map_contigs+=[contig_name]
			if length1>=min_size_assemblyerror:
				largecontigs+=[contig_name]
				largecontiglength[contig_name]=length1
				totallength_large+=length1

		halpnum+=1

	f.close()

	length.sort(reverse=True)
	f=open(outpath+'contig_length_info','w')
	for c in contig_length_info:
		f.write(c+'\n')
	f.close()
	f=open(outpath+'summary_statistics','w')
	f.write('Statics of contigs:\n')

	iii=0
	total=sum(length)/2
	for c in length:
		iii+=c
		if iii>=total:
			n50=c; break
	length_ae=[c for c in length if c > min_size_assemblyerror]


	f.write('Number of contigs\t'+str(len(length))+'\n')
	f.write('Number of contigs > '+str(min_size)+' bp\t'+str(len(map_contigs))+'\n')
	f.write('Number of contigs >'+str(min_size_assemblyerror)+' bp\t'+str(len(length_ae))+'\n')
	f.write('Total length\t'+str(sum(length))+'\n')
	f.write('Total length of contigs > '+str(min_size)+' bp\t'+str(totallength)+'\n')
	f.write('Total length of contigs >'+str(min_size_assemblyerror)+'bp\t'+str(sum(length_ae))+'\n')

	if len(length)==0:
		logf=open(outpath+'Inspector.log','a')
		logf.write('Error: No contigs found. Check if input file is empty or if --min_contig_length is too high.\n')
		f.close()
		logf.close()
		quit()

	if len(length_ae)==0:
		logf=open(outpath+'Inspector.log','a')
		logf.write('Warning: No contigs larger than '+str(min_size_assemblyerror)+'bp. No structural errors will be reported. Check if --min_contig_length_assemblyerror is too high.\n')
		logf.close()

	f.write('Longest contig\t'+str(length[0])+'\n')
	if len(length)>1:
		f.write('Second longest contig length\t'+str(length[1])+'\n')
	f.write('N50\t'+str(n50)+'\n')

	

	iii=0; total=sum(length)/2; n50=0
	for c in length:
		iii+=c
		if iii>total:
			n50=c; break
	f.write('N50 of contigs >1Mbp\t'+str(n50)+'\n\n\n')
	f.close()
	return [all_contigs,map_contigs,largecontigs,totallength,totallength_large,maxcontig,maxlen,largecontiglength]


def mapping_info_ctg(outpath,largechrom,smallchrom,contiglength,contiglength_large):

	f=open(outpath+'summary_statistics','a')
	f.write('Read to Contig alignment:\n')

	os.system('touch '+outpath+'map_depth/maplength_large_null')
	os.system('touch '+outpath+'map_depth/readnum_large_null')
	os.system('touch '+outpath+'map_depth/splitread_large_null')
	

	os.system('cat '+outpath+'map_depth/maplength_large_* > '+outpath+'map_depth/all_maplength_large')
	os.system('cat '+outpath+'map_depth/maplength_* > '+outpath+'map_depth/all_maplength_total')
	os.system('cat '+outpath+'map_depth/readnum_large_* > '+outpath+'map_depth/all_readnum_large')
	os.system('cat '+outpath+'map_depth/readnum_* > '+outpath+'map_depth/all_readnum_total')
	os.system('cat '+outpath+'map_depth/splitread_large_* > '+outpath+'map_depth/all_splitread_large')
	os.system('cat '+outpath+'map_depth/splitread_* > '+outpath+'map_depth/all_splitread_total')

	unmapped=int(pysam.AlignmentFile(outpath+'read_to_contig.bam','rb').unmapped)
	
	info=open(outpath+'map_depth/all_readnum_total','r').read().split('\n')[:-1]
	mapped=sum([int(ccc) for ccc in info])
	totalread=mapped+unmapped
	if totalread==0:
		logf=open(outpath+'Inspector.log','a')
		logf.write('Warning: No reads found in read_to_contig alignment.\n')
		logf.close()
		return 0
	mapprate=round(10000*float(mapped)/(totalread))/100.0
	f.write('Mapping rate /%\t'+str(mapprate)+'\n')

	info=open(outpath+'map_depth/all_splitread_total','r').read().split('\n')[:-1]
	splitread=sum([int(ccc) for ccc in info])
	splrate=round(10000*float(splitread)/mapped)/100.0
	f.write('Split-read rate /%\t'+str(splrate)+'\n')

	info=open(outpath+'map_depth/all_maplength_total','r').read().split('\n')[:-1]
	mappedlen=sum([int(ccc) for ccc in info])
	cov=round(10000*float(mappedlen)/contiglength)/10000.0
	f.write('Depth\t'+str(cov)+'\n')

	try:
		info=open(outpath+'map_depth/all_readnum_large','r').read().split('\n')[:-1]
		mapped=sum([int(ccc) for ccc in info])
		mapprate=round(10000*float(mapped)/(totalread))/100.0
		f.write('Mapping rate in large contigs /%\t'+str(mapprate)+'\n')

		info=open(outpath+'map_depth/all_splitread_large','r').read().split('\n')[:-1]
		splitread=sum([int(ccc) for ccc in info])
		splrate=round(10000*float(splitread)/mapped)/100.0
		f.write('Split-read rate in large contigs /%\t'+str(splrate)+'\n')

		info=open(outpath+'map_depth/all_maplength_large','r').read().split('\n')[:-1]
		mappedlen=sum([int(ccc) for ccc in info])
		cov=round(10000*float(mappedlen)/contiglength_large)/10000.0
		f.write('Depth in large conigs\t'+str(cov)+'\n\n\n')
		f.close()

	except:
		logf=open(outpath+'Inspector.log','a')
		logf.write('Warning: Failed to characterize read alignment in large contigs. \n')
		logf.close()

	return cov



def sort_sv(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]


def assembly_info_cluster(outpath,min_size,max_size):
	os.system("cat "+outpath+"ae_merge_workspace/del_merged_* > "+outpath+"ae_merge_workspace/deletion-merged")
	os.system("cat "+outpath+"ae_merge_workspace/ins_merged_* > "+outpath+"ae_merge_workspace/insertion-merged")
	os.system("cat "+outpath+"ae_merge_workspace/inv_merged_* > "+outpath+"ae_merge_workspace/inversion-merged")
	f=open(outpath+'assembly_errors.bed','w')
	alldel=open(outpath+'ae_merge_workspace/deletion-merged','r').read().split('\n')[:-1]
	alldel=[c for c in alldel if min_size<=int(c.split('\t')[2])<=max_size]
	for c in alldel:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tExpansion\tSize='+c[2]+'\t'+c[7]+'\n')
	allins=open(outpath+'ae_merge_workspace/insertion-merged','r').read().split('\n')[:-1]
	allins=[c for c in allins if min_size<=int(c.split('\t')[2])<=max_size]
	for c in allins:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+1)+'\t'+c[3]+'\tCollapse\tSize='+c[2]+'\t'+c[7]+'\n')
	allinv=open(outpath+'ae_merge_workspace/inversion-merged','r').read().split('\n')[:-1]
	allinv=[c for c in allinv if min_size<=int(c.split('\t')[2])<=max_size]
	for c in allinv:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tInversion\t'+c[7]+'\n')
	f.close()
	return 0


def assembly_info(outpath):
	
	os.system("cat "+outpath+"del-merged-* > "+outpath+"deletion-merged")
	os.system("cat "+outpath+"ins-merged-* > "+outpath+"insertion-merged")
	os.system("cat "+outpath+"dup-merged-* > "+outpath+"duplication-merged")
	os.system("cat "+outpath+"inv-merged-* > "+outpath+"inversion-merged")
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
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tExpansion\tSize='+c[2]+'\t'+c[6]+'\n')
		if 'Ins' in c :
			c=c.split('\t')
			f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+1)+'\t'+c[3]+'\tCollapse\tSize='+c[2]+'\t'+c[6]+'\n')
	allinv=open(outpath+'inversion-merged','r').read().split('\n')[:-1]
	for c in allinv:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+str(int(c[1])+int(c[2]))+'\t'+c[3]+'\tInversion\t'+c[6]+'\n')
	f.close()

	f=open(outpath+'summary_statistics','a')
	f.write('Number of assembly collapse\t'+str(len(allins))+'\n')
	f.write('Number of assembly expansion\t'+str(len(alldel))+'\n')
	f.write('Number of assembly inversion\t'+str(len(allinv))+'\n')
	f.close()
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
	f=open(outpath+'structural_errors_ref.bed','w')
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
	f=open(outpath+'small_scale_error_ref.bed','w')
	for c in svs:
		f.write(c+'\n')
	f.close()
	return 0

def sortblock(a):
	return [a[0],a[1]]

def count_ref_coverage(refcoveredall):
	allchrom=list(set([c[0] for c in refcoveredall]))
	new=[]
	for chrom in allchrom:
		refcovered=[c for c in refcoveredall if c[0]==chrom] 
		refcovered.sort(key=sortblock)
		ifovlp=0
		while len(refcovered)>1:
			if refcovered[0][2]<=refcovered[1][1]:
				new+=[refcovered[0]]
				refcovered=refcovered[1:]
			else:
				i=0
				ovlpstart=refcovered[i+1][1]; ovlpend=min(refcovered[i][2],refcovered[i+1][2])
				newblock=[]
				if refcovered[i+1][1] > refcovered[i][1]:
					newblock+=[[refcovered[i][0],refcovered[i][1],refcovered[i+1][1],refcovered[i][3]]]
				newblock+=[[refcovered[i][0],ovlpstart,ovlpend,refcovered[i][3]+refcovered[i+1][3]]]
				if refcovered[i+1][2] > refcovered[i][2]:
					newblock+=[[refcovered[i][0],refcovered[i][2],refcovered[i+1][2],refcovered[i+1][3]]]
				if refcovered[i+1][2]<refcovered[i][2]:
					newblock+=[[refcovered[i][0],refcovered[i+1][2],refcovered[i][2],refcovered[i][3]]]
				refcovered=newblock+refcovered[2:]
				refcovered.sort(key=sortblock)

		new+=refcovered

	b1=[c[2]-c[1] for c in new if c[3]==1]
	b2=[c[2]-c[1] for c in new if c[3]==2]
	b3=[c[2]-c[1] for c in new if c[3]>2]
	base1=sum(b1)
	base2=sum(b2)
	base3=sum(b3)
	return (base1,base2,base3)
	

def get_ref_align_info(path,totallength):
	f=pysam.AlignmentFile(path+'contig_to_ref.sam','r')
	allali=f.fetch()
	maplen=[]
	refcovered=[]
	for c in allali:
		if c.flag==4:
			continue
		readlen=c.query_alignment_length
		refcovered+=[[c.reference_name,c.reference_start,c.reference_end,1]]
		if c.flag in [0,2048]:
			leftclipinfo=c.cigartuples[0]
			leftclip = leftclipinfo[1] if leftclipinfo[0]==5 or leftclipinfo[0]==4 else 0
			leftclip = leftclipinfo[1] if leftclipinfo[0]==5 or leftclipinfo[0]==4 else 0
			maplen+=[[c.query_name,leftclip,leftclip+readlen]]
		if c.flag in [16,2064]:
			leftclipinfo=c.cigartuples[-1]
			leftclip = leftclipinfo[1] if leftclipinfo[0]==5 or leftclipinfo[0]==4 else 0
			maplen+=[[c.query_name,leftclip,leftclip+readlen]]
		
	n50info=[c[2]-c[1] for c in maplen]
	n50info.sort(reverse=True)
	lenacc=0
	na50=0
	info=sum(n50info)
	for c in n50info:
		lenacc+=c
		if lenacc>=0.5*info:
			na50=c;break

	assembly_maplenratio=float(info)/totallength
	(base1,base2,base3)=count_ref_coverage(refcovered)

	totalrefbase=sum(f.lengths)
	allrefchrom=list(f.references)

	base0=totalrefbase-base1-base2-base3
	f=open(path+'summary_statistics','a')
	f.write('\n\n\nReference-based mode:\n')
	f.write('Genome Coverage /% '+str(float(base1+base2+base3)/totalrefbase)+'\nReference base with Depth=0 (including Ns): '+str(base0)+';\t'+str(base0/float(totalrefbase)*100)+'%\n')
	f.write('Reference base with Depth=1 '+str(base1)+';\t'+str(base1/float(totalrefbase)*100)+'%\n')
	f.write('Reference base with Depth=2 '+str(base2)+';\t'+str(base2/float(totalrefbase)*100)+'%\n')
	f.write('Reference base with Depth>2 '+str(base3)+';\t'+str(base3/float(totalrefbase)*100)+'%\n')
	f.write('Assembly contig mapping ratio (length) /%\t'+str(assembly_maplenratio)+'\n')
	f.write('Assembly contig NA50\t'+str(na50)+'\n')
	f.close()

	return allrefchrom


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
	
	f=open(outpath+'summary_statistics','a')
	f.write('#BP with cov=0\t'+str(cov0*100.00/total)+'\n')
	f.write('#BP with cov=1\t'+str(cov1*100.00/total)+'\n')
	f.write('#BP with cov=2\t'+str(cov2*100.00/total)+'\n')
	f.write('#BP with cov>2\t'+str(cov3*100.00/total)+'\n')
	f.write('Coverage\t'+str(1-round(10000*float(cov0)/total)/10000.0)+'\n')
	f.close()
	
	return 0

