import os
import time
import sys
import pysam

def mergeinfolengthsort(a):
	return int(a.split('\t')[2])

def mergeinfo_insertion(candi,min_support):
	candi.sort(key=mergeinfolengthsort)

	if len(candi)>=1.5*min_support:
		upper=int(candi[len(candi)*3//4].split('\t')[2])
		lower=int(candi[len(candi)//4].split('\t')[2])
		if upper>1.75*lower and upper-lower>50:
			svgroups=assign_candi_insertion(candi,upper,lower)
			svgroups=assign_candi_insertion(candi,svgroups[2],svgroups[3])
			svgroups=assign_candi_insertion(candi,svgroups[2],svgroups[3])
			mergedsv=[]
			if len(svgroups[0])>=min_support:
				mergedsv+=mergeinfo_insertion_oneevent(svgroups[0],min_support)
			if len(svgroups[1])>=min_support:
				mergedsv+=mergeinfo_insertion_oneevent(svgroups[1],min_support)
			if len(mergedsv)==2:
				mergedsv=[c+'\tCompoundSV' for c in mergedsv]
			if len(mergedsv)==1:
				mergedsv=[mergedsv[0]+'\tUnique']
			return mergedsv
	mergedsv=mergeinfo_insertion_oneevent(candi,min_support)
	if len(mergedsv)==1:
		return [mergedsv[0]+'\tUnique']
	else:
		return []

def assign_candi_insertion(candi,mean1,mean2):
	group1=[]
	group2=[]
	for c in candi:
		if abs(int(c.split('\t')[2])-mean1)<=abs(mean2-int(c.split('\t')[2])):
			group1+=[c]
		else:
			group2+=[c]
	mean1_new=sum([int(c.split('\t')[2]) for c in group1])//len(group1)
	mean2_new=sum([int(c.split('\t')[2]) for c in group2])//len(group2)
	return [group1,group2,mean1_new,mean2_new]


def mergeinfo_insertion_oneevent(candi,min_support):
	candi.sort(key=mergeinfolengthsort)
	min_support=max(2,min_support)
	while len(candi)>max(2,min_support-2):
		if int(candi[-1].split('\t')[2]) > 2* int(candi[len(candi)//2].split('\t')[2]) and  int(candi[-1].split('\t')[2]) -int(candi[len(candi)//2].split('\t')[2]) >30:
			candi.remove(candi[-1])
			continue
		if int(candi[len(candi)//2].split('\t')[2]) >  2*int(candi[0].split('\t')[2]) and int(candi[len(candi)//2].split('\t')[2]) -int(candi[0].split('\t')[2]) >30:
			candi.remove(candi[0])
			continue
		break
	if len(candi)>=max(2,min_support):
		chrom=candi[0].split('\t')[0]
		position=[int(c.split('\t')[1]) for c in candi]
		length=[int(c.split('\t')[2]) for c in candi]
		quality=[float(c.split('\t')[6]) for c in candi]
		position=sum(position)//len(position)
		quality=sum(quality)/float(len(quality))
		length=sum(length)//len(length)
		readnames=''
		for c in candi:
			readnames+=c.split('\t')[4]+';'
		readnames=readnames[:-1]
		numread=len(readnames.split(';'))
		return[chrom+'\t'+str(position)+'\t'+str(length)+'\t'+str(len(candi))+'\t'+str(numread)+'\t'+str(quality)+'\t'+'\t'+readnames]
	else:
		return []

def counttimesort(a):
	return [int(a.split('\t')[1]),int(a.split('\t')[2])]

def cluster(outpath,chrom,contiglength,mins,maxdepth):
	allsv=open(outpath+'debreak_workspace/read_to_contig_'+chrom+'.debreak.temp','r').read().split('\n')[:-1]

	# Large DEL
	largesv=[c for c in allsv if 'D-' in c and int(c.split('\t')[2])>2000]
	window=1600
	largesv.sort(key=counttimesort)
	largedel=[]
	start=0
	candi=[]
	for event in largesv:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=mins:
			largedel+=mergeinfo_insertion(candi,mins)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=mins:
		largedel+=mergeinfo_insertion(candi,mins)
		candi=[]

	#smaller DEL
	allsv=[c for c in allsv if 'D-' in c and int(c.split('\t')[2])<=3000]
	genomeposition=[0]*contiglength

	for c in allsv:
		start=int(c.split('\t')[1])
		end=int(c.split('\t')[1])+int(c.split('\t')[2])
		original=genomeposition[start-1:end-1]
		new=[mm+1 for mm in original]
		genomeposition[start-1:end-1]=new
	svregion=[]
	inblock=False
	threshold=3

	for i in range(len(genomeposition)):
		if inblock:
			if genomeposition[i]>=max(maxdep/10.0,threshold):
				localdep+=[genomeposition[i]]
				if genomeposition[i]>maxdep:
					maxdep=genomeposition[i]
			else:
				inblock=False
				end=i
				if maxdep<=maxdepth:
					peakpos=localdep.index(maxdep)
					peakleftsize=0
					for i in range(peakpos):
						if localdep[peakpos-i-1]>=maxdep/10.0:
							peakleftsize+=1
						else:
							break
					svregion+=[(start+peakpos-peakleftsize,end,maxdep)]
				start=0
				end=0
				maxdep=0

		else:
			if genomeposition[i] > threshold:
				inblock=True
				localdep=[genomeposition[i]]
				start=i
				maxdep=genomeposition[i]

	svregion=[c for c in svregion if c[2] < maxdepth]
	allsvinfo={}
	for c in svregion:
		allsvinfo[c]=[]

	for c in allsv:
		start=int(c.split('\t')[1])
		end=start+int(c.split('\t')[2])
		for d in svregion:
			if min(end,d[1])-max(d[0],start)>0:
				allsvinfo[d]+=[c]

	sv=[]
	for c in svregion:
		svinfo=allsvinfo[c]
		sv+=mergeinfo_insertion(svinfo,mins)

	newsv=[]
	for c in largedel:
		testif=0
		for d in sv:
			if min(int(c.split('\t')[1])+int(c.split('\t')[2]), int(d.split('\t')[1])+int(d.split('\t')[2])) - max(int(c.split('\t')[1]),int(d.split('\t')[1]))>0 and 0.8*int(d.split('\t')[2])<=int(c.split('\t')[2])<=int(d.split('\t')[2])/0.8:
				testif=1; break
		if testif==0:
			newsv+=[c]
	newsv+=sv
	newsv.sort(key=counttimesort)


	if newsv==[]:
		return 0

	f=open(outpath+'ae_merge_workspace/del_merged_'+chrom,'w')
	for c in newsv:
		f.write(c+'\n')
	f.close()

	return 0



def cluster_ins(outpath,chrom,contiglength,mins,maxdepth,svtype):
	allsv=open(outpath+'debreak_workspace/read_to_contig_'+chrom+'.debreak.temp','r').read().split('\n')[:-1]
	
	# Large INS
	if svtype=='ins':
		largesv=[c for c in allsv if 'I-' in c and int(c.split('\t')[2])>2000]
	else:
		largesv=[c for c in allsv if 'INV-' in c and int(c.split('\t')[2])>2000]

	window=1600
	largesv.sort(key=counttimesort)
	largedel=[]
	start=0
	candi=[]
	for event in largesv:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=mins:
			largedel+=mergeinfo_insertion(candi,mins)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=mins:
		largedel+=mergeinfo_insertion(candi,mins)
		candi=[]
	
	# Small INS
	if svtype=='ins':
		allsv=[c for c in allsv if 'I-' in c and int(c.split('\t')[2])<=3000]
	else:
		allsv=[c for c in allsv if 'INV-' in c and int(c.split('\t')[2])<=3000]

	genomeposition=[0]*contiglength

	for c in allsv:
		start=int(c.split('\t')[1])-100
		end=int(c.split('\t')[1])+100
		original=genomeposition[start-1:end-1]
		new=[mm+1 for mm in original]
		genomeposition[start-1:end-1]=new
	
	svregion=[]
	inblock=False
	threshold=3

	for i in range(len(genomeposition)):
		if inblock:
			if genomeposition[i]>=max(maxdep/10.0,threshold):
				localdep+=[genomeposition[i]]
				if genomeposition[i]>maxdep:
					maxdep=genomeposition[i]
			else:
				inblock=False
				end=i
				if maxdep<=maxdepth:
					peakpos=localdep.index(maxdep)
					peakleftsize=0
					for i in range(peakpos):
						if localdep[peakpos-i-1]>=maxdep/10.0:
							peakleftsize+=1
						else:
							break
					svregion+=[(start+peakpos-peakleftsize,end,maxdep)]
				start=0;end=0;maxdep=0

		else:
			if genomeposition[i] > threshold:
				inblock=True
				localdep=[genomeposition[i]]
				start=i
				maxdep=genomeposition[i]

	svregion=[c for c in svregion if c[2] < maxdepth]
	allsvinfo={}
	for c in svregion:
		allsvinfo[c]=[]

	for c in allsv:
		start=int(c.split('\t')[1])-50
		end=start+100
		for d in svregion:
			if min(end,d[1])-max(d[0],start)>0:
				allsvinfo[d]+=[c]
	sv=[]
	for c in svregion:
		svinfo=allsvinfo[c]
		mergedins=mergeinfo_insertion(svinfo,mins)
		for m in mergedins:
			sv+=[m+'\t'+chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])]

	newsv=[]
	for c in largedel:
		testif=0
		for d in sv:
			if min(int(c.split('\t')[1])+int(c.split('\t')[2]), int(d.split('\t')[1])+int(d.split('\t')[2])) - max(int(c.split('\t')[1]),int(d.split('\t')[1]))>0 and 0.8*int(d.split('\t')[2])<=int(c.split('\t')[2])<=int(d.split('\t')[2])/0.8:
				testif=1; break
		if testif==0:
			newsv+=[c]
	newsv+=sv
	newsv.sort(key=counttimesort)

	if newsv==[]:
		return 0

	if svtype=='ins':
		f=open(outpath+'ae_merge_workspace/ins_merged_'+chrom,'w')
	else:
		f=open(outpath+'ae_merge_workspace/inv_merged_'+chrom,'w')

	for c in newsv:
		f.write(c+'\n')
	f.close()
	return 0


def genotype(depth,outpath):
	highcov=depth*2
	allae=open(outpath+'assembly_errors.bed','r').read().split('\n')[:-1]
	samfile=pysam.AlignmentFile(outpath+'read_to_contig.bam',"rb")
	f=open(outpath+'assembly_errors.bed-gt','w')
	coll=0;expan=0;inv=0

	for c in allae:
		chrom=c.split('\t')[0]
		start=int(c.split('\t')[1])
		stop=int(c.split('\t')[2])
		if start<0:
			continue	
		if 'Expansion' in c or 'Inversion' in c:
			leftcov=samfile.count(chrom,max(start-100,0),start,read_callback='all')
			rightcov=samfile.count(chrom,stop,stop+100,read_callback='all')
		if 'Collapse' in c:
			leftcov=samfile.count(chrom,max(0,start-100),start,read_callback='all')
			rightcov=samfile.count(chrom,stop,stop+100,read_callback='all')
		gtinfo='./.'
		if int(c.split('\t')[3])>=0.6*min(leftcov,rightcov):
			gtinfo='1/1'
		else:
			gtinfo='1/0'
		f.write(c+'\t'+gtinfo+'\t'+str(leftcov)+'\t'+str(rightcov)+'\t'+str(min(rightcov,leftcov))+'\n')
	f.close()

def filterae(depth,outpath,min_size,datatype):
	allsv=open(outpath+'assembly_errors.bed-gt','r').read().split('\n')[:-1]
	if datatype=='hifi':
		rat=0.8
	else:
		rat=0.7

	highcov=depth*2
	lowcov=depth/2
	exp=[c for c in allsv if 'Exp' in c]
	col=[c for c in allsv if 'Col' in c]
	inv=[c for c in allsv if 'Inv' in c]
	new=[]
	exponly=[]
	for i in range(len(exp)):
		c=exp[i].split('\t')
		testif=0
		for d in col:
			if c[0]==d.split('\t')[0] and int(c[1])-250<=int(d.split('\t')[1])<=250+int(c[2]) and int(c[5].split('=')[1])<20*int(d.split('\t')[5].split('=')[1]):
				testif=1
				expread=c[6].split(';');goodexp=len(list(dict.fromkeys((expread))))
				colread=d.split('\t')[6].split(';'); goodcol=len(list(dict.fromkeys((colread))))
				totaln=len(list(dict.fromkeys((expread+colread))))
				if 0.33<=int(c[3])/float(d.split('\t')[3])<=3:
					if totaln<min(goodexp+goodcol//2,goodcol+goodexp//2):
						expsize=int(c[5].split('=')[1])
						colsize=int(d.split('\t')[5].split('=')[1])
						if expsize>colsize+min_size:
							new+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+str(goodexp)+'\t'+c[4]+'\tSize='+str(expsize-colsize)+'\t'+c[7]+'\t'+c[8]+'\t'+c[9]+'\t'+c[10]+'\t'+';'.join(expread)]
						if expsize<colsize-min_size:
							dd=d.split('\t')
							new+=[c[0]+'\t'+dd[1]+'\t'+dd[2]+'\t'+str(goodcol)+'\tCollapse\tSize='+str(colsize-expsize)+'\t'+dd[7]+'\t'+dd[8]+'\t'+dd[9]+'\t'+dd[10]+'\t'+';'.join(colread)]
						col.remove(d);break
					else:
						new+=[c[0]+'\t'+c[1]+';'+d.split('\t')[1]+'\t'+c[2]+';'+d.split('\t')[2]+'\t'+str(totaln)+'\tHaplotypeSwitch\tSize='+str(int(c[2])-int(c[1]))+';'+d.split('\t')[5].split('=')[1]+'\t-/-\t'+c[8]+'\t'+c[9]+'\t'+c[10]+'\t'+';'.join(expread)+':'+';'.join(colread)+'\t'+str(goodexp)+';'+str(goodcol)]
				if 0.33>int(c[3])/float(d.split('\t')[3]):
					dd=d.split('\t')
					new+=[c[0]+'\t'+dd[1]+'\t'+dd[2]+'\t'+str(goodcol)+'\tCollapse\tSize='+str(int(c[2])-int(c[1]))+';'+d.split('\t')[5].split('=')[1]+'\t-/-\t'+dd[8]+'\t'+dd[9]+'\t'+dd[10]+'\t'+';'.join(colread)]
				if int(c[3])/float(d.split('\t')[3])>3:
					new+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+str(goodexp)+'\tExpansion\tSize='+str(int(c[2])-int(c[1]))+';'+d.split('\t')[5].split('=')[1]+'\t-/-\t'+c[8]+'\t'+c[9]+'\t'+c[10]+'\t'+';'.join(expread)]
				col.remove(d);break
		if testif==0:
			exponly+=[exp[i]]
	allsv=new
	for c in exponly+col:
		c=c.split('\t')
		expread=c[6].split(';');goodexp=len(list(dict.fromkeys((expread))))
		allsv+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+str(goodexp)+'\t'+c[4]+'\t'+c[5]+'\t'+c[7]+'\t'+c[8]+'\t'+c[9]+'\t'+c[10]+'\t'+';'.join(expread)]
	
	for c in inv:
		c=c.split('\t')
		expread=c[5].split(';');goodexp=len(list(dict.fromkeys((expread))))
		allsv+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+str(goodexp)+'\t'+c[4]+'\tSize='+str(int(c[2])-int(c[1]))+'\t'+c[6]+'\t'+c[7]+'\t'+c[8]+'\t'+c[9]+'\t'+';'.join(expread)]
	new=[]
	for c in allsv:
		if  max([int(mm) for mm in c.split('\t')[5].split('=')[1].split(';')])<min_size:
			continue
		if int(c.split('\t')[3]) >=10 and int(c.split('\t')[3])>=rat*int(c.split('\t')[9]) and lowcov<=int(c.split('\t')[9])<highcov:
			new+=[c]; continue	
	f=open(outpath+'structural_error.bed','w')
	f.write('#Contig_Name\tStart_Position\tEnd_Position\tSupporting_Read\tType\tSize\tHaplotype_Info\tDepth_Left\tDepth_Right\tDepth_Min\tSupporting_Read_Name\tHaplotype_Switch_Info\n')
	for c in new:
		f.write(c+'\n')
	f.close()

	exp=len([mm for mm in new if 'Exp' in mm])
	col=len([mm for mm in new if 'Col' in mm])
	het=len([mm for mm in new if 'Haplo' in mm])
	inv=len([mm for mm in new if 'Inv' in mm])
	f=open(outpath+'summary_statistics','a')
	f.write('Structural error\t'+str(len(new))+'\nExpansion\t'+str(exp)+'\nCollapse\t'+str(col))
	f.write('\nHaplotype switch\t'+str(het)+'\nInversion\t'+str(inv)+'\n')
	f.close()

	os.system('rm '+outpath+'assembly_errors.bed*')
	os.system('rm '+outpath+'read_to_contig.debreak.temp')
	totalbase=0
	for c in new:
		if 'Inv' in c:
			totalbase+=200;continue
		size=min([int(mm) for mm in c.split('\t')[5].split('=')[1].split(';')]+[10000])
		totalbase+=size

	return totalbase



