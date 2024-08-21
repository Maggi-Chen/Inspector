import pysam
import time
import os

def cigardeletion_ref(flag,chrom,position,cigar,min_size,max_size): 
	pos=int(position)
	numbers='1234567890'
	num=''
	reflen=0
	readlen=0
	leftclip=0
	rightclip=0
	deletions=[]
	insertions=[]
	for c in cigar:
		if c in numbers:
			num+=c; continue
		if c in 'MNP=X':
			readlen+=int(num); reflen+=int(num);  num='';  continue
		if c=='I':
			if  int(num)>=min_size and int(num)<=max_size:
				insertions+=[[chrom,pos+reflen,int(num),'I-cigar',readlen+leftclip]]
			readlen+=int(num)
			num=''; continue
		if c == 'D':
			if  int(num)>=min_size and int(num)<=max_size:
				deletions+=[[chrom,pos+reflen,int(num),'D-cigar',readlen+leftclip]]
			reflen+=int(num);  num='';  continue

		if c in 'SH':
			if readlen==0:
				leftclip=int(num)
			else:
				rightclip=int(num)
			num=''; continue
	testif=1
	window=500
	while testif==1:
		testif=0
		if len(deletions)==1:
			break
		i=len(deletions)-1
		while i>0:
			gaplength=deletions[i][1]-deletions[i-1][1]-deletions[i-1][2]
			if gaplength <= window:
				deletions=deletions[:i-1]+[[chrom,deletions[i-1][1],deletions[i-1][2]+deletions[i][2],'D-cigar',deletions[i-1][4]]]+deletions[i+1:]
				testif=1;break
			else:
				i-=1
	testif=1
	while testif==1:
		testif=0
		if len(insertions)==1:
			break
		i=len(insertions)-1
		while i>0:
			gaplength=insertions[i][1]-insertions[i-1][1]
			l1=insertions[i][2]
			l2=insertions[i-1][2]
			window=200 if max(l1,l2)<100 else 400
			window=400 if window==400 and max(l1,l2) <500 else 600
			if gaplength >window :
				i-=1
			else:
				insertions=insertions[:i-1]+[[chrom,insertions[i-1][1],l1+l2,'I-cigar',insertions[i-1][4]]]+insertions[i+1:]
				testif=1;break

	svcallset=deletions+insertions
	return [svcallset,reflen,[leftclip,readlen,rightclip]]


def cigardeletion(flag,chrom,position,cigar,min_size,max_size):	#input a read line, return list of deletions
	flag=int(flag)
	if flag<=16:
		detect_cigar_sv=True
	else:
		detect_cigar_sv=True
	pos=int(position)
	numbers='1234567890'
	num=''
	reflen=0
	readlen=0
	leftclip=0
	rightclip=0
	deletions=[]
	insertions=[]
	for c in cigar:
		if c in numbers:
			num+=c
			continue
		if c in 'MNP=X':
			readlen+=int(num); reflen+=int(num);  num='';  continue
		if c=='I':
			if detect_cigar_sv and int(num)>=min_size and int(num)<=max_size:
				insertions+=[[chrom,pos+reflen,int(num),'I-cigar',readlen]]
			readlen+=int(num)
			num=''; continue
		if c == 'D':
			if detect_cigar_sv and  int(num)>=min_size and int(num)<=max_size:
				deletions+=[[chrom,pos+reflen,int(num),'D-cigar']]
			reflen+=int(num);  num='';  continue
		if c in 'SH':
			if readlen==0:
				leftclip=int(num)
			else:
				rightclip=int(num)
			num=''; continue
	#merge deletions
	if detect_cigar_sv:
		testif=1
		window=500
		while testif==1:
			testif=0
			if len(deletions)==1:
				break
			i=len(deletions)-1
			while i>0:
				gaplength=deletions[i][1]-deletions[i-1][1]-deletions[i-1][2]
				if gaplength <= window:
					deletions=deletions[:i-1]+[[chrom,deletions[i-1][1],deletions[i-1][2]+deletions[i][2],'D-cigar']]+deletions[i+1:]
					testif=1
					break
				else:
					i-=1
		#merge insertions
		testif=1
		while testif==1:
			testif=0
			if len(insertions)==1:
				break
			i=len(insertions)-1
			while i>0:
				l1=insertions[i][2]
				l2=insertions[i-1][2]
				gaplength=insertions[i][1]-insertions[i-1][1]
				window=200 if max(l1,l2)<100 else 400
				window=400 if window==400 and max(l1,l2) <500 else 600
				if gaplength >window :
					i-=1
				else:
					insertions=insertions[:i-1]+[[chrom,insertions[i-1][1],l1+l2,'I-cigar',insertions[i-1][4]]]+insertions[i+1:]
					testif=1
					break
	
	svcallset=deletions+insertions
	return [svcallset,reflen,[leftclip,readlen,rightclip]]


def segmentdeletion_ref(segments,min_size,max_size,if_contig):
	segments=[c for c in segments if c[5][1]>=min(0.05*(c[5][0]+c[5][1]+c[5][2]),20000)]
	if len(segments)<=1:
		return []
	svcallset=[]
	for iii  in range(len(segments)-1):
		primary=segments[iii]
		chrom=primary[2]
		priflag=(int(primary[1])%32)>15
		samedirchr=[];samechr=[];diffchr=[]
		restsegments=segments[iii+1:]
		for c in restsegments:
			ch=c[2]; f=int(c[1])%32>15
			if c[5][1]<300:
				continue
			if ch!=chrom:
				diffchr+=[c]
			elif f!=priflag:
				samechr+=[c]
			else:
				samedirchr+=[c]
		for c in samedirchr:
			if c[3]>primary[3] and c[4]-primary[4]>-200:
				leftread=primary; rightread=c
			elif c[3]<primary[3] and primary[4]-c[4]>-200:
				leftread=c; rightread=primary
			else:
				continue
			leftinfo=leftread[5]
			rightinfo=rightread[5]
			window=300
			if if_contig:
				window=min(2000,leftinfo[1]/2,rightinfo[1]/2)
			if abs(rightread[3]-leftread[4])<=window:
				overlap=rightread[3]-leftread[4]
				ins_size=rightinfo[0]-leftinfo[1]-leftinfo[0]-overlap
				if min_size<=ins_size<=max_size:
					if not priflag:
						svcallset+=[chrom+'\t'+str(min(rightread[3],leftread[4]))+'\t'+str(ins_size)+'\t'+'I-segment'+'\t'+primary[0]+'\t'+str(int(c[1])+int(primary[1]))+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[0]+leftinfo[1])]
					else:
						svcallset+=[chrom+'\t'+str(min(rightread[3],leftread[4]))+'\t'+str(ins_size)+'\t'+'I-segment'+'\t'+primary[0]+'\t'+str(int(c[1])+int(primary[1]))+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[2])]

			overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
			window_max=2000 #Test for rescue FN
			overlap_window=-200
			if if_contig:
				overlap_window=-5000
				window_max=5000
			if overlap_window<overlapmap<window_max:
				del_size=rightread[3]-leftread[4]+overlapmap
				if min_size<=del_size<=max_size:
					if not priflag:
						svcallset+=[chrom+'\t'+str(leftread[4]-max(0,overlapmap))+'\t'+str(del_size)+'\t'+'D-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[0]+leftinfo[1])]
					else:
						svcallset+=[chrom+'\t'+str(leftread[4]-max(0,overlapmap))+'\t'+str(del_size)+'\t'+'D-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[2])]

		for c in samechr:
			if c[3]>primary[3] and c[4]-primary[4]>-200:
				leftread=primary; rightread=c
			elif c[3]<primary[3] and primary[4]-c[4]>-200:
				leftread=c; rightread=primary
			else:
				continue
			leftinfo=leftread[5]
			rightinfo=rightread[5]
			window_max=500
			overlap_window=-200
			if if_contig:
				overlap_window=-2000
			overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[2]
			if overlap_window<overlapmap<window_max and (rightread[4]-leftread[4])>=max(100,overlapmap):
				inv_size=rightread[4]-leftread[4]-overlapmap
				if min_size<=inv_size<=max_size:
					if int(leftread[1]) % 32 <16:
						svcallset+=[chrom+'\t'+str(leftread[4])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[0]+leftinfo[1])]
					else:
						svcallset+=[chrom+'\t'+str(leftread[4])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[2])]
					continue
			overlapmap=rightinfo[1]+rightinfo[2]-leftinfo[0]
			if overlap_window<overlapmap<window_max and (rightread[3]-leftread[3])>=max(100,overlapmap):
				inv_size=rightread[3]-leftread[3]-overlapmap
				if min_size<=inv_size<=max_size:
					if int(leftread[1]) % 32 <16:
						svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[0])]
					else:
						svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)+'\t'+str(leftinfo[1]+leftinfo[2])]
					continue
	return svcallset






def segmentdeletion(segments,min_size,max_size,if_contig):  #input a list of segments,return list of deletions
	if len([c for c in segments if int(c[1])<=16])==0:
		return []
	segments=[c for c in segments if c[5][1]>=min(0.05*(c[5][0]+c[5][1]+c[5][2]),20000)]
	if len(segments)<=1:
		return []
	svcallset=[]
	for iii  in range(len(segments)-1):
		primary=segments[iii]
		chrom=primary[2]
		priflag=(int(primary[1])%32)>15
		samedirchr=[]
		samechr=[]
		diffchr=[]
		restsegments=segments[iii+1:]
		for c in restsegments:
			ch=c[2]
			f=int(c[1])%32>15
			if c[5][1]<300:
				continue
			if ch!=chrom:
				diffchr+=[c]
			elif f!=priflag:
				samechr+=[c]
			else:
				samedirchr+=[c]
		for c in samedirchr:
			if c[3]>primary[3] and c[4]-primary[4]>-200:
				leftread=primary
				rightread=c
			elif c[3]<primary[3] and primary[4]-c[4]>-200:
				leftread=c
				rightread=primary
			else:
				continue	
			leftinfo=leftread[5]
			rightinfo=rightread[5]
			#insertion:
			window=300
			if if_contig:
				window=min(2000,leftinfo[1]/2,rightinfo[1]/2)
			if abs(rightread[3]-leftread[4])<=window:
				overlap=rightread[3]-leftread[4]
				ins_size=rightinfo[0]-leftinfo[1]-leftinfo[0]-overlap
				if min_size<=ins_size<=max_size:
					svcallset+=[chrom+'\t'+str(min(rightread[3],leftread[4]))+'\t'+str(ins_size)+'\t'+'I-segment'+'\t'+primary[0]+'\t'+str(int(c[1])+int(primary[1]))+'\t'+str((int(c[6])+int(primary[6]))/2)]

			#deletion:
			overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
			#window_max=1500
			window_max=2000 #Test for rescue FN
			overlap_window=-200
			if if_contig:
				overlap_window=-5000
				window_max=5000
			if overlap_window<overlapmap<window_max:
				del_size=rightread[3]-leftread[4]+overlapmap
				if min_size<=del_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[4]-max(0,overlapmap))+'\t'+str(del_size)+'\t'+'D-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
			'''
			#duplication:
			if not if_contig:
				overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
				window_max=500
				if -200<overlapmap<window_max and leftread[4]-rightread[3]>=max(50,overlapmap):
					dup_size=leftread[4]-rightread[3]-max(overlapmap,0)
					if min_size<=dup_size<=max_size:
						svcallset+=[chrom+'\t'+str(rightread[3])+'\t'+str(dup_size)+'\t'+'I-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
				overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[0]
				if -200<overlapmap<window_max and (rightread[4]-leftread[3])>=max(1000,overlapmap):
					dup_size=rightread[4]-leftread[3]-overlapmap
					if min_size<=dup_size<=max_size:
						svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(dup_size)+'\t'+'I-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
			'''
		#inversion:
		for c in samechr:
			if c[3]>primary[3] and c[4]-primary[4]>-200:
				leftread=primary
				rightread=c
			elif c[3]<primary[3] and primary[4]-c[4]>-200:
				leftread=c
				rightread=primary
			else:
				continue
			leftinfo=leftread[5]
			rightinfo=rightread[5]
			window_max=500
			overlap_window=-200
			if if_contig:
				overlap_window=-2000
			overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[2]
			if overlap_window<overlapmap<window_max and (rightread[4]-leftread[4])>=max(100,overlapmap):
				inv_size=rightread[4]-leftread[4]-overlapmap
				if min_size<=inv_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[4])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					continue
			overlapmap=rightinfo[1]+rightinfo[2]-leftinfo[0]
			if overlap_window<overlapmap<window_max and (rightread[3]-leftread[3])>=max(100,overlapmap):
				inv_size=rightread[3]-leftread[3]-overlapmap
				if min_size<=inv_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					continue
	return svcallset




def detect_sortbam(workpath,min_size,max_size,chrom):
	f=pysam.AlignmentFile(workpath+'read_to_contig.bam', "rb")
	segmentreads={}
	tempfile=open(workpath+'debreak_workspace/read_to_contig_'+chrom+'.debreak.temp','w')
	totalmaplength=0
	number_read=0
	split_num=0

	for align in f.fetch(chrom,):
		if align.is_secondary:
			continue
		readname=align.query_name
		flag=align.flag
		position=align.reference_start+1
		refend=align.reference_end+1
		cigar_info=[0,0,0]
		if align.cigar[0][0] in [4,5]:
			cigar_info[0]=align.cigar[0][1]
		if align.cigar[-1][0] in [4,5]:
			cigar_info[2]=align.cigar[-1][1]
		cigar_info[1]=align.query_alignment_length
		mappingquality=align.mapping_quality
		readinfo=[readname,flag,chrom,position,refend,cigar_info,mappingquality]
		if align.is_supplementary:
			cigar=align.cigarstring
			cigarinfo=cigardeletion(flag,chrom,position,cigar,5,max_size)
			cigarsv=[mm for mm in cigarinfo[0] if int(mm[2])>=min_size]
			for d in cigarsv:
				tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+str(flag)+'\t'+str(mappingquality)+'\n')

			pri_chrom=align.get_tag("SA").split(',')[0]
			if pri_chrom!=chrom:
				continue
			else:
				if readname not in segmentreads:
					segmentreads[readname]=[readinfo]
				else:
					segmentreads[readname]+=[readinfo]

		else:
			totalmaplength+=align.query_length
			number_read+=1
			cigar=align.cigarstring
			cigarinfo=cigardeletion(flag,chrom,position,cigar,5,max_size)
			cigarsv=[mm for mm in cigarinfo[0] if int(mm[2])>=min_size]
			for d in cigarsv:
				tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+str(flag)+'\t'+str(mappingquality)+'\n')

			if align.has_tag("SA") :
				if align.mapping_quality > 50:
					split_num+=1
				if chrom in [c.split(',')[0] for c in align.get_tag("SA").split(';')[:-1]]:
					if readname not in segmentreads:
						segmentreads[readname]=[readinfo]
					else:
						segmentreads[readname]+=[readinfo]

	for readgroup in segmentreads:
		if len(segmentreads[readgroup])<2 or len(segmentreads[readgroup])>20:
			continue
		segmentsv=segmentdeletion(segmentreads[readgroup],min_size,max_size,False)
		for d in segmentsv:
			tempfile.write(d+'\n')
	tempfile.close()
	f.close()
	if totalmaplength!=0:
		f=open(workpath+'map_depth/maplength_large_'+chrom,'w')
		f.write(str(totalmaplength)+'\n')
		f.close()
		f=open(workpath+'map_depth/readnum_large_'+chrom,'w')
		f.write(str(number_read)+'\n')
		f.close()
		f=open(workpath+'map_depth/splitread_large_'+chrom,'w')
		f.write(str(split_num)+'\n')
		f.close()
	return 0


def detect_sortbam_nosv(writepath,chrom,contig_type):
	print 'Collect info from '+chrom
	samfile=pysam.AlignmentFile(writepath+'read_to_contig.bam',"rb")
	allreads=samfile.fetch(chrom,)
	totalmaplength=0
	number_read=0
	split_num=0
	for align in allreads:
		if align.is_secondary or align.is_supplementary:
			continue
		totalmaplength+=align.query_length
		number_read+=1
		if align.has_tag("SA"):
			if align.mapping_quality > 50:
				split_num+=1

	if totalmaplength!=0:
		f=open(writepath+'map_depth/maplength_'+contig_type+'_'+chrom,'w')
		f.write(str(totalmaplength)+'\n')
		f.close()
		f=open(writepath+'map_depth/readnum_'+contig_type+'_'+chrom,'w')
		f.write(str(number_read)+'\n')
		f.close()
		f=open(writepath+'map_depth/splitread_'+contig_type+'_'+chrom,'w')
		f.write(str(split_num)+'\n')
		f.close()
	return 0



def detect_sam_ref(filename,readpath,writepath,min_size,max_size):
	f=open(readpath+filename,'r')
	c=f.readline()
	g=open(writepath+filename[:-4]+'.debreak.temp','w')
	lastname=''
	segments=['']
	unmapped=0
	mapped=0
	multimap=[]
	totalmappedlength=0
	while c!='':
		#remove headerlines, secondary alignments, alignment on scallfolds
		if c[0]=='@' or c.split('\t')[1] not in ['0','16','4','256','272','2048','2064']:
			c=f.readline(); continue
		if c.split('\t')[1]=='4':
			unmapped+=1; c=f.readline(); continue
		if c.split('\t')[1] in ['256','272'] :
			readname=c.split('\t')[0]
			if readname not in multimap:
				 multimap+=[readname]
			c=f.readline(); continue
		#detect the deletion from cigar
		readname=c.split('\t')[0]
		flag=c.split('\t')[1]
		chrom=c.split('\t')[2]
		position=int(c.split('\t')[3])
		mappingquality=c.split('\t')[4]
		cigar=c.split('\t')[5]
		cigarinfo=cigardeletion_ref('0',chrom,position,cigar,min_size,max_size)

		if flag in ['0','16']:
			mapped+=1
		totalmappedlength+=cigarinfo[2][1]


		if cigarinfo[2][1]<10000 or cigarinfo[2][1]<0.01*(sum(cigarinfo[2])):
		#if cigarinfo[2][1]<100000 or cigarinfo[2][1]<0.01*(sum(cigarinfo[2])):
		#if cigarinfo[2][1]<500000 and cigarinfo[2][1]<0.05*(sum(cigarinfo[2])):

			c=f.readline()
			continue

		refend=position+cigarinfo[1]
		cimplecigar=str(cigarinfo[2][0])+'\t'+str(cigarinfo[2][1])+'\t'+str(cigarinfo[2][2])
		# if primary: write deletions from cigar string
		cigarsv=cigarinfo[0]
		if int(flag)%32<16:
			for d in cigarsv:
				g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\t'+str(d[4])+'\n')
		else:
			totalreadlength=sum(cigarinfo[2])
			for d in cigarsv:
				if 'I-cigar' in d:
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\t'+str(totalreadlength-d[4]-d[2])+'\n')
				else:
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\t'+str(totalreadlength-d[4])+'\n')

		readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality]
		if readname!=lastname:
			if 1<len(segments):
				segmentd=segmentdeletion_ref(segments,min_size,max_size,True)
				for d in segmentd:
					g.write(d+'\n')
			lastname=readname
			segments=[readinfo]
		else:
			segments+=[readinfo]
		c=f.readline()
	if 1<len(segments):
		segmentd=segmentdeletion_ref(segments,min_size,max_size,True)
		for d in segmentd:
			g.write(d+'\n')
	segments=[]
	f.close()
	g.close()
	
	return [unmapped,mapped,len(multimap),totalmappedlength]

		

