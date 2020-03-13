import pysam
import time
import os

def cigardeletion(flag,chrom,position,cigar,min_size,max_size):	#input a read line, return list of deletions
	flag=int(flag)
	if flag<=16:
		detect_cigar_sv=True
	else:
		detect_cigar_sv=False
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
		
			#duplication:
			if not if_contig:
				overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
				window_max=500
				if -200<overlapmap<window_max and leftread[4]-rightread[3]>=max(50,overlapmap):
					dup_size=leftread[4]-rightread[3]-max(overlapmap,0)
					if min_size<=dup_size<=max_size:
						svcallset+=[chrom+'\t'+str(rightread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
				overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[0]
				if -200<overlapmap<window_max and (rightread[4]-leftread[3])>=max(1000,overlapmap):
					dup_size=rightread[4]-leftread[3]-overlapmap
					if min_size<=dup_size<=max_size:
						svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
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
		'''
		#translocation:
		for c in diffchr:
			pinfo=primary[5]
			cinfo=c[5]
			window_max=200
			bp1=''
			bp2=''
			if abs(pinfo[0]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[0]-cinfo[1]-cinfo[2])<=window_max:
				chrom1=primary[2]
				bp1=primary[3]
			elif abs(pinfo[2]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[2]-cinfo[1]-cinfo[2])<=window_max:
				chrom1=primary[2]
				bp1=primary[4]
			if abs(cinfo[0]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[0]-pinfo[1]-pinfo[2])<=window_max:
				chrom2=c[2]
				bp2=c[3]
			elif abs(cinfo[2]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[2]-pinfo[1]-pinfo[2])<=window_max:
				chrom2=c[2]
				bp2=c[4]
			if bp1!='' and bp2!='':
				if 'chr' in chrom1 and 'chr' in chrom2:
					if 'X' in chrom1:
						svcallset+=[chrom2+'\t'+str(bp2)+'\t'+chrom1+'\t'+str(bp1)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					elif 'X' in chrom2:
						svcallset+=[chrom1+'\t'+str(bp1)+'\t'+chrom2+'\t'+str(bp2)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					elif int(chrom1.split('hr')[1]) > int(chrom2.split('hr')[1]):
						svcallset+=[chrom2+'\t'+str(bp2)+'\t'+chrom1+'\t'+str(bp1)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					elif int(chrom1.split('hr')[1])< int(chrom2.split('hr')[1]):
						svcallset+=[chrom1+'\t'+str(bp1)+'\t'+chrom2+'\t'+str(bp2)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
				else:
					if chrom1 > chrom2:
						svcallset+=[chrom2+'\t'+str(bp2)+'\t'+chrom1+'\t'+str(bp1)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
					if chrom1 < chrom2:
						svcallset+=[chrom1+'\t'+str(bp1)+'\t'+chrom2+'\t'+str(bp2)+'\tTRA-segment\t'+primary[0]+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))/2)]
		'''
	return svcallset




def detect_sam(filename,readpath,writepath,chromosomes,min_size,max_size):
	f=open(readpath+filename,'r')
	c=f.readline()
	g=open(writepath+filename[:-4]+'.debreak.temp','w')
	lastname=''
	segments=['']
	unmapped=0
	mapped=0
	multimap=[]
	totalmappedlength=0
	totalmappedlength_smallchr=0
	mapped_smallchr=0
	multimap_smallchr=[]
	while c!='':
		#remove headerlines, secondary alignments, alignment on scallfolds
		if c[0]=='@' or c.split('\t')[1] not in ['0','16','4','256','272','2048','2064']:
			c=f.readline(); continue
		if c.split('\t')[1]=='4':
			unmapped+=1
			c=f.readline(); continue
		if c.split('\t')[1] in ['256','272'] :
			readname=c.split('\t')[0]
			if c.split('\t')[2] in chromosomes:
				if readname not in multimap:
					multimap+=[readname]
			else:
				if readname not in multimap_smallchr:
					multimap_smallchr+=[readname]
			c=f.readline(); continue
		#detect the deletion from cigar 
		readname=c.split('\t')[0]
		flag=c.split('\t')[1]
		chrom=c.split('\t')[2]
		position=int(c.split('\t')[3])
		mappingquality=c.split('\t')[4]
		readseq=c.split('\t')[9]
		cigar=c.split('\t')[5]
		cigarinfo=cigardeletion(flag,chrom,position,cigar,min_size,max_size)

		if int(cigarinfo[2][1])<500:
			c=f.readline(); continue

		if chrom not in chromosomes:
			if flag in ['0','16']:
				mapped_smallchr+=1
			totalmappedlength_smallchr+=cigarinfo[2][1]
		else:
			if flag in ['0','16']:
				mapped+=1
			totalmappedlength+=cigarinfo[2][1]

		refend=position+cigarinfo[1]
		cimplecigar=str(cigarinfo[2][0])+'\t'+str(cigarinfo[2][1])+'\t'+str(cigarinfo[2][2])
		# if primary: write deletions from cigar string
		if int(c.split('\t')[1])%4096<2048:
			cigarsv=cigarinfo[0]
			for d in cigarsv:
				if 'I-cigar' in d:
					insertseq=readseq[d[4]:d[4]+d[2]]
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\t'+insertseq+'\n')
				else:
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\n')
		readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality]
		if readname!=lastname:
			if 1<len(segments)<=20:
				segmentd=segmentdeletion(segments,min_size,max_size,False)
				for d in segmentd:
					g.write(d+'\n')
			lastname=readname
			segments=[readinfo]
		else:
			segments+=[readinfo]



		c=f.readline()
	if 1<len(segments)<=20:
		segmentd=segmentdeletion(segments,min_size,max_size,False)
		for d in segmentd:
			g.write(d+'\n')
	segments=[]
	f.close()
	g.close()

	mapped_smallchr+=mapped
	multimap_smallchr+=multimap
	multimap_smallchr=list(dict.fromkeys(multimap_smallchr))
	totalmappedlength_smallchr+=totalmappedlength
	return [unmapped,mapped,len(multimap),totalmappedlength,mapped_smallchr,len(multimap_smallchr),totalmappedlength_smallchr]





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
		readseq=c.split('\t')[9]
		cigar=c.split('\t')[5]
		cigarinfo=cigardeletion('0',chrom,position,cigar,min_size,max_size)

		if flag in ['0','16']:
			mapped+=1
		totalmappedlength+=cigarinfo[2][1]

		if cigarinfo[2][1]<20000 or cigarinfo[2][1]<0.05*(sum(cigarinfo[2])):
			c=f.readline()
			continue

		refend=position+cigarinfo[1]
		cimplecigar=str(cigarinfo[2][0])+'\t'+str(cigarinfo[2][1])+'\t'+str(cigarinfo[2][2])
		# if primary: write deletions from cigar string
		cigarsv=cigarinfo[0]
		for d in cigarsv:
			if 'I-cigar' in d:
				insertseq=readseq[d[4]:d[4]+d[2]]
				g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\t'+insertseq+'\n')
			else:
				g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'\t'+flag+'\t'+mappingquality+'\n')
		readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality]
		if readname!=lastname:
			if 1<len(segments):
				segmentd=segmentdeletion(segments,min_size,max_size,True)
				for d in segmentd:
					g.write(d+'\n')
			lastname=readname
			segments=[readinfo]
		else:
			segments+=[readinfo]
		c=f.readline()
	if 1<len(segments):
		segmentd=segmentdeletion(segments,min_size,max_size,True)
		for d in segmentd:
			g.write(d+'\n')
	segments=[]
	f.close()
	g.close()
	
	return [unmapped,mapped,len(multimap),totalmappedlength]

		

if __name__ =="__main__":
	readpath='/data/scratch/maggic/denovo_evaluation/hg002_clr/inspector_canu/'
	writepath=readpath
	filename='contig_to_ref.sam'
	detect_sam_ref(filename,readpath,writepath,50,400000000)

	#filename='read_to_contig.sam'
	#chrom for haploid sim
	#chromosomes=['ctg1','ctg2','ctg3','ctg4','ctg5','ctg6','ctg7','ctg8','ctg9','ctg10','ctg11','ctg12','ctg13','ctg14','ctg15','ctg16','ctg17','ctg18','ctg19','ctg20','ctg21','ctg22','ctg23','ctg24','ctg25']
	#filename='read_to_large_contig.sam'
	#chromosomes=['tig00000035','tig00000109','tig00000141','tig00000148','tig00000154','tig00000190','tig00000224','tig00000264','tig00000359','tig00000419','tig00000446','tig00000506','tig00000589','tig00000676','tig00021687','tig00021690','tig00021694','tig00021695','tig00043385','tig00043387']
	#chromosomes=['scaffold_5','contig_9','contig_10','contig_12','contig_15','contig_17','contig_20','contig_22','contig_27','contig_28','contig_33','scaffold_36','contig_54','contig_55','contig_62','contig_69','contig_98','scaffold_111','contig_190','contig_225','scaffold_239']

	#chrom for diploid sim
	#chromosomes=['tig00022707','tig00022710','tig00022728','tig00022729','tig00022720','tig00000463','tig00022716','tig00000369','tig00022725','tig00000565','tig00022714','tig00000613','tig00022709','tig00022723','tig00000747','tig00000902','tig00000943','tig00000959','tig00001050','tig00022736','tig00022734','tig00001232','tig00022721','tig00001347']
	#chromosomes=['contig_10','contig_14','contig_18','contig_2','contig_20','contig_24','contig_26','contig_28','contig_35','contig_4','contig_44','contig_6','contig_8','contig_82','contig_86','contig_9','scaffold_21','scaffold_30','scaffold_69']
	#chromosomes=['ctg1','ctg2','ctg3','ctg4','ctg5','ctg6','ctg7','ctg8','ctg9','ctg10','ctg11','ctg12','ctg13','ctg14','ctg15','ctg16','ctg17','ctg18','ctg19','ctg20','ctg21','ctg22','ctg23','ctg24','ctg25','ctg26','ctg27']

	#detect_sam(filename,readpath,writepath,chromosomes,50,400000000)

