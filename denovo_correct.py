import os
import pysam
import sys
import random
import time
import multiprocessing
import subprocess

def sort_snp(a):
	return int(a.split('\t')[1])

def get_snpcut_start(snp):
	if snp.split('\t')[7] in ['BaseSubstitution', 'SmallCollapse']:
		return int(snp.split('\t')[1])+1
	else:
		return int(snp.split('\t')[2])+1

def get_snpcut_end(snp):
	if 'BaseSubstitution' == snp.split('\t')[7]:
		return int(snp.split('\t')[1])
	else:
		return int(snp.split('\t')[1])+1

def base_correction(ctgseq,snpset,ctg):
	t1=time.time()
	snpset.sort(key=sort_snp)
	bad=[]
	for i in range(len(snpset)-1):
		if 'BaseSubstitution' in snpset[i+1] and 'SmallExpansion' in snpset[i] and min(int(snpset[i+1].split('\t')[2]),int(snpset[i].split('\t')[2])+1)-max(int(snpset[i+1].split('\t')[1]),int(snpset[i].split('\t')[1])+1)>0:
			bad+=[snpset[i],snpset[i+1]]
	snpset=[c for c in snpset if c not in bad]
	cutposinfo=[]
	if snpset==[]:
		return (ctgseq,snpset)
	for i in range(len(snpset)):
		cutinfo=[0,0,'']
		if i>0:
			cutinfo[0]=get_snpcut_start(snpset[i-1])
		cutinfo[1]=get_snpcut_end(snpset[i])
		if 'SmallExpansion' == snpset[i].split('\t')[7]:
			cutinfo[2]=''
		elif snpset[i].split('\t')[7]== 'BaseSubstitution' or snpset[i].split('\t')[7]== 'SmallCollapse':
			cutinfo[2]=snpset[i].split('\t')[4]
		else:
			print ('Warning: Possible error in small-error correction.')
		cutposinfo+=[cutinfo]
	newseq=''
	for cutinfo in cutposinfo:
		newseq+=ctgseq[cutinfo[0]:cutinfo[1]]+cutinfo[2]
	newseq+=ctgseq[get_snpcut_start(snpset[-1]):]
	t2=time.time()
	print ('Base error correction for ',ctg,' finished. Time cost: ',t2-t1)
	return (newseq,snpset)

def call_flye_timeout(datatype,outpath,aeinfo,outtime):
	testp = multiprocessing.dummy.Pool(1)
	testres = testp.apply_async(call_flye, args=(datatype,outpath,aeinfo))
	try:
		testout = testres.get(outtime)  # Wait timeout seconds for func to complete.
		return testout
	except multiprocessing.TimeoutError:
		print ('Flye assembly time out for ',aeinfo)
		raise

def call_flye(datatype,outpath,aeinfo):
	tt0=time.time()
	os.system('flye --'+datatype+' '+outpath+'assemble_workspace/read_ass_'+aeinfo+'.fa -o '+outpath+'assemble_workspace/flye_out_'+aeinfo+'/ -t 4  ')
	tt1=time.time()
	print ('FLYETIME for ',aeinfo,tt1-tt0)
	return 0


def findpos(aeset,snpset,bamfile,outpath,datatype,thread,outtime):
	snpsetshift=[c for c in snpset if 'Small' in c]
	snpsetshift.sort(key=sort_snp)
	new=[]
	bam=pysam.AlignmentFile(bamfile,'rb')
	ctg=aeset[0].split('\t')[0]
	aeinfolist={}
	for c in aeset:
		if 'Inversion' in c:
			continue
		if 'HaplotypeSwitch' in c:
			if int(c.split('\t')[11].split(';')[0])>=int(c.split('\t')[11].split(';')[1]):
				readgroup=c.split('\t')[10].split(':')[0].split(';')
				aestart=int(c.split('\t')[1].split(';')[0])
				aeend=int(c.split('\t')[2].split(';')[0])
				aesize=c.split('\t')[5].split('=')[1].split(';')[0]
				aeinfo=ctg+'__'+str(aestart)+'__'+str(aeend)+'__'+str(aesize)+'__exp'
				aeinfolist[c]=aeinfo
			else:
				readgroup=c.split('\t')[10].split(':')[1].split(';')
				aestart=int(c.split('\t')[1].split(';')[1])
				aeend=int(c.split('\t')[2].split(';')[1])
				aesize=c.split('\t')[5].split('=')[1].split(';')[1]
				aeinfo=ctg+'__'+str(aestart)+'__'+str(aeend)+'__'+str(aesize)+'__col'
				aeinfolist[c]=aeinfo
		else:
			readgroup=c.split('\t')[10].split(';')
			aestart=int(c.split('\t')[1])
			aeend=int(c.split('\t')[2])
			aesize=c.split('\t')[5].split('=')[1].split(';')[0]
			aeinfo=ctg+'__'+str(aestart)+'__'+str(aeend)+'__'+str(aesize)+'__exp' if 'Exp' in c else ctg+'__'+str(aestart)+'__'+str(aeend)+'__'+str(aesize)+'__col'
			aeinfolist[c]=aeinfo

		f=open(outpath+'assemble_workspace/read_ass_'+aeinfo+'.fa','w')
		allread=bam.fetch(ctg,max(0,aestart-2000),aeend+2000)
		iii=0
		for read in allread:
			if read.query_name not in readgroup or read.flag>16:
				continue
			f.write('>'+read.query_name+'\n'+read.query_sequence+'\n')
			iii+=1
		f.close()
	
	flyerun=multiprocessing.Pool(thread)	
	for c in aeinfolist:
		aeinfo=aeinfolist[c]
		flyerun.apply_async(call_flye_timeout,args=(datatype,outpath,aeinfo,outtime))
	flyerun.close()
	flyerun.join()

	for c in aeinfolist:
		aeinfo=aeinfolist[c]
		try:
			allctg=open(outpath+'assemble_workspace/flye_out_'+aeinfo+'/assembly.fasta','r').read().split('>')[1:]
		except:
			allctg=[]	
			print ('Inspector Assembly Fail ' ,aeinfo)
			os.system('rm -rf '+outpath+'assemble_workspace/flye_out_'+aeinfo+'/')
			continue
		if len(allctg)==1:
			f=open(outpath+'assemble_workspace/new_contig_'+ctg+'.fa','a')
			newassseq=''.join(allctg[0].split('\n')[1:-1])
			f.write('>'+aeinfo+'__newctg\n'+newassseq+'\n')
			f.close()
		else:
			print ('Inspector Multi/No Alignment ' ,aeinfo)
		os.system('rm -rf '+outpath+'assemble_workspace/flye_out_'+aeinfo+'/')
		
		shiftpos=0
		for d in snpsetshift:
			if int(d.split('\t')[2])<=aestart:
				if 'SmallCollapse' in d:
					shiftpos+=len(d.split('\t')[4])
				else:
					shiftpos-=len(d.split('\t')[3])
			else:
				break
		c=c.split('\t')
		new+=[ctg+'\t'+str(aestart+shiftpos)+'\t'+str(aeend+shiftpos)+'\t'+c[4]+'\t'+aesize+'\t'+aeinfo]
	return new


def substitute_seq(ctgseq,newseq,ctgstart,ctgend,newstart,newend,diffsize):
	newpart=newseq[newstart:newend]
	oldpart=ctgseq[ctgstart-1000-10:ctgend+1000+10]
	realdiff=len(newpart)-len(oldpart)
	if (realdiff/float(diffsize)>2 or realdiff/float(diffsize)<0.5) and (abs(diffsize-realdiff)>300 and realdiff*diffsize>0):
		return (ctgseq,False)
	leftside= check_same(newpart[:100],oldpart[10:110])
	rightside= check_same(newpart[-100:],oldpart[-110:-10])
	shift1=0
	if leftside<90:
		for shift1 in range(-5,5):
			leftside= check_same(newpart[:100],oldpart[10+shift1:110+shift1])
			if leftside>=90:
				break
	shift2=0
	if rightside <90:
		for shift2 in range(-5,5):
			rightside= check_same(newpart[-100:],oldpart[-110+shift2:-10+shift2])
			if rightside>=90:
				break
	if leftside>=90 and rightside>=90:
		ctgseq=ctgseq[:ctgstart-1000+shift1]+newpart+ctgseq[ctgend+1000+shift2:]
		return (ctgseq,True)
	else:
		ctgseq=ctgseq[:ctgstart-1000]+newpart+ctgseq[ctgend+1000:]
		return (ctgseq,True)

def ae_correct_within(seq,read,start,end,size):
	mapping=read.get_aligned_pairs()
	readstart=0;readend=0
	for c in mapping:
		if c[1]==start-1000:
			readstart=c[0]
		if c[1]==end+1000:
			readend=c[0]
			break
	if readstart!=0 and readend!=0:
		(seq,ifcorr)=substitute_seq(seq,read.query_sequence,start,end,readstart,readend,size)
		return (seq,ifcorr)
	else:
		return (seq,False)


def ae_correct_between(seq,align,start,end,size):
	readstart=0;readend=0
	for c in align:
		if c.reference_start < start-1000 and c.reference_end > start-100 :
			mapping=c.get_aligned_pairs()
			for m in mapping:
				if m[1]==start-1000:
					readstart=m[0]
		if c.reference_start < end+100 and c.reference_end  >  end+1000:
			mapping=c.get_aligned_pairs()
			for m in mapping:
				if m[1]==end+1000:
					readend=m[0]
	if readstart!=0 and readend!=0:
		(seq,ifcorr)=substitute_seq(seq,align[0].query_sequence,start,end,readstart,readend,size)
		return (seq,ifcorr)
	else:
		return (seq,False)


def ae_correct_expcol(seq,align,aetype):
	aeinfo=align[0].query_name
	start=int(aeinfo.split('__')[1])
	end=int(aeinfo.split('__')[2])
	if aetype=='exp':
		size=0-int(aeinfo.split('__')[3])
	else:
		size=int(aeinfo.split('__')[3])

	for read in align:
		if read.reference_start < start-1000 and end+1000<read.reference_end:
			(seq,ifcorr)=ae_correct_within(seq,read,start,end,size)
			return (seq,ifcorr)
	if len(align)<2:
		return (seq,False)
	(seq,ifcorr)=ae_correct_between(seq,align,start,end,size)
	return (seq,ifcorr)

def check_same(a,b):
	a=list(a)
	b=list(b)
	numsame=0
	for i in range(len(a)):
		if a[i]==b[i]:
			numsame+=1
	return numsame

def sortctg(a):
	return int(a.split('__')[1])


def ae_correction(ctgseq,aeset,outpath):	
	ctg=aeset[0].split('\t')[0]
		
	f=open(outpath+'assemble_workspace/old_contig_'+ctg+'.fa','w')
	f.write('>old_ctg_'+ctg+'\n'+ctgseq+'\n')
	f.close()
	try:
		allctg=open(outpath+'assemble_workspace/new_contig_'+ctg+'.fa','r').read().split('>')[1:]
	except:
		return (ctgseq,0)
	newcontig={}
	ctgname=[]
	for c in allctg:
		newcontig[c.split('\n')[0]]=c.split('\n')[1]
		ctgname+=[c.split('\n')[0]]
	ctgname.sort(key=sortctg,reverse=True)
	f=open(outpath+'assemble_workspace/new_contig_'+ctg+'.fa','w')
	for c in ctgname:
		f.write('>'+c+'\n'+newcontig[c]+'\n')
	f.close()
	
	os.system('minimap2 -a '+outpath+'assemble_workspace/old_contig_'+ctg+'.fa '+outpath+'assemble_workspace/new_contig_'+ctg+'.fa --MD --eqx -t 6 --secondary=no  -Y  > '+outpath+'assemble_workspace/ctgalignment_'+ctg+'.sam')
	
	
	alignfile=pysam.AlignmentFile(outpath+'assemble_workspace/ctgalignment_'+ctg+'.sam','r')
	aeset.sort(key=sort_snp,reverse=True)

	allalign=alignfile.fetch(until_eof=True)
	lastreadname=''
	samectg=[]
	numcorr=0
	for aligninfo in allalign:
		if aligninfo.flag==4:
			print (aligninfo.query_name,' contig not aligned.')
			continue
		if aligninfo.query_name == lastreadname:
			samectg+=[aligninfo]
			continue
		if samectg!=[]:
			if 'exp' in lastreadname:
				(ctgseq,ifcorr)=ae_correct_expcol(ctgseq,samectg,'exp')
				if ifcorr:
					numcorr+=1
			if 'col' in lastreadname:
				(ctgseq,ifcorr)=ae_correct_expcol(ctgseq,samectg,'col')
				if ifcorr:
					numcorr+=1
		lastreadname=aligninfo.query_name
		samectg=[aligninfo]
	if samectg!=[]:
		if 'exp' in lastreadname:
			(ctgseq,ifcorr)=ae_correct_expcol(ctgseq,samectg,'exp')
			if ifcorr:
				numcorr+=1
		if 'col' in lastreadname:
			(ctgseq,ifcorr)=ae_correct_expcol(ctgseq,samectg,'col')
			if ifcorr:
				numcorr+=1
	logf=open(outpath+'Inspector_correct.log','a')
	logf.write('total ae'+str(len(aeset))+', corrected error '+str(numcorr)+'\n')
	logf.close()
	return (ctgseq,numcorr)

	mapinfo={}
	for c in aeset:
		mapinfo[c.split('\t')[5]]=int(c.split('\t')[1])
	allread=alignfile.fetch(until_eof=True)
	for aligninfo in allread:
		if aligninfo.is_secondary:
			continue
		if type( mapinfo[aligninfo.query_name[:-8]])==int  and aligninfo.reference_start+1000 <mapinfo[aligninfo.query_name[:-8]] < aligninfo.reference_end-1000 :
			mapinfo[aligninfo.query_name[:-8]]=aligninfo

	correctedstructural=0
	for c in aeset:
		aeinfo=c.split('\t')[5]
		aestart=int(c.split('\t')[1])
		aeend=int(c.split('\t')[2])
		aligninfo=mapinfo[aeinfo]
		if type(aligninfo)==int:
			continue
		cigar=aligninfo.cigarstring
		refpos=aligninfo.reference_start
		readpos=0
		num=''
		readstart=-1
		readend=-1
		for m in cigar:
			if m in '1234567890':
				num+=m;continue
			if m in 'M=XD':
				refpos+=int(num)
				if m!='D':
					readpos+=int(num)
				if readstart==-1 and refpos>= aestart-1000:
					readstart=readpos-(refpos-aestart+1000)
				if readend==-1 and refpos>=aeend+1000:
					readend=readpos-(refpos-aeend-1000)
				if readstart>0 and readend>0: 
					break
				num='';continue
			if m in 'IS':
				readpos+=int(num);num='';continue
			if m=='H':
				num='';continue
		newseq=aligninfo.query_sequence[readstart:readend]
		oldseq=ctgseq[aestart-1000-10:aeend+1000+10]
		leftside= check_same(newseq[:100],oldseq[10:110])
		rightside= check_same(newseq[-100:],oldseq[-110:-10])
		shift1=0
		if leftside<90:
			for shift1 in range(-5,5):
				leftside= check_same(newseq[:100],oldseq[10+shift1:110+shift1])
				if leftside>=90:
					break
		shift2=0
		if rightside <90:
			for shift2 in range(-5,5):
				rightside= check_same(newseq[-100:],oldseq[-110+shift2:-10+shift2])
				if rightside>=90:
					break
		if leftside>=90 and rightside>=90:
			ctgseq=ctgseq[:aestart-1000+shift1]+newseq+ctgseq[aeend+1000+shift2:]
			correctedstructural+=1
	return (ctgseq,correctedstructural)



def error_correction_large(ctg,oldseq,aeset,snpset,bamfile,outpath,datatype,thread,flyeouttime):
	t0=time.time()
	(newseq,snpset)=base_correction(oldseq,snpset,ctg)
	if aeset!=[]:
		aeset=findpos(aeset,snpset,bamfile,outpath,datatype,thread,flyeouttime)
	if aeset!=[]:
		(newseq,numcorr)=ae_correction(newseq,aeset,outpath)
	ff=open(outpath+'contig_corrected_'+ctg+'.fa','w')
	ff.write('>'+ctg+'\n'+newseq+'\n')
	ff.close()
	t1=time.time()
	logf=open(outpath+'Inspector_correct.log','a')
	logf.write('TIME used for structural error correction of '+ctg+': '+str(t1-t0)+'\n')
	logf.close()
	return 0

def error_correction_small(ctg,oldseq,snpset,bamfile,outpath,datatype):
	t0=time.time()
	#snpset=[c for c in snpset if c.split('\t')[0]==ctg]
	(newseq,snpset)=base_correction(oldseq,snpset,ctg)
	ff=open(outpath+'contig_corrected_'+ctg+'.fa','w')
	ff.write('>'+ctg+'\n'+newseq+'\n')
	ff.close()
	t1=time.time()
	logf=open(outpath+'Inspector_correct.log','a')
	logf.write('TIME used for small error correction of '+ctg+':'+str(t1-t0)+'\n')
	logf.close()
	return 0



