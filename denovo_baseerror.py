import re
import os
import statsmodels.stats.proportion
import statsmodels.stats.multitest

def find2(li):
	num=0
	val=''
	for c in  li:
		if li.count(c)>num:
			val=c
	return val

def getsnv(path,chrom,mincount,maxcov,mindepth):
	logf=open(path+'Inspector.log','a')
	logf.write('Start small-scale error detection for '+chrom+'\n')
	logf.close()
	g=open(path+'base_error_workspace/baseerror_'+chrom+'.bed','w')
	os.system('samtools mpileup -Q 0 '+path+'read_to_contig.bam -r '+chrom+' -o '+path+'base_error_workspace/base_'+chrom+'.pileup -f '+path+'valid_contig.fa')
	f=open(path+'base_error_workspace/base_'+chrom+'.pileup','r')
	a=f.readline()
	numbaseerror=0
	validctgbase=0
	if mindepth==False and type(mindepth)==bool:
		mindepth=maxcov/10.0

	while a!='':
		if a.split('\t')[2]!='N' and mindepth<=int(a.split('\t')[3]) <=maxcov:
			validctgbase+=1
		if int(a.split('\t')[3]) <mincount or int(a.split('\t')[3])-a.split('\t')[4].count('*') > maxcov:
			a=f.readline(); continue
		info=a.split('\t')[4]
		info=info.replace(',','.')
		info=re.sub('\^.','',info)
		info=info.replace('a','A')
		info=info.replace('t','T')
		info=info.replace('c','C')
		info=info.replace('g','G')
		depth=int(a.split('\t')[3])-info.count('*')
		min_supp=max(mincount,depth*0.2)
		ins=info.count('+')
		dels=info.count('-')
		ifindel=False
		if ins>=min_supp:
			ifindel=True;
			insinfp=info.split('+')[1:]
			insseq=[]
			for m in insinfp:
				num='';inum=0
				for dd in m:
					if dd in '1234567890':
						num+=dd; inum+=1
					else:
						break
				if int(num)<=mincount//2:
					insseq+=[m[inum:][:int(num)]]
				else:
					ins-=1
			if ins>=min_supp :
				mostf1=find2(insseq)
				numbaseerror+=1
				g.write(a.split('\t')[0]+'\t'+str(int(a.split('\t')[1])-1)+'\t'+a.split('\t')[1]+'\t-\t'+mostf1+'\t'+str(ins)+'\t'+str(depth)+'\tSmallCollapse\n')
		if dels>=min_supp:
			ifindel=True;
			insinfp=info.split('-')[1:]
			insseq=[]
			for m in insinfp:
				num='';inum=0
				for dd in m:
					if dd in '1234567890':
						num+=dd; inum+=1
					else:
						break
				if int(num)<=mincount//2:
					insseq+=[m[inum:][:int(num)]]
				else:
					dels-=1
			if dels>=min_supp:
				mostf1=find2(insseq)
				numbaseerror+=1
				g.write(a.split('\t')[0]+'\t'+str(int(a.split('\t')[1])-1)+'\t'+str(int(a.split('\t')[1])+len(mostf1)-1)+'\t'+mostf1+'\t-\t'+str(dels)+'\t'+str(depth)+'\tSmallExpansion\n')

		if info.count('.')+info.count('*')>0.8*int(a.split('\t')[3]) :
			a=f.readline(); continue
		acount=info.count('A')
		tcount=info.count('T')
		ccount=info.count('C')
		gcount=info.count('G')

		if '+'  in a or '-'  in info:
			insseq=''
			if '+' in info:
				insinfp=info.split('+')[1:]
				for m in insinfp:
					num=''
					inum=0
					for dd in m:
						if dd in '1234567890':
							num+=dd; inum+=1
						else:
							break
					insseq+=m[inum:][:int(num)]
			if '-' in info:
				insinfp=info.split('-')[1:]
				for m in insinfp:
					num=''
					inum=0
					for dd in m:
						if dd in '1234567890':
							num+=dd; inum+=1
						else:
							break
					insseq+=m[inum:][:int(num)]

			insacount=insseq.count('A')
			instcount=insseq.count('T')
			insccount=insseq.count('C')
			insgcount=insseq.count('G')

			acount-=insacount; tcount-=instcount; ccount-=insccount; gcount-=insgcount

		if max(acount,tcount,ccount,gcount) >=min_supp:
			if max(acount,tcount,ccount,gcount)==acount:
				altbase='A'
			if max(acount,tcount,ccount,gcount)==tcount:
				altbase='T'
			if max(acount,tcount,ccount,gcount)==ccount:
				altbase='C'
			if max(acount,tcount,ccount,gcount)==gcount:
				altbase='G'
			numbaseerror+=1
			g.write(a.split('\t')[0]+'\t'+str(int(a.split('\t')[1])-1)+'\t'+a.split('\t')[1]+'\t'+a.split('\t')[2]+'\t'+altbase+'\t'+str(max(acount,tcount,ccount,gcount))+'\t'+str(depth)+'\tBaseSubstitution\n')

		a=f.readline()
	f.close()
	os.system('rm '+path+'base_error_workspace/base_'+chrom+'.pileup')
	g.close()
	if numbaseerror==0:
		os.system('rm '+path+'base_error_workspace/baseerror_'+chrom+'.bed')
	f=open(path+'base_error_workspace/validbase','a')
	f.write(str(validctgbase)+'\n')
	f.close()
	return 0


def count_baseerrror(path,ctgtotallen,datatype,ave_depth):
	os.system('cat '+path+'base_error_workspace/baseerror_*bed > '+path+'base_error_workspace/allbaseerror.bed')
	allsnv=open(path+'base_error_workspace/allbaseerror.bed','r').read().split('\n')[:-1]
	snv=0;indelins=0;indeldel=0

	baseerror=[]
	iii=0
	if datatype=='hifi':
		propvalue=0.5
		pcutoff=0.01
		readcutoff=0.75
		if ave_depth<25:
			pcutoff=0.02
		if ave_depth<15:
			pcutoff=0.1
	else:
		propvalue=0.4
		pcutoff=0.05
		readcutoff=0.5
		if ave_depth<25:
			pcutoff=0.1


	allpvalue=[]

	for c in allsnv:
		p=0
		nread=int(c.split('\t')[5])
		depth=int(c.split('\t')[6])
		if nread<readcutoff*depth:
			continue
		for i in range(nread,depth+1):
			p+=statsmodels.stats.proportion.binom_test(i, depth, prop=propvalue, alternative='larger')

		
		if p<pcutoff :
			iii+=1
			baseerror+=[c+'\t'+str(p)]
			if 'BaseSubstitution' in c:
				snv+=1
			if 'SmallExpansion' in c:
				indeldel+=1
			if 'SmallCollapse' in c:
				indelins+=1
		
	per=float(iii)/ctgtotallen*1000000
	f=open(path+'small_scale_error.bed','w')
	f.write('#Contig_Name\tStart_Position\tEnd_Position\tBase_Contig\tBase_Read\tSupporting_Read\tDepth\tType\tPvalue\n')
	for c in baseerror:
		f.write(c+'\n')
	f.close()
	f=open(path+'summary_statistics','a')
	f.write('\n\nSmall-scale assembly error /per Mbp\t'+str(per)+'\nTotal small-scale assembly error\t'+str(iii)+'\nBase substitution\t'+str(snv)+'\nSmall-scale expansion\t'+str(indeldel)+'\n')
	f.write('Small-scale collapse\t'+str(indelins)+'\n')

	return iii	

	

