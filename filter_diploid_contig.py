import sys
def test(size1,size2,svtype,svtype2):
	readpath='/data/scratch/maggic/denovo_evaluation/hg002_ccs/inspector_canu/'
	f=open(readpath+'contig_to_ref.debreak.temp.filtered','r')
	#f=open(readpath+'missed_assembly_error_merged','r')
	#f=open(readpath+'assembly_errors.bed-gt.refpos_filtered','r')
	new=f.read().split('\n')[:-1]
	f.close()
	vcf=[]
	for c in new:  
		if svtype2 in c :#and int(c.split('\t')[4].split('tg')[1])<=25: 
			vcf+=[c]

	#f=open(readpath+'test_map_to_large_contig/assembly_errors.bed-gt.refpos_filtered','r')
	#f=open('/data/scratch/maggic/simulation/denovosimulation/diploid_chr1_homohetero/makeref/sv_simulated','r')
	#f=open('/data/scratch/maggic/HG002/'+svtype,'r')
	#f=open('/data/scratch/maggic/debreak_result_11.15.19/HG002_pacbio_ccs/'+svtype+'-merged','r')
	f=open('/data/scratch/maggic/debreak_result_11.15.19/HG002_SVs_union/'+svtype+'-union','r')

	new=f.read().split('\n')[:-1]
	pac=[c for c in new if c.split('\t')[0]=='chr1']
	#pac=new

	f.close()
	#sensitivity / number of detected vcf
	groups=[]
	for c in pac:
		if groups==[]:
			groups=[[c]]
		else:
			i=0
			for group in groups:
				chrom=group[0].split('\t')[0]
				if c.split('\t')[0]==chrom:
					group+=[c]
					i=1
					break
			if i==0:
				groups+=[[c]]
	detected=[]
	detected2=[]
	no=[]

	correct=[]
	for c in vcf:
		iii=0
		chrom=c.split('\t')[0]
		s1=int(c.split('\t')[1])
		s2=int(c.split('\t')[2])+s1
		length=int(c.split('\t')[2])
		window=1000
		sizeratio=0.5
		for group in groups:
			ch=group[0].split('\t')[0]
			if chrom!=ch:
				continue
			for d in group:
				d1=int(d.split('\t')[1])
				d2=float(d.split('\t')[2])+d1
				if s1-window<=d1<=s1+window  and (sizeratio*float(d.split('\t')[2])<=length<=float(d.split('\t')[2])/sizeratio or abs(float(d.split('\t')[2])-length)<=200):
					detected+=[c]
					correct+=[d]

					
					#group.remove(d)
					if group==[]:
						groups.remove([])	
					
					break
	#precision / number of correct pac
	correct=[]
	groups=[]
	for c in vcf:
		if groups==[]:
			groups=[[c]]
		else:
			i=0
			for group in groups:
				chrom=group[0].split('\t')[0]
				if c.split('\t')[0]==chrom:
					group+=[c]
					i=1
					break
			if i==0:
				groups+=[[c]]
	for c in pac:
		iii=0
		chrom=c.split('\t')[0]
		s1=int(c.split('\t')[1])
		s2=int(c.split('\t')[2])+s1
		length=int(c.split('\t')[2])
		window=1000
		sizeratio=0.5
		for group in groups:
			ch=group[0].split('\t')[0]
			if chrom!=ch:
				continue
			for d in group:
				d1=int(d.split('\t')[1])
				d2=float(d.split('\t')[2])+d1
				if s1-window<=d1<=s1+window  and (sizeratio*float(d.split('\t')[2])<=length<=float(d.split('\t')[2])/sizeratio or abs(length-float(d.split('\t')[2]))<=200):
					correct+=[c]
					#group.remove(d)
					if group==[]:
						groups.remove([])
					break
	#done




	false=[c for c in pac if c not in correct]
	missed=[c for c in vcf if c not in detected]
	if len(vcf)==0:
		sen=0.0
	else:
		sen=round(float(len(detected))/len(vcf)*10000)/100.0
	if len(pac)==0:
		pre=0.0
	else:
		pre=round(float(len(correct))/len(pac)*10000)/100.0
	if sen==0 and pre==0:
		f1=0
	else:
		f1=2*sen*pre/(sen+pre)/100
	
	print(str(size1)+'-'+str(size2)+'\t'+str(len(detected))+'\t'+str(len(vcf))+'\t'+str(sen)+'\t'+str(len(correct))+'\t'+str(len(pac))+'\t'+str(pre))
	print f1
	#print len([c for c in detected if '1/0' in c])

	
	f=open(readpath+'contig_to_ref.debreak.temp.filtered_union','a')
	for c in missed:
		f.write(c+'\t'+svtype+'\n')
	f.close()
	
	return True


print 'Expansion:'
test(50,100000000000000000,'insertion','I-')
#test(50,100000000000000000,'insertion','Exp')

print 'Collapse:'
test(50,100000000000000000,'deletion','D-')
#test(50,100000000000000000,'deletion','Col')

