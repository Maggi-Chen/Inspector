import matplotlib
matplotlib.use('Agg')
import pysam
import matplotlib.pyplot as plt

def plot_n100(outpath,minlen):
	ctglen=open(outpath+'contig_length_info','r').read().split('\n')[:-1]
	ctglen=[int(c.split('\t')[1]) for c in ctglen if int(c.split('\t')[1]) >= minlen]

	n100=[]
	x100=[]
	ctglen.sort(reverse=True)
	totallen=sum(ctglen)
	addlen=0
	lastlen=0
	for i in range(100):
		x100+=[i+1]

		while addlen < (i+1)/100.0*totallen:
			try:
				lastlen=ctglen.pop(1)
				addlen+=lastlen
			except:
				break
		n100+=[lastlen]
	plt.plot(x100,n100,linewidth=2)
	plt.xlabel('N1-N100')
	plt.ylabel('Contig Length /bp')
	plt.savefig(outpath+'plot_n1n100.pdf')
	print ('end n100')
	return 0


def plot_na100(outpath):
	samfile=pysam.AlignmentFile(outpath+'contig_to_ref.sam','r')
	allread=samfile.fetch()
	alignlen=[]
	for align in allread:
		if align.flag==4:
			continue
		alignlen+=[align.query_alignment_length]
	n100=[]
	x100=[]
	alignlen.sort(reverse=True)
	totallen=sum(alignlen)
	addlen=0
	lastlen=0
	for i in range(100):
		x100+=[i+1]
		while addlen < (i+1)/100.0*totallen:
			try:
				lastlen=alignlen.pop(1)
				addlen+=lastlen
			except:
				break
		n100+=[lastlen]
	plt.plot(x100,n100,linewidth=2)
	plt.xlabel('NA1-NA100')
	plt.ylabel('Contig Length /bp')
	plt.savefig(outpath+'plot_na1na100.pdf')
	print ('end na100')
	return 0



def findpos(c,ctglength,step,startrefpos,ctgstartpos):
	temppos=[]
	ctgname=c.query_name
	refpos=c.reference_start
	cigar=c.cigarstring
	if 'S' not in cigar.split('M')[0].split('=')[0] and 'H' not in cigar.split('M')[0].split('=')[0]:
		ctgpos=0
	else:
		ctgpos=int(cigar.split('M')[0].split('=')[0].split('S')[0].split('H')[0])
	currpos=ctgpos

	num=''

	for m in cigar:
		num+=m; continue
		if m in 'M=X':
			ctgpos+=int(num); refpos+=int(num); num=''
		if m=='I':
			ctgpos+=int(num);num=''
		if m =='D':
			refpos+=int(num); num='';continue
		if m in 'SH':
			if refpos>c.reference_start:
				break
			else:
				num=''
		while  ctgpos>=currpos+step:
			if ctgpos>=currpos+step*2:
				temppos+=[[currpos+step,refpos+startrefpos,ctgname]]
			else:
				temppos+=[[ctgpos,refpos+startrefpos,ctgname]]
			currpos+=step
	if c.flag in [16,2064]:
		temppos=[ [ctglength-mm[0],mm[1],mm[2]] for mm in temppos]
	updatestart=[]
	for mm in temppos :
		updatestart+=[[mm[0]+ctgstartpos,mm[1],mm[2]]]

	return updatestart



def plot_dotplot(outpath):
	print ('start dot plot')
	samfile=pysam.AlignmentFile(outpath+'contig_to_ref.bam','rb')
	allchrom=samfile.references
	allchromlen=samfile.lengths
	maxreflen=max(allchromlen)
	idex=allchromlen.index(maxreflen)
	maxchrom=allchrom[idex]
	print (maxchrom)
	allread=samfile.fetch(maxchrom)
	if maxreflen >= 10000000:
		step=10000
	elif maxreflen>=1000000:
		step=1000
	else:
		step=100

	ctgleninfo=open(outpath+'contig_length_info','r').read().split('\n')[:-1]
	ctglen={}
	for c in ctgleninfo:
		ctglen[c.split('\t')[0]]=int(c.split('\t')[1])
	alignedctg={}
	for align in allread:
		if align.query_name not in alignedctg:
			alignedctg[align.query_name]=align.query_alignment_length
		else:
			alignedctg[align.query_name]+=align.query_alignment_length
	longalignctg=[c for c in alignedctg if alignedctg[c] >= maxreflen/100.0]

	print (len(longalignctg))

	startpos=0
	contig_startpos={}
	for c in longalignctg:
		contig_startpos[c]=startpos
		startpos+=ctglen[c]
	allpos=[]

	allread=samfile.fetch(maxchrom)

	for align in allread:
		if align.query_name not in longalignctg:
			continue
		temppos=[]
		ctglength=ctglen[align.query_name]
		ctgstart=contig_startpos[align.query_name]
		temppos=findpos(align,ctglength,step,0,ctgstart)
		allpos+=temppos


	plotx=[c[0] for c in allpos]
	ploty=[c[1] for c in allpos]
	plotcolor=[c[2] for c in allpos]

	allcolor=set(plotcolor)
	colors={}
	i=1
	for c in allcolor:
		colors[c]=i
		i+=10
	plotcolor2=[colors[c] for c in plotcolor]
	size=[1]*len(plotcolor2)
	plt.scatter(plotx,ploty,s=size,c=plotcolor2)
	plt.xlabel('Reference Position')
	plt.ylabel('Contig Position')
	plt.savefig(outpath+'plot_synteny.pdf')
	return 0








