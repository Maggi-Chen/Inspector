#!/usr/bin/env python
import os
import argparse
import denovo_static
import debreak_detect_test as debreak_detect
import debreak_merge_clustering as debreak_cluster
import debreak_merge
#import denovo_baseerror
import multiprocessing
import math 
import time

t0=time.time()

parser=argparse.ArgumentParser(description='de novo assembly evaluator', usage='inspector.py [-h] -c contig.fa -r raw_reads.fa -o output_dict/')
parser.add_argument('--version', action='version', version='Inspector_v1.0.1')
parser.add_argument('-c','--contig',action='append', dest='contigfile',default=[],help='assembly contigs in .fa format',required=True)
#parser.add_argument('-r','--read',type=str,default=False,help='sequencing reads in .fa format',required=True)
parser.add_argument('-r','--read',type=str,default=False,help='sequencing reads in .fa format',required=True,nargs='+')
parser.add_argument('-d','--datatype',type=str,default='clr',help='Input read type. (clr, hifi, nanopore) [clr]')
parser.add_argument('-o','--outpath',type=str,default='./adenovo_evaluation-out/',help='output directory')
parser.add_argument('--ref',type=str,default=False,help='OPTIONAL reference genome in .fa format')

#parser.add_argument('--depth_plot',action='store_true',help='plot depth at all contigs.')
parser.add_argument('-t','--thread',type=int,default=8,help='number of threads. [8]')
parser.add_argument('--min_contig_length',type=int,default=10000,help='minimal length for a contig to be evaluated. [10000]')
parser.add_argument('--min_contig_length_assemblyerror',type=int,default=1000000,help='minimal contig length for assembly error detection. [1000000]')
parser.add_argument('--min_assembly_error_size',type=int,default=50,help='minimal size for assembly errors. [50]')
parser.add_argument('--max_assembly_error_size',type=int,default=4000000,help='maximal size for assembly errors. [4000000]')
parser.add_argument('--pvalue',type=float,default=False,help='p-value cut off for small-scale error identification. [0.01 for HiFi, 0.05 for others]')
parser.add_argument('--skip_read_mapping',action='store_true',default=False,help='skip the step of mapping reads to contig.')


parser.add_argument('--skip_structural_error',action='store_true',default=False,help='skip the step of identifying large structural errors.')
parser.add_argument('--skip_structural_error_detect',action='store_true',default=False,help='skip the step of detecting large structural errors.')
parser.add_argument('--skip_base_error',action='store_true',default=False,help='skip the step of identifying small-scale errors.')
parser.add_argument('--skip_base_error_detect',action='store_true',default=False,help='skip the step of detecting small-scale errors from pileup.')

denovo_args=parser.parse_args()

#print denovo_args.contig
print "Start Assembly evaluation with contigs: " + str(denovo_args.contigfile)
validate_read=[]
for inputfile in denovo_args.read:
	try:
		open(inputfile,'r')
		validate_read+=[inputfile]
	except:
		print 'Warning: cannot open input file \"'+inputfile+'\". Removed from list.'
if len(validate_read)==0:
	print 'Error: No valida input read file. Abort.'
	quit()

if denovo_args.datatype not in ['clr','hifi','nanopore']:
	print 'Warning: Invalid input datatype (--datatype/-d). Should be one of the following: clr, ccs, nanopore. Use clr as default.'
	denovo_args.datatype='clr'

# check input arguments
if denovo_args.outpath[-1]!='/':
	denovo_args.outpath+='/'

if not os.path.exists(denovo_args.outpath):
	os.mkdir(denovo_args.outpath)

if len(denovo_args.contigfile)==1:
	singlecontig=True
elif len(denovo_args.contigfile)==2:
	singlecontig=False
else:
	print 'Error: Input contig file should be either 1 fasta file or two halploid.fa files. Check input -c/--contig.'
	quit()

# Simple statistics of contigs	
contiginfo=denovo_static.simple(denovo_args.contigfile,denovo_args.outpath,denovo_args.min_contig_length,denovo_args.min_contig_length_assemblyerror)
chromosomes=contiginfo[0]
chromosomes_map=contiginfo[1]
chromosomes_large=contiginfo[2]
largecontig_length=contiginfo[7]
chromosomes_small=[mmm for mmm in chromosomes_map if mmm not in chromosomes_large]
totalcontiglen=contiginfo[3]
totalcontiglen_large=contiginfo[4]


t1=time.time()
print 'TIME: Before read mapping ',t1-t0

if not denovo_args.skip_read_mapping:
	inputfileid=1
	for inputfile in validate_read:
		print "minimap2 -a -Q -N 1 -I 10G -t " + str(denovo_args.thread) + "  "+denovo_args.outpath+"valid_contig.fa " + inputfile + " | samtools sort -@ " + str(denovo_args.thread) + " -o  "+denovo_args.outpath+"read_to_contig_"+str(inputfileid)+".bam"
		os.system("minimap2 -a -Q  -N 1 -I 10G -t " + str(denovo_args.thread) + "  "+denovo_args.outpath+"valid_contig.fa " + inputfile + " | samtools sort -@ " + str(denovo_args.thread) + " -o  "+denovo_args.outpath+"read_to_contig_"+str(inputfileid)+".bam")
		inputfileid+=1
	if len(validate_read)>1:
		os.system("samtools merge  "+denovo_args.outpath+"read_to_contig.bam  "+denovo_args.outpath+"read_to_contig_*.bam")
	else:
		os.system("mv "+denovo_args.outpath+"read_to_contig_1.bam "+denovo_args.outpath+"read_to_contig.bam")
	os.system("samtools index "+str(denovo_args.outpath)+"read_to_contig.bam")
	os.system("rm "+str(denovo_args.outpath)+"read_to_contig_*.bam")

t2=time.time()
print 'TIME: Read Alignment: ',t2-t1
# Structural assembly error detection


if not denovo_args.skip_structural_error_detect:
	os.system("mkdir "+denovo_args.outpath+"map_depth/")
	if not denovo_args.skip_structural_error:
		os.system("mkdir "+denovo_args.outpath+"debreak_workspace/")
		debreak_det=multiprocessing.Pool(denovo_args.thread)
		for i in range(len(chromosomes_large)):
			debreak_det.apply_async(debreak_detect.detect_sortbam,args=(denovo_args.outpath,denovo_args.min_assembly_error_size,denovo_args.max_assembly_error_size,chromosomes_large[i]))
		for i in range(len(chromosomes_small)):
			debreak_det.apply_async(debreak_detect.detect_sortbam_nosv,args=(denovo_args.outpath,chromosomes_small[i],'small'))
		debreak_det.close()
		debreak_det.join()
		os.system("cat "+denovo_args.outpath+"debreak_workspace/read_to_contig_*debreak.temp > "+denovo_args.outpath+"read_to_contig.debreak.temp")
	else:
		debreak_det=multiprocessing.Pool(denovo_args.thread)
		for i in range(len(chromosomes_large)):
			debreak_det.apply_async(debreak_detect.detect_sortbam_nosv,args=(denovo_args.outpath,chromosomes_large[i],'large'))
		for i in range(len(chromosomes_small)):
			debreak_det.apply_async(debreak_detect.detect_sortbam_nosv,args=(denovo_args.outpath,chromosomes_small[i],'small'))
		debreak_det.close()
		debreak_det.join()


cov=denovo_static.mapping_info_ctg(denovo_args.outpath,chromosomes_large,chromosomes_small,totalcontiglen,totalcontiglen_large)
minsupp=round(cov/10.0)

t3=time.time()
print 'TIME: Raw sv detection: ',t3-t2

aelen_structuralerror=0
if not denovo_args.skip_structural_error:
	os.system('mkdir '+denovo_args.outpath+'ae_merge_workspace')
	for chrom in largecontig_length:
		contiglength=largecontig_length[chrom]
		debreak_cluster.cluster(denovo_args.outpath,chrom,contiglength,minsupp,cov*2)
		debreak_cluster.cluster_ins(denovo_args.outpath,chrom,contiglength,minsupp,cov*2,'ins')
		debreak_cluster.cluster_ins(denovo_args.outpath,chrom,contiglength,minsupp,cov*2,'inv')
	denovo_static.assembly_info_cluster(denovo_args.outpath,denovo_args.min_assembly_error_size,denovo_args.max_assembly_error_size)
	debreak_cluster.genotype(cov,denovo_args.outpath)
	
	aelen_structuralerror=debreak_cluster.filterae(cov,denovo_args.outpath,denovo_args.min_assembly_error_size,denovo_args.datatype)

t4=time.time()
print 'TIME: SV clustering : ',t4-t3

# SNP & indel detection
aelen_baseerror=0
if not denovo_args.skip_base_error:
	import denovo_baseerror
	if not denovo_args.skip_base_error_detect:
		os.system('samtools faidx '+denovo_args.outpath+'valid_contig.fa')
		debreak_det=multiprocessing.Pool(denovo_args.thread)
		os.system('mkdir '+denovo_args.outpath+'base_error_workspace')
		for chrom in chromosomes_map:
			debreak_det.apply_async(denovo_baseerror.getsnv,args=(denovo_args.outpath,chrom,cov*2/5,cov*2))
		#denovo_baseerror.getsnv(denovo_args.outpath,chrom,cov*2/5,cov*2)
		debreak_det.close()
		debreak_det.join()

	aelen_baseerror=denovo_baseerror.count_baseerrror(denovo_args.outpath,totalcontiglen,denovo_args.datatype)

t5=time.time()
print 'TIME: SNV calling: ',t5-t4

#QV
if aelen_structuralerror+aelen_baseerror>0:
	print aelen_baseerror, aelen_structuralerror,totalcontiglen
	qv=-10 * math.log10(float(aelen_baseerror+aelen_structuralerror)/totalcontiglen)
	f=open(denovo_args.outpath+'summary_statistics','a')
	f.write('\nQV: '+str(qv)+'\n')
	f.close()

t6=time.time()
print 'TIME: QV calculation: ',t6-t5

# if reference is also provided:
if denovo_args.ref:
	mapinfo=os.system("minimap2 -a -I 10G --eqx -x asm5 -t " + str(denovo_args.thread/2) + " "+denovo_args.ref+" " + denovo_args.outpath + "valid_contig.fa  --secondary=no > "+ denovo_args.outpath+"contig_to_ref.sam")
	#print 'os.system return value: '+str(mapinfo)
	#quit()
	os.system("samtools sort -@ " + str(denovo_args.thread/2) + " "+ denovo_args.outpath+"contig_to_ref.sam -o  " + denovo_args.outpath+"contig_to_ref.bam")
	os.system("samtools index "+ denovo_args.outpath+"contig_to_ref.bam")
	chromosomes=denovo_static.get_ref_align_info(denovo_args.outpath,totalcontiglen)
	mapping_info=debreak_detect.detect_sam_ref("contig_to_ref.sam",denovo_args.outpath,denovo_args.outpath,denovo_args.min_assembly_error_size,denovo_args.max_assembly_error_size)
	
	minsupp=1

	allsvsignal=open(denovo_args.outpath+'contig_to_ref.debreak.temp','r').read().split('\n')[:-1]
	rawdelcalls={}; rawinscalls={};rawdupcalls={};rawinvcalls={}
	for chrom in chromosomes:
		rawdelcalls[chrom]=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[6]+'\t'+c.split('\t')[4] for c in allsvsignal if 'D-' in c and c.split('\t')[0]==chrom]
		rawinscalls[chrom]=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[6]+'\t'+c.split('\t')[4] for c in allsvsignal if 'I-' in c and c.split('\t')[0]==chrom]
		rawdupcalls[chrom]=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[6]+'\t'+c.split('\t')[4] for c in allsvsignal if 'DUP-' in c and c.split('\t')[0]==chrom]
		rawinvcalls[chrom]=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[6]+'\t'+c.split('\t')[4] for c in allsvsignal if 'INV-' in c and c.split('\t')[0]==chrom]
	for chrom in chromosomes:
		debreak_merge.merge_insertion(minsupp,0,denovo_args.outpath,rawinscalls[chrom],chrom,'ins',True,)
		debreak_merge.merge_deletion(minsupp,0,denovo_args.outpath,rawdelcalls[chrom],chrom,'del',True,)
		debreak_merge.merge_deletion(minsupp,0,denovo_args.outpath,rawdupcalls[chrom],chrom,'dup',True,)
		debreak_merge.merge_insertion(minsupp,0,denovo_args.outpath,rawinvcalls[chrom],chrom,'inv',True,)

	denovo_static.assembly_info_ref(denovo_args.outpath)
	
	denovo_static.basepair_error_ref(denovo_args.outpath,contiginfo[5])

t7=time.time()
print 'TIME: ref based mode: ',t7-t6

