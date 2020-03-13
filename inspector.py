#!/usr/bin/env python
import os
import argparse
import denovo_static
import debreak_detect
import debreak_merge

parser=argparse.ArgumentParser(description='de novo assembly evaluator', usage='denovoeval.py [-h] -c contig.fa -r raw_reads.fa -o output_dict/')
parser.add_argument('-c','--contig',type=str,default=False,help='assembly contigs in .fa format',required=True)
parser.add_argument('-r','--read',type=str,default=False,help='sequencing reads in .fa format',required=True)
parser.add_argument('-o','--outpath',type=str,default='./adenovo_evaluation-out/',help='output directory')

parser.add_argument('--ref',type=str,default=False,help='optional reference genome in .fa format')
parser.add_argument('--ngs',type=str,default=False,help='optional NGS reads for small assembly erroes detection')
parser.add_argument('--depth_plot',action='store_true',help='plot depth at all contigs')
parser.add_argument('-t','--thread',type=int,default=8,help='number of threads')
parser.add_argument('--min_contig_length',type=int,default=10000,help='minimal length for a contig to be considered')
parser.add_argument('--min_contig_length_assemblyerror',type=int,default=1000000,help='minimal contig length for assembly error detection')
parser.add_argument('--min_assembly_error_size',type=int,default=50,help='minimal size for assembly errors')
parser.add_argument('--max_assembly_error_size',type=int,default=4000000,help='maximal size for assembly errors')
parser.add_argument('--skip_structural_error',action='store_true',default=False,help='skip the step of identifying large structural errors')
parser.add_argument('--call_base_error',action='store',help='Force call small assembly errors at basepair level. Results may NOT be confident.')



denovo_args=parser.parse_args()

# check input arguments
if denovo_args.outpath[-1]!='/':
	denovo_args.outpath+='/'
if not os.path.exists(denovo_args.outpath):
	os.mkdir(denovo_args.outpath)


# Simple statistics of contigs
contiginfo=denovo_static.simple(denovo_args.contig,denovo_args.outpath,denovo_args.min_contig_length,denovo_args.min_contig_length_assemblyerror)
chromosomes=contiginfo[0]
chromosomes_large=contiginfo[1]
totalcontiglen=contiginfo[2]
totalcontiglen_large=contiginfo[3]

# map the reads to contigs
os.system("minimap2 -a -t " + str(denovo_args.thread) + "  "+denovo_args.outpath+"valid_contig.fa  " + denovo_args.read + "  > " + denovo_args.outpath+"read_to_contig.sam")
os.system("samtools sort -@ "+str(denovo_args.thread)+" "+denovo_args.outpath+"read_to_contig.sam -o "+str(denovo_args.outpath)+"read_to_contig.sort.bam")
os.system("samtools index "+str(denovo_args.outpath)+"read_to_contig.sort.bam")

# Structural assembly error detection
mapping_info=debreak_detect.detect_sam("read_to_contig.sam",denovo_args.outpath,denovo_args.outpath,chromosomes_large,denovo_args.min_assembly_error_size,denovo_args.max_assembly_error_size)

f=open(denovo_args.outpath+'summary_statistics','a')
f.write('Read to Contig alignment:\n')
f.write('Mapping rate:\t'+str(round(10000*float(mapping_info[4])/(int(mapping_info[4])+int(mapping_info[0])))/10000.0)+'\n')
f.write('Unique mapping rate\t'+str(1-round(10000*float(mapping_info[5])/int(mapping_info[4]))/10000.0)+'\n')
f.write('Depth\t'+str(round(10000*float(mapping_info[6])/totalcontiglen)/10000.0)+'\n')
f.write('Mapping rate in large contigs\t'+str(round(10000*float(mapping_info[1])/(int(mapping_info[1])+int(mapping_info[0])))/10000.0)+'\n')
f.write('Unique mapping rate in large contigs\t'+str(1-round(10000*float(mapping_info[2])/int(mapping_info[1]))/10000.0)+'\n')
f.write('Depth in large conigs\t'+str(round(10000*float(mapping_info[3])/totalcontiglen)/10000.0)+'\n\n\n')
f.close()

depth=round(10000*float(mapping_info[3])/totalcontiglen)/10000.0
minsupp=round(depth/10.0)

#depth=103
#minsupp=round(depth/10.0)  # Should remove

if not denovo_args.skip_structural_error:
	allsvsignal=open(denovo_args.outpath+'read_to_contig.debreak.temp','r').read().split('\n')[:-1]
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

	denovo_static.assembly_info(denovo_args.outpath,chromosomes_large)
	debreak_merge.genotype(depth,denovo_args.outpath,chromosomes_large)

# SNP & indel detection

if denovo_args.call_base_error or denovo_args.ngs:
	if not denovo_args.ngs:
		os.system("freebayes -f "+denovo_args.contig+" "+str(denovo_args.outpath)+"read_to_contig.sort.bam -r "+contiginfo[4]+":"+str(contiginfo[5]/4)+"-"+str(contiginfo[5]/4+100000)+"| vcffilter -f \"QUAL > 20\" > "+str(denovo_args.outpath)+"assembly_basepair_error.vcf")
		denovo_static.basepair_error(denovo_args.outpath)
	else:
		pass


# if reference is also provided:
if denovo_args.ref:
	os.system("minimap2 -a --eqx -t " + str(denovo_args.thread) + " "+denovo_args.ref+" " + denovo_args.contig + "  > " + denovo_args.outpath+"contig_to_ref.sam")
	os.system("samtools sort -@ "+ str(denovo_args.thread) + " " + denovo_args.outpath+"contig_to_ref.sam  -o  " + denovo_args.outpath+"contig_to_ref.sort.bam")
	os.system("samtools index "+ denovo_args.outpath+"contig_to_ref.sort.bam")
	os.system("samtools depth -a "+ denovo_args.outpath+"contig_to_ref.sort.bam  > "+ denovo_args.outpath+"contig_to_ref.depth")

	(chromosomes,totalcontiglen,largestchr,largestlen,coveredlen)=denovo_static.get_ref_chroms(denovo_args.outpath)
	mapping_info=debreak_detect.detect_sam_ref("contig_to_ref.sam",denovo_args.outpath,denovo_args.outpath,denovo_args.min_assembly_error_size,denovo_args.max_assembly_error_size)
	f=open(denovo_args.outpath+'summary_statistics','a')
	f.write('Contig to Reference alignment:\n')
	f.write('Mapping rate\t'+str(round(10000*float(mapping_info[1])/(int(mapping_info[1])+int(mapping_info[0])))/10000.0)+'\n')
	f.write('Unique mapping rate\t'+str(1-round(10000*float(mapping_info[2])/int(mapping_info[1]))/10000.0)+'\n')
	f.write('Depth\t'+str(round(10000*float(mapping_info[3])/totalcontiglen)/10000.0)+'\n')
	f.write('Coverage\t'+str(round(10000*float(coveredlen)/totalcontiglen)/10000.0)+'\n\n\n')
	f.close()
	denovo_static.check_depth_ref(denovo_args.outpath,denovo_args.ref)
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
	
	#contiginfo=['','','','','ctg1']
	denovo_static.basepair_error_ref(denovo_args.outpath,contiginfo[4])

