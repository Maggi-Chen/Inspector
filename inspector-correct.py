#!/usr/bin/env python
import argparse
import multiprocessing
import sys
import denovo_correct as inspector_correct
import os
import time 


t0=time.time()
parser=argparse.ArgumentParser(description='Assembly error correction based on Inspector assembly evaluation', usage='inspector-correct.py [-h] -i inspector_out/ --datatype pacbio-raw ')
parser.add_argument('-v','--version', action='version', version='Inspector_correct_v1.0')
parser.add_argument('-i','--inspector',type=str,default=False,help='Inspector evaluation directory. Original file names are required.',required=True)
parser.add_argument('--datatype',type=str,default=False,help='Type of read used for Inspector evaluation. This option is required for structural error correction when performing local assembly with Flye. (pacbio-raw, pacbio-hifi, nano-raw,pacbio-corr, nano-corr)',required=True)
parser.add_argument('-o','--outpath',type=str,default=False,help='output directory')
parser.add_argument('--skip_structural',action='store_true',default=False,help='Do not correct structural errors. Local assembly will not be performed.')
parser.add_argument('--skip_baseerror',action='store_true',default=False,help='Do not correct base errors.')
parser.add_argument('-t','--thread',type=int,default=8,help='number of threads')

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
inscor_args=parser.parse_args()

if not inscor_args.skip_structural and not inscor_args.datatype:
	print 'Error:  No data type (--datatype) given!\nFor Debreak usage, use -h or --help'
	sys.exit(1)

if inscor_args.datatype not in ['pacbio-raw','pacbio-hifi', 'pacbio-corr', 'nano-raw',' nano-corr']:
	print 'Error:  Data type (--datatype) not valid. Supported read types are: pacbio-raw, pacbio-hifi, pacbio-corr, nano-raw, nano-corr.'
	sys.exit(1)


if inscor_args.inspector[-1]!='/':
	readpath=inscor_args.inspector+'/'
else:
	readpath=inscor_args.inspector

if not inscor_args.outpath:
	outpath=readpath
else:
	if inscor_args.outpath[-1]!='/':
		outpath=inscor_args.outpath+'/'
	else:
		outpath=inscor_args.outpath
if not os.path.exists(outpath):
	os.mkdir(outpath)

t1=time.time()
print 'TIME for validating parameter',t1-t0
try:
	allctg=open(readpath+'valid_contig.fa','r').read().split('>')[1:]
except:
	print 'Error: Contig file not valid. Please keep original file name in the inspector output directory.\nCheck if file is valid: '+readpath+'valid_contig.fa'
	sys.exit(1)
ctginfo={}
for c in allctg:
	ctginfo[c.split('\n')[0]]=c.split('\n')[1]


t2=time.time()
print 'TIME for reading contig and length',t2-t1
newsnplist=[]

if not inscor_args.skip_baseerror:
	try:
		allsnplist=open(readpath+'small_scale_error.bed','r').read().split('\n')[1:-1]
	except:
		print 'Warning: small-scale eror bed file not found. Check file \'small_scale_error.bed\' in Inspector evaluation directory. Continue without small-scale error correction.'
		allsnplist=[]
else:
	allsnplist=[]

snpctg={}

if not inscor_args.skip_structural:
	os.system('mkdir '+outpath+'assemble_workspace/')
	try:
		allaelist=open(readpath+'structural_error.bed','r').read().split('\n')[1:-1]
	except:
		print 'Warning: structural eror bed file not found. Check file \'structural_error.bed\' in Inspector evaluation directory. Continue without structural error correction.'
		allaelist=[]
else:
	allaelist=[]

snpctg={}
aectg={}
for ctgname in ctginfo:
	snpctg[ctgname]=[]
	aectg[ctgname]=[]
for aeinfo in allsnplist:
	snpctg[aeinfo.split('\t')[0]]+=[aeinfo]
for aeinfo in allaelist:
	aectg[aeinfo.split('\t')[0]]+=[aeinfo]

allsnplist=[]
allaelist=[]
bamfile=readpath+'read_to_contig.bam'

t3=time.time()
print 'TIME for reading assembly errors',t3-t2

for chrominfo in ctginfo:
	inspector_correct.error_correction_large(chrominfo,ctginfo[chrominfo],aectg[chrominfo],snpctg[chrominfo],bamfile,outpath,inscor_args.datatype,inscor_args.thread/3)

t4=time.time()
print 'TIME for correcting all contigs',t4-t3

f=open(outpath+'contig_corrected.fa','w')
for chrominfo in ctginfo:
	try:
		correctedinfo=open(outpath+'contig_corrected_'+chrominfo+'.fa','r').read()
		f.write(correctedinfo)
	except:
		print 'Warning: corrected contig ',chrominfo,'not found.'

f.close()
t5=time.time()
print 'TIME for writing corrected contig',t5-t4
os.system('rm '+outpath+'contig_corrected_*fa')

print 'Error correction DONE.'


