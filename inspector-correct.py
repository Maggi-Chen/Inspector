#!/usr/bin/env python
import argparse
import multiprocessing
import sys
import denovo_correct as inspector_correct
import os
from datetime import datetime
import time 


t0=time.time()
parser=argparse.ArgumentParser(description='Assembly error correction based on Inspector assembly evaluation', usage='inspector-correct.py [-h] -i inspector_out/ --datatype pacbio-raw ')
parser.add_argument('-v','--version', action='version', version='Inspector_correct_v1.0')
parser.add_argument('-i','--inspector',type=str,default=False,help='Inspector evaluation directory. Original file names are required.',required=True)
parser.add_argument('--datatype',type=str,default=False,help='Type of read used for Inspector evaluation. This option is required for structural error correction when performing local assembly with Flye. (pacbio-raw, pacbio-hifi, nano-raw,pacbio-corr, nano-corr)',required=True)
parser.add_argument('-o','--outpath',type=str,default=False,help='output directory')
parser.add_argument('--flyetimeout',type=int,default=1200,help='Maximal runtime for local assembly with Flye. Unit is second. [1200]')
parser.add_argument('--skip_structural',action='store_true',default=False,help='Do not correct structural errors. Local assembly will not be performed.')
parser.add_argument('--skip_baseerror',action='store_true',default=False,help='Do not correct base errors.')
parser.add_argument('-t','--thread',type=int,default=8,help='number of threads')

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

inscor_args=parser.parse_args()
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


logf=open(outpath+'Inspector_correct.log','a')
logf.write('Inspector assembly error correction starting... '+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+'\n')

if not inscor_args.skip_structural and not inscor_args.datatype:
	logf.write('Error:  No data type (--datatype) given!\nFor Debreak usage, use -h or --help\n')
	sys.exit(1)

if inscor_args.datatype not in ['pacbio-raw','pacbio-hifi', 'pacbio-corr', 'nano-raw','nano-corr']:
	logf.write('Error:  Data type (--datatype) not valid. Supported read types are: pacbio-raw, pacbio-hifi, pacbio-corr, nano-raw, nano-corr.\n')
	sys.exit(1)


t1=time.time()
logf.write('TIME for validating parameter'+str(t1-t0)+'\n')

try:
	allctg=open(readpath+'valid_contig.fa','r').read().split('>')[1:]
except:
	logf.write('Error: Contig file not valid. Please keep original file name in the inspector output directory.\nCheck if file is valid: '+readpath+'valid_contig.fa\n')
	sys.exit(1)
ctginfo={}
for c in allctg:
	ctginfo[c.split('\n')[0]]=c.split('\n')[1]


t2=time.time()
logf.write('TIME for reading contig and length'+str(t2-t1)+'\n')
newsnplist=[]

if not inscor_args.skip_baseerror:
	try:
		allsnplist=open(readpath+'small_scale_error.bed','r').read().split('\n')[1:-1]
	except:
		logf.write('Warning: small-scale eror bed file not found. Check file \'small_scale_error.bed\' in Inspector evaluation directory. Continue without small-scale error correction.\n')
		allsnplist=[]
else:
	allsnplist=[]

snpctg={}

if not inscor_args.skip_structural:
	os.system('mkdir '+outpath+'assemble_workspace/')
	try:
		allaelist=open(readpath+'structural_error.bed','r').read().split('\n')[1:-1]
	except:
		logf.write('Warning: structural eror bed file not found. Check file \'structural_error.bed\' in Inspector evaluation directory. Continue without structural error correction.\n')
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
logf.write('TIME for reading assembly errors'+str(t3-t2)+'\n')
logf.close()


for chrominfo in ctginfo:
	inspector_correct.error_correction_large(chrominfo,ctginfo[chrominfo],aectg[chrominfo],snpctg[chrominfo],bamfile,outpath,inscor_args.datatype,inscor_args.thread//3,inscor_args.flyetimeout)

t4=time.time()
logf=open(outpath+'Inspector_correct.log','a')
logf.write('TIME for correcting all contigs'+str(t4-t3)+'\n')
logf.close()

f=open(outpath+'contig_corrected.fa','w')
for chrominfo in ctginfo:
	try:
		correctedinfo=open(outpath+'contig_corrected_'+chrominfo+'.fa','r').read()
		f.write(correctedinfo)
	except:
		logf=open(outpath+'Inspector_correct.log','a')
		logf.write('Warning: corrected contig ',chrominfo,'not found.\n')
		logf.close()
f.close()
t5=time.time()
logf=open(outpath+'Inspector_correct.log','a')
logf.write('TIME for writing corrected contig'+str(t5-t4)+'\n')
os.system('rm '+outpath+'contig_corrected_*fa')
logf.write('Inspector error correction finished. Bye.\n')
logf.close()


