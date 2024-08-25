# Inspector

A reference-free assembly evaluator.

Author: Maggi Chen

Email: maggic@uab.edu

Current update: Aug. 25, 2024

## Quick Start
```sh
git clone https://github.com/ChongLab/Inspector.git
cd Inspector/
./inspector.py -h

# Evaluate assembly with raw reads
inspector.py -c contig.fa -r rawreads.1.fastq rawreads.2.fastq -o inspector_out/ --datatype clr 
# Evaluate assembly with hifi reads
inspector.py -c contig.fa -r ccsreads.1.fastq ccsreads.2.fastq -o inspector_out/ --datatype hifi

# With reference-based evaluation
inspector.py -c contig.fa -r rawreads.1.fastq --ref reference.fa -o inspector_out/ --datatype clr

# Reference-based only evaluation
inspector.py -c contig.fa -r emptyfile --ref reference.fa -o inspector_out/ 

# Error correction
inspector-correct.py -i inspector_out/ --datatype pacbio-hifi -o inspector_out/

```



## Description

Inspector is a tool for assembly evaluation with long read data. The input includes a contig file, long reads (PacBio CLR, PacBio HiFi, Oxford Nanopore, or mixed platform), and a reference genome (optional). The output includes A summary report, read-to-contig alignment file, a list of structrual errors and small-scale errors. This program was tested on a x86_64 Linux system with a 128GB physical memory.

## Depencency

Dependencies for Inspector:

* python  
* pysam
* statsmodels (tested with version 0.10.1)

* minimap2  (tested with version 2.10 and 2.15)
* samtools  (tested with version 1.9)


Dependencies for Inspector error correction module:
* flye  (tested with version 2.8.3)


## Installation

To create an environment with conda (recommended):
```
conda create --name ins
conda activate ins
conda install -c bioconda inspector

```
Git install after installing all the dependencies. 
```
git clone https://github.com/ChongLab/Inspector.git
export PATH=$PWD/Inspector/:$PATH
```



A subset of human genome assembly is available as testing dataset to validate successful installation. The contig_test.fa includes two contigs (1.4Mbp and 10Kbp). The read_test.fastq.gz includes ~60X PacBio HiFi reads belonging to these two contigs. There are 3 structural errors and 281 small-scale errors present in the testing dataset.
```
cd Inspector/
./inspector.py -c testdata/contig_test.fa -r testdata/read_test.fastq.gz -o test_out/ --datatype hifi 
./inspector-correct.py -i test_out/ --datatype pacbio-hifi 
```
(The Inspector evaluation on testing dataset should finish within several minutes with 4 CPUs and 400MB memory.
The Inspector error correction should finish within 10-15 minutes with 4 CPUs and 500MB memory.)


## General usage


```

inspector.py [-h] -c contig.fa -r raw_reads.fa -o output_dict/
  required arguments:
  --contig,-c           FASTA/FASTQ file containing contig sequences to be evaluated
  --read,-r             A list of FASTA/FASTQ files containing long read sequences

  optional arguments:
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
  --datatype,-d         Input read type. (clr, hifi, nanopore, mixed) [clr]
  --ref                 OPTIONAL reference genome in .fa format
  --thread,-t           Number of threads. [8]
  --min_contig_length   Minimal length for a contig to be evaluated [10000]
  --min_contig_length_assemblyerror    Minimal contig length for assembly error detection. [1000000]
  --pvalue              Maximal p-value for small-scale error identification [0.01 for HiFi, 0.05 for others]
  --skip_read_mapping   Skip the step of mapping reads to contig
  --skip_structural_error       Skip the step of identifying large structural errors
  --skip_structural_error_detect       Skip the step of detecting large structural errors
  --skip_base_error     Skip the step of identifying small-scale errors
  --skip_base_error_detect      Skip the step of detecting small-scale errors from pileup



inspector-correct.py [-h] -i inspector_out/ --datatype pacbio-raw 
  required arguments:
  --inspector,-i        Inspector evaluation directory with original file names
  --datatype            Type of read used for Inspector evaluation. Required for structural error correction
  --outpath,-o          Output directory
  --flyetimeout         Maximal runtime for local assembly with Flye
  --thread,-t           Number of threads
  --skip_structural     Do not correct structural errors. Local assembly will not be performed
  --skip_baseerror      Do not correct small-scale errors
  

```

## Use cases
Inspector evaluates the contigs and identifies assembly errors with sequencing reads. You can use reads from single platform:
```
inspector.py -c contig.fa -r rawreads.1.fastq rawreads.2.fastq -o inspector_out/ --datatype clr
```
Or use a mixed data type:
```
inspector.py -c contig.fa -r rawreads.fastq nanopore.fastq -o inspector_out/ --datatype mixed
```
Reference-based evaluation is also supported:
(Note that reported assembly error from reference-based mode will contain genetic variants)
```
inspector.py -c contig.fa -r rawreads.fastq --ref reference.fa -o inspector_out/ --datatype clr
```
If only the continuity analysis is needed, simply provide an empty file for --read:
```
inspector.py -c contig.fa -r emptyfile -o inspector_out/ --skip_base_error --skip_structural_error
```
For the '--skip' options, do not use unless you are repeating the evaluation with same contig and read files in the same output directory. These may help save time when testing different options for error detection.



Inspector provides an error-correction module to improve assembly accuracy. High-accuracy reads are recommended, especially for small-scale error correction:
(Note that only reads from single platform are supported for error correction.)
```
inspector-correct.py -i inspector_out/ --datatype pacbio-hifi -o inspector_out/ 
```


## Output file descriptions
Inspector writes its evaluation reports into a new directory to avoid conflicts in file names. Inspector error correction uses the evaluationary results to generate corrected assembly. The output directory of Inspector includes:

### summary_statistics
An evaluation report of the input assembly. This file includes the contig continuity statistics reports, the read mapping summary, number of structural and small-scale errors, the QV score, and the contig alignment summary from reference-based mode when available. An assembly with expected total length, high read-to-contig mapping rate, low number of structural and small-scale errors, and high QV score indicates a high assembly quality. When the reference genome is provided, a higher genome coverage and NA50 also indicates more complete assembly.

### structural_error.bed
This file includes all structural errors identified in the assembly. <br />
The first, second and third column indicate the contig, start and end position of structural error. For Expansions and Inversions, the size of error equals the distance between start and end position. For Collapses, the collapsed sequences are missing in the contigs, therefore the EndPosition is StartPosition+1. The length of collapsed sequence should be inferred from the Size column. For HaplotypeSwitches, Inspector normally considers the haplotype containing Expansion-like pattern as haplotype 1, and considers the haplotype with Collapse-like pattern as haplotype 2. The position of error in haplotype1 and haplotype2 are separated by ";" in the second and third columns. <br />
The fourth column indicates the number of error-supporting reads. A high number of error-supporting reads indicates a confident error call. <br />
The fifth and sixth column indicates the type and size of the erorr. For HaplotypeSwitch, the error sizes in the two haplotypes are usually different. Inspector normally lists the size in haplotype1 and haplotype2, corresponding to the position columns. <br />
Column seven to twelve include other information about the structural errors. These are kept for developmental purpose. 

### small_scale_error.bed
This file includes all small-scale errors identified in the assembly. <br />
The first, second and third column indicate the contig, start and end position of small-scale errors. Similar to the structural errors, the distance between StartPosition and EndPosition equals error size for small expansions and equals 1 for small collapses. <br />
The fourth and fifth column indicate the base in the contig and in the reads. <br />
The sixth and seventh column indicate the number of error-supporting reads and the local sequencing depth. A high supporting-read-to-depth ratio means a confident error call. <br />
The eighth column indicates the type of error. <br />
The nineth column indicates the p-value from binominal test. 

### contig_corrected.fa 
The output corrected assembly of Inspector error-correction module. <br />
Only contigs contained in the valid_contig.fa file (longer than --min_contig_length) are corrected. The small-scale errors listed in small_scale_error.bed should all be fixed. The structural errors in structural_error.bed are fixed if the local de novo assembly generates a full-length contig that can be confidently aligned to the original error region. Otherwise, the original sequence will be remained. 



## Citing Inspector
If you use Inspector, please cite
> Chen, Y., Zhang, Y., Wang, A.Y. et al. Accurate long-read de novo assembly evaluation with Inspector. Genome Biol 22, 312 (2021). https://doi.org/10.1186/s13059-021-02527-4
