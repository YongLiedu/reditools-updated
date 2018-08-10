# REDItools: python scripts for RNA editing detection by RNA-Seq data

<!-- TOC START min:2 max:4 link:true update:true -->
- [Introduction](#introduction)
- [Publications](#publications)
- [Download](#download)
  - [Download REDItools](#download-reditools)
  - [Download testing dataset](#download-testing-dataset)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input files](#input-files)
  - [Preparing input files](#preparing-input-files)
    - [FASTA](#fasta)
    - [BAM](#bam)
    - [GTF](#gtf)
    - [TAB](#tab)
    - [SpliceSites](#splicesites)
- [Output files](#output-files)
- [Usage Information](#usage-information)
  - [REDItoolDnaRna.py](#reditooldnarnapy)
  - [REDItoolKnown.py](#reditoolknownpy)
  - [REDItoolDenovo.py](#reditooldenovopy)
- [Additional scripts:](#additional-scripts)
  - [REDItoolBlatCorrection.py](#reditoolblatcorrectionpy)
  - [FilterTable.py](#filtertablepy)
  - [AnnotateTable.py](#annotatetablepy)
  - [SearchInTable.py](#searchintablepy)
  - [selectPositions.py](#selectpositionspy)
  - [SortTable.py (new in version 1.0.3)](#sorttablepy-new-in-version-103)
  - [tableToTabix.py (new in version 1.0.3)](#tabletotabixpy-new-in-version-103)
  - [SortGFF.py (new in version 1.0.3)](#sortgffpy-new-in-version-103)
  - [GFFtoTabix.py (new in version 1.0.3)](#gfftotabixpy-new-in-version-103)
  - [TableToGFF.py (new in version 1.0.3)](#tabletogffpy-new-in-version-103)
- [Usage Example](#usage-example)
  - [Test dataset](#test-dataset)
- [Contact](#contact)

<!-- TOC END -->

## Introduction
REDItools are python scripts developed with the aim to study RNA editing at genomic scale by next generation sequencing data. RNA editing is a post-transcriptional phenomenon involving the insertion/deletion or substitution of specific bases in precise RNA localizations.
In human, RNA editing occurs by deamination of cytosine to uridine (C-to-U) or mostly by the adenosine to inosine (A-to-I) conversion through ADAR enzymes. A-to-I substitutions may have profound functional consequences and have been linked to a variety of human diseases including neurological and neurodegenerative disorders or cancer. Next generation sequencing technologies offer the unique opportunity to investigate in depth RNA editing even though no dedicated software has been released up to now.

REDItools are simple python scripts conceived to facilitate the investigation of RNA editing at large-scale and devoted to research groups that would to explore such phenomenon in own data but don’t have sufficient bioinformatics skills.
They work on main operating systems (although unix/linux-based OS are preferred), can handle reads from whatever platform in the standard BAM format and implement a variety of filters.

## Publications
* *Picardi E, D'Erchia AM, Gallo A and Pesole G* **Detection of post-transcriptional RNA editing events.** *Methods Mol Biol. 2015;1269:189-205* PubMed PMID: [25577380](http://www.ncbi.nlm.nih.gov/pubmed/25577380)

* *Picardi E, D'Erchia AM, Gallo A, Montalvo A and Pesole G.* **Uncovering RNA Editing Sites in Long Non-Coding RNAs.** *Front Bioeng Biotechnol. 2014 Dec 5;2:64* PubMed PMID: [25538940](http://www.ncbi.nlm.nih.gov/pubmed/25538940)

* *Picardi E and Pesole G.* **REDItools: high-throughput RNA editing detection made easy.** *Bioinformatics. 2013 Jul 15;29(14):1813-4* PubMed PMID: [23742983](http://www.ncbi.nlm.nih.gov/pubmed/23742983)

## Download
### Download REDItools
* [REDItools v1.0.4](http://sourceforge.net/projects/reditools/files/REDItools-1.0.4.tar.gz/download)
* [REDItools v1.0.3](http://sourceforge.net/projects/reditools/files/REDItools-1.0.3.tar.gz/download)
* [REDItools v1.0.2](http://sourceforge.net/projects/reditools/files/REDItools-1.0.2.tar.gz/download)

### Download testing dataset
* [testREDItools.tar.gz](http://srv00.ibbe.cnr.it/reditools/data/testREDItools.tar.gz)

## Prerequisites
REDItools require python 2.6 or 2.7 available at [python web-site](http://www.python.org).
All scripts have not been tested on python 3.
In addition, two external modules need to be installed to run and use REDItools, even though the second module is optional:

* pysam module version 0.6 or 0.7 available on [Google Code](http://code.google.com/p/pysam/)
* fisher module version 0.1.4 (optional) available at [python web site](http://pypi.python.org/pypi/fisher/)

To perform Blat correction and format alignment exchange (SAM to BAM and vice versa) the
following packages should be installed or already present in your path:

* Blat package including gfServer and gfClient executables, available from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/admin/exe/)
* Samtools and tabix from [sourceforge](http://samtools.sourceforge.net/)

## Installation
REDItools require pysam module. To check its presence in your system try the
command:
```bash
 python -c 'import pysam'
```
If no errors are raised, the module is already installed. Alternatively, download a copy
from the above link and follow instructions inside.

REDItools are stand-alone scripts and each one can be simply launched by:
```bash
 python REDItoolscript.py
```
Alternatively you can easily install REDItools by the following commands:
```bash
 gunzip reditools-1.0.2.tar.gz
 tar –xvf reditools-1.0.2.tar
 cd reditools-1.0.2
 python setup.py install
```
To install scripts in a specific location the --prefix option can be used:
```bash
 python setup.py install --prefix=/my/path
```
All scripts will be located in /my/path/bin.

## Input files
REDItools have been designed to handle mapped reads from next generation sequencing
technologies in the standard BAM format.
The following input files are accepted:

* BAM files containing RNA-Seq or DNA-Seq reads aligned onto a reference genome (mandatory);
* FASTA file of reference genome (mandatory);
* GTF file containing specific genomic positions (optional);
* TAB file containing individual genomic positions (optional);
* SpliceSites file containing splice sites information (optional).

### Preparing input files
#### FASTA
Reference genome in fasta format needs to be indexed using samtools:
```bash
 samtools faidx reference.fa
```
#### BAM
Starting from alignments in SAM format, you first need to convert SAM to BAM:
```bash
 samtools view –bS –T reference.fa myfile.sam > myfile.bam
```
then sort your BAM:
```bash
 samtools sort myfile.bam myfile.sorted
```
and finally index your sorted BAM:
```bash
 samtools index myfile.sorted.bam
```
#### GTF
GTF files for specific annotations and genomic regions can be downloaded from [UCSC](http://genome.ucsc.edu/).
Alternatively, UCSC tables in psl format may be used and converted in GTF by the following awk/gawk commands:

* For RepeakMask table:
```bash
  $ gunzip rmsk.txt.gz
  $ gawk 'OFS="\t"{print $6,"rmsk_hg19",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt > rmsk.gtf
```
* For simpleRepeat table:
```bash
  $ gunzip simpleRepeat.txt.gz
  $ gawk 'OFS="\t"{print $2,"trf_hg19","trf",$3+1,$4,".","+",".","gene_id \""$5"\"; transcript_id \""$17"\";"}' simpleRepeat.txt > simpleRepeat.gtf
```
* For RefSeq annotations you can download from UCSC and use the [genePredToGtf](http://hgdownload.soe.ucsc.edu/admin/exe/) utility:
```bash
  $ gunzip refGene.txt.gz
  $ cut -f 2- refGene.txt | genePredToGtf -utr -source=hg19_refseq file stdin refGene.gtf
```
Other annotations can be generated accordingly. Please, consider that the GTF format includes 9 Tab-delimited fields:

* chromosome/region name
* source (for example hg19_refseq)
* feature (for example CDS, exon, snp ...)
* start coordinate (1-based)
* end coordinate (1-based)
* score (use a dot if unknown)
* strand (+ or - use + if unknown)
* frame for CDS only or a dot (a dot also if unknown)
* attributes. For REDItools there are two mandatory attributes: gene_id and transcript_id

Following lines show a GTF example for RefSeq annotations:
```bash
 $ head -n6 refGene.gtf
 chr7    hg19_refseq     5UTR    139025878       139026130       .       +       .       gene_id "C7orf55"; transcript_id "NM_197964";
 chr7    hg19_refseq     CDS     139026131       139026268       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
 chr7    hg19_refseq     CDS     139030247       139030447       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
 chr7    hg19_refseq     3UTR    139030451       139031065       .       +       .       gene_id "C7orf55"; transcript_id "NM_197964";
 chr7    hg19_refseq     start_codon     139026131       139026133       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
 chr7    hg19_refseq     stop_codon      139030448       139030450       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
```
GTF files can be sorted by the command:
```bash
 $ sort -k1,1 -k4,4n mygft > mygtf_sorted
```
#### TAB
TAB files are simple textual files with at least three tabulated columns including:

* genomic region (generally the chromosome name according to the reference genome)
* coordinate of the position (1-based)
* strand (+ or -). You can also indicate strand by 0 (strand -), 1 (strand +) or 2 (+ and - or unknown)

==============  ==========  ======
genomic region  coordinate  strand
==============  ==========  ======
chr21           10205589    \-
chr21           10205629    \-
chr21           15411496    \+
chr21           15412990    \+
chr21           15414553    \+
chr21           15415901    \+
chr21           15417667    \+
chr21           15423330    \+
==============  ==========  ======

TAB files must be coordinate sorted. In unix/linux environment they can be sorted by the sort command:
```bash
 sort -k1,1 -k2,2n mytable.txt > mytable.sorted.txt
```
#### SpliceSites
SpliceSites files are simple textual file including 5 columns separated by space:

* genomic region (chromosome/region name according to reference and input BAMs)
* start coordinate of splice site (1-based)
* end coordinate of splice site (1-based)
* splice site type, A for acceptor and D for donor
* strand (+ or -)

==============  ==========  =====   ====  ======
genomic region  start       end     type  strand
==============  ==========  =====   ====  ======
chr1            13224       13225   A     \+
chr1            12227       12228   D     \+
chr1            12594       12595   A     \+
chr1            12721       12722   D     \+
chr1            13402       13403   A     \+
chr1            13655       13656   D     \+
chr1            29320       29321   D     \-
chr1            24891       24892   A     \-
chr1            24737       24738   D     \-
chr1            18379       18380   A     \-
==============  ==========  =====   ====  ======

SpliceSites files can be obtained starting from UCSC annotations in psl format and using
the perl script psl_splicesites included in the [GMAP/GSNAP](http://research-pub.gene.com/gmap/) package:
```bash
 gunzip -c refGene.txt.gz | psl_splicesites -s 1 > mysplicesites
```
The format of mysplicesites is:

 >NM_005236.exon11/11 chr16:14041470..14041471 acceptor 2778
 >NM_005235.exon1/28 chr2:213403173..213403172 donor 413544
 >NM_005235.exon2/28 chr2:212989629..212989628 acceptor 413544
 >NM_005235.exon2/28 chr2:212989477..212989476 donor 177135
 >NM_005235.exon3/28 chr2:212812342..212812341 acceptor 177135
 >NM_005235.exon3/28 chr2:212812155..212812154 donor 159270
 >NM_005235.exon4/28 chr2:212652885..212652884 acceptor 159270
 >NM_005235.exon4/28 chr2:212652750..212652749 donor 37320
 >NM_005235.exon5/28 chr2:212615430..212615429 acceptor 37320

Then, the following awk/gawk command line can be used to get the final SpliceSite file:
```bash
	$ gawk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' mysplicesites > mysplicesites.ss
	$ more mysplicesites.ss
	chr16 14041470 14041471 A +
	chr2 213403172 213403173 D -
	chr2 212989628 212989629 A -
	chr2 212989476 212989477 D -
	chr2 212812341 212812342 A -
	chr2 212812154 212812155 D -
	chr2 212652884 212652885 A -
	chr2 212652749 212652750 D -
	chr2 212615429 212615430 A -
```
## Output files
REDItools print out results in simple textual tables:

======  ========  =========  ======  ============  =====  ==================  =======  =========
Region  Position  Reference  Strand  Coverage-q25  MeanQ  BaseCount[A,C,G,T]  AllSubs  Frequency
======  ========  =========  ======  ============  =====  ==================  =======  =========
chr21   15412990  A          1       18            37.22   [3, 0, 15, 0]      AG       0.83
chr21   15415901  A          1       13            37.15   [2, 0, 11, 0]      AG       0.85
chr21   15423330  A          1       11            38.27   [4, 0, 7, 0]       AG       0.64
chr21   15425640  A          1       8             36.12   [0, 0, 8, 0]       AG       1.00
chr21   15456434  T          1       90            34.96   [0, 6, 1, 83]      TC TG    0.07
chr21   15461406  A          1       83            37.27   [73, 0, 10, 0]     AG       0.12
chr21   15461417  A          1       90            36.26   [72, 0, 18, 0]     AG       0.20
chr21   15461444  A          1       64            37.22   [26, 0, 38, 0]     AG       0.59
chr21   15461479  A          1       70            36.96   [66, 0, 4, 0]      AG       0.06
chr21   15461486  A          1       68            37.06   [61, 0, 7, 0]      AG       0.10
chr21   15461503  A          1       76            37.26   [69, 0, 7, 0]      AG       0.09
chr21   15461511  A          1       81            37.68   [55, 0, 26, 0]     AG       0.32
======  ========  =========  ======  ============  =====  ==================  =======  =========

where:
 * Region: is the genomic region according to reference
 * Position: is the exact genomic coordinate (1-based)
 * Reference: is the nucleotide base in reference genome
 * Strand: is strand information with notation 1 for + strand, 0 for - strand and 2 for unknown or not defined strand
 * Coverage-qxx: is the depth per site at a given xx quality score (min. value)
 * MeanQ: is the mean quality score per site
 * BaseCount[A,C,G,T]: is the base distribution per site in the order A, C, G and T
 * AllSubs: is the list of observed substitutions at a given site, separated by a space. A character "-" is included in case of invariant sites.
 * Frequency: is the observed frequency of substitution. In case of multiple substitutions, it refers to the first in the AllSubs field.
   In the table above, for example, at line 5 the frequency 0.07 is linked to TC substitution and not to TG that will be always lower than this value.

REDItoolDnaRna.py includes five additional columns to take into account information from DNA-Seq reads. Such columns are indicated as:
 * gCoverage-q25: is the depth per site at a given xx quality score (min. value) in DNA-Seq
 * gMeanQ: is the mean quality score per site in DNA-Seq
 * gBaseCount[A,C,G,T]: is the base distribution per site in the order A, C, G and T
 * gAllSubs: is the list of observed substitutions at a given site, separated by a space. A character "-" is included in case of invariant sites.
 * gFrequency: is the observed frequency of substitution. In case of multiple substitutions, it refers to the first in the gAllSubs field.

In case of positions not supported by DNA-Seq reads, a character "-" will be added to each extra column.

REDItoolDenovo.py includes one additional column concerning the reliability of editing prediction.
 * Pvalue: is the pvalue per site calculated according to Fisher exact test. It indicates how much the observed base distribution for a change is different from the expected,
   calculated by the empirical base substitution for the entire RNA-Seq experiment.

**Note** that BaseCount and Reference column are modified according to Strand column. In case of Strand value of 2, the genomic Reference is used.

## Usage Information
### REDItoolDnaRna.py
REDItoolDnaRna.py is the main script devoted to the identification of RNA editing events
taking into account the combined information from RNA-Seq and DNA-Seq data in BAM format.
To look at potential RNA editing candidates, RNA-Seq data alone can be used.

Options:
|   option   |   description  |
| :--------: | :------------- |
| -i	|	RNA-Seq BAM file
| -j	|	DNA-Seq BAM files separated by comma or folder containing BAM files. **Note** that each chromosome/region must be present in a single BAM file only.
| -I	|	Sort input RNA-Seq BAM file
| -J	|	Sort input DNA-Seq BAM file
| -f	|	Reference file in fasta format. **Note** that chromosome/region names in the reference must match chromosome/region names in BAMs files.
| -C	|	Base interval to explore [100000]. It indicates how many bases have to be loaded during the run.
| -k	|	List of chromosomes to skip separated by comma or file (each line must contain a chromosome/region name).
| -t	|	Number of threads [1]. It indicates how many processes should be launched. Each process will work on an individual chromosome/region.
| -o	|	Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.
| -F	|	Internal folder name [null] is the main folder containing output tables.
| -M	|	Save a list of columns with quality scores. It produces at most two files in the pileup-like format.
| -c	|	Minimum read coverage (dna,rna) [10,10]
| -Q	|	Fastq offset value (dna,rna) [33,33]. For Illumina fastq 1.3+ 64 should be used.
| -q	|	Minimum quality score (dna,rna) [25,25]
| -m	|	Minimum mapping quality score (dna,rna) [25,25]
| -O	|	Minimum homoplymeric length (dna,rna) [5,5]
| -s	|	Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA; 2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.
| -g	|	Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)
| -x	|	Strand confidence [0.70]
| -S	|	Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.
| -G	|	Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.
| -K	|	GFF File with positions to exclude (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.
| -T	|	Work only on given GFF positions (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.
| -X	|	Sort annotation files. It requires grep and sort unix executables.
| -e	|	Exclude multi hits in RNA-Seq
| -E	|	Exclude multi hits in DNA-Seq
| -d	|	Exclude duplicates in RNA-Seq
| -D	|	Exclude duplicates in DNA-Seq
| -p	|	Use paired concardant reads only in RNA-Seq
| -P	|	Use paired concardant reads only in DNA-Seq
| -u	|	Consider mapping quality in RNA-Seq
| -U	|	Consider mapping quality in DNA-Seq
| -a	|	Trim x bases up and y bases down per read [0-0] in RNA-Seq
| -A	|	Trim x bases up and y bases down per read [0-0] in DNA-Seq
| -b	|	Blat folder for correction in RNA-Seq
| -B	|	Blat folder for correction in DNA-Seq
| -l	|	Remove substitutions in homopolymeric regions in RNA-Seq
| -L	|	Remove substitutions in homopolymeric regions in DNA-Seq
| -v	|	Minimum number of reads supporting the variation [3] for RNA-Seq
| -n	|	Minimum editing frequency [0.1] for RNA-Seq
| -N	|	Minimum variation frequency [0.1] for DNA-Seq
| -z	|	Exclude positions with multiple changes in RNA-Seq
| -Z	|	Exclude positions with multiple changes in DNA-Seq
| -W	|	Select RNA-Seq positions with defined changes (separated by comma ex: AG,TC) [default all]
| -R	|	Exclude invariant RNA-Seq positions
| -V	|	Exclude sites not supported by DNA-Seq
| -w	|	File containing splice sites annotations (SpliceSite file format see above for details)
| -r	|	Num. of bases near splice sites to explore [4]
| --gzip |	Gzip output files
| -h, --help	|	Print the help

Example:
```bash
 REDItoolDnaRna.py -i rnaseq.bam -j dnaseq.bam -f myreference.fa -o myoutputfolder
```
### REDItoolKnown.py
REDItoolKnown.py has been developed to explore the RNA editing potential of RNA-Seq data
sets using known editing events. Such events can be downloaded from DARNED database or
generated from supplementary materials of a variety of publications. Known RNA editing events
have to be stored in TAB files (see above for details).

Options:
|   option   |   description  |
| :--------: | :------------- |
| -i	|	BAM file
| -I	|	Sort input BAM file
| -f	|	Reference in fasta file
| -l	|	List of known RNA editing events
| -C	|	Base interval to explore [100000]
| -k	|	List of chromosomes to skip separated by comma or file
| -t	|	Number of threads [1]
| -o	|	Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.
| -F	|	Internal folder name [null] is the main folder containing output tables.
| -c	|	Min. read coverage [10]
| -Q	|	Fastq offset value [33]
| -q	|	Minimum quality score [25]
| -m	|	Minimum mapping quality score [25]
| -O	|	Minimum homoplymeric length [5]
| -s	|	Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA; 2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.
| -g	|	Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)
| -x	|	Strand confidence [0.70]
| -S	|	Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.
| -G	|	Infer strand by GFF annotation (must be sorted, otherwise use -X). Sorting requires grep and sort unix executables.
| -X	|	Sort annotation files. It requires grep and sort unix executables.
| -K	|	File with positions to exclude (chromosome_name coordinate)
| -e	|	Exclude multi hits
| -d	|	Exclude duplicates
| -p	|	Use paired concardant reads only
| -u	|	Consider mapping quality
| -T	|	Trim x bases up and y bases down per read [0-0]
| -B	|	Blat folder for correction
| -U	|	Remove substitutions in homopolymeric regions
| -v	|	Minimum number of reads supporting the variation [3]
| -n	|	Minimum editing frequency [0.1]
| -E	|	Exclude positions with multiple changes
| -P	|	File containing splice sites annotations (SpliceSite file format see above for details)
| -r	|	Num. of bases near splice sites to explore [4]
| -h	|	Print the help

Example:
```bash
 REDItoolKnown.py -i rnaseq.bam -f reference.fa -l knownEditingSites.tab
```
### REDItoolDenovo.py
REDItoolDenovo.py has been conceived to predict potential RNA editing events using RNA-Seq
data alone and without any a priori knowledge about genome information and biological properties
of RNA editing phenomenon.

Options:
|   option   |   description  |
| :--------: | :------------- |
| -i	|	BAM file
| -I	|	Sort input BAM file
| -f	|	Reference in fasta file
| -k	|	List of chromosomes to skip separated by comma or file
| -t	|	Number of threads [1]
| -o	|	Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.
| -F	|	Internal folder name [null] is the main folder containing output tables.
| -b	|	Use input distribution file
| -a	|	Fisher Tail [l, r, t] [default l] [l=left_tail, r=right_tail, t=two_tail]
| -c	|	Min. read coverage [10]
| -Q	|	Fastq offset value [33]
| -q	|	Min. quality score [25]
| -m	|	Min. mapping quality score [25]
| -O	|	Min. homoplymeric length [5]
| -s	|	Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA; 2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.
| -g	|	Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)
| -x	|	Strand confidence [0.70]
| -S	|	Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.
| -G	|	Infer strand by gff annotation (must be sorted, otherwise use -X)
| -X	|	Sort annotation files
| -K	|	File with positions to exclude
| -e	|	Exclude multi hits
| -d	|	Exclude duplicates
| -l	|	Select significant sites
| -V	|	Significant value [0.05]
| -w	|	Statistical test [BH, BO, NO] [default BH] [BH=Benjamini, BO=Bonferroni, NO=No Correction]
| -U	|	Use specific substitutions separated by comma [example: AG,TC]
| -p	|	Use paired concardant reads only
| -u	|	Consider mapping quality
| -T	|	Trim x bases up and y bases down per read [0-0]
| -B	|	Blat folder for correction
| -W	|	Remove substitutions in homopolymeric regions
| -v	|	Minimum number of reads supporting the variation [3]
| -n	|	Minimum editing frequency [0.1]
| -E	|	Exclude positions with multiple changes
| -P	|	File containing splice sites annotations (SpliceSite file format see above for details)
| -r	|	Number of bases near splice sites to explore [4]
| -h	|	Print the help

Example:
```bash
 REDItoolDenovo.py -i rnaseq.bam -f reference.fa
```
## Additional scripts:
### REDItoolBlatCorrection.py
REDItoolBlatCorrection.py requires gfServer and gfClient programs from [Blat package](http://hgdownload.cse.ucsc.edu/admin/exe/).
Reference fasta file can be converted in .2bit format using the utility [faToTwoBit](http://hgdownload.cse.ucsc.edu/admin/exe/).

Options:
|   option   |   description  |
| :--------: | :------------- |
|   -i	|	BAM file
|   -I	|	Sort input BAM file
|   -f	|	Genomic fasta file
|   -F	|	Genomic fasta file in 2bit format for gfServer
|   -t	|	Num. of working threads [1]
|   -o	|	Output Folder [BlatCorrection_XXXX]
|   -k	|	List of chromosomes to skip separated by comma or file
|   -r	|	Regions in GFF in which Blat will be launched
|   -s	|	Sort GFF (if unsorted). It requires grep and sort unix executables.
|   -q	|	Minimum quality score [25]
|   -Q	|	Fastq offset value [33]
|   -V	|	Verify if gfServer is alive
|   -T	|	Stop gfServer at script end
|   -h	|	Print the help

Example:
```bash
 REDItoolBlatCorrection.py -i rnaseq.bam -f reference.fa -F reference.2bit -o BlatCorrection -V -T
```
At the end of the run, in the BlatCorrection folder, there will be different .bad files, one per chromosome/region.
Each .bad file includes two space separated columns:

* read name
* a number from 0 to 2 indicating the read mate; 1: if read is first in pair, 2: if read is second in pair, 0: if read is a single end

And the content looks like:

 TUPAC_0006:2:41:10999:19783#0 1
 TUPAC_0006:2:41:10999:19783#0 2
 TUPAC_0006:3:54:8655:18923#0 2
 PRESLEY_0005:1:73:13079:17565#0 1
 TUPAC_0006:2:9:12695:14552#0 1
 PRESLEY_0005:1:73:13079:17565#0 1
 TUPAC_0006:2:9:12695:14552#0 1
 PRESLEY_0005:1:12:1152:17918#0 2
 SINATRA_0006:7:15:8730:10887#0 1
 SINATRA_0006:7:67:12736:11713#0 1
 SINATRA_0006:7:50:9125:3151#0 1
 SINATRA_0006:8:49:4592:20505#0 1
 SINATRA_0006:7:14:10225:3766#0 1
 PRESLEY_0005:1:118:4480:20145#0 1
 SINATRA_0006:7:6:19272:9901#0 1

### FilterTable.py
FilterTable.py filters positions of a input table according to specific annotations indexed by tabix tool.
Filtered out positions will be marked with "#" at the beginning of each line. To exclude such lines the
option -E should be used. Features are the same as indicated in the third field of GTF annotation file.

Options:
|   option   |   description  |
| :--------: | :------------- |
|  -i        |      Table file
|  -f        |      Sorted file with positions to filter in (GTF sorted file as above)
|  -s        |      Sorted file with positions to filter out (GTF sorted file as above)
|  -F        |      Features to filter in (separated by comma)
|  -S        |      Features to filter out (separated by comma)
|  -E        |      Exclude positions filtered out
|  -o        |      Save filtered lines to a file [stdout]
|  -p        |      Print simple statistics
|  -h        |      Print the help

Example:
```bash
 FilterTable.py -i mytable -s dbsnp137.gtf.gz -S snp -o mytable_filtered -E -p
```
### AnnotateTable.py
AnnotateTable.py annotates individual positions of a table file according to annotations indexed by tabix tool.
It adds from 1 to 3 additional columns according to -c option.

Options:
|   option   |   description  |
| :--------: | :------------- |
| -a	|	Sorted Annotation file in GTF format
| -i	|	Annotate a file of positions [column1=region, column2=coordinate (1 based)] or a single position [region:coordinate (1 based)]
| -s	|	Strand column in annotation file [4]
| -u	|	Not use table strand info (fix it to 2)
| -c	|	Add columns separated by comma (feature:1, gene_id:2, transcript_id:3) [1,2]
| -n	|	Column name [Col]
| -S	|	Correct strand by annotation
| -C	|	Columns with base distribution [7,12] (in combination with -S)
| -o	|	Save lines to a file
| -h	|	Print the help

Example:
```bash
 AnnotateTable.py -i mytable -a rmsk.gtf.gz -u -c1,2,3 -n rmsk -o mytable_rmsk
```
### SearchInTable.py
SearchInTable.py looks for individual positions in a list or table of positions.

Options:
|   option   |   description  |
| :--------: | :------------- |
|  -i        |      Sorted table file (first col=reference; second col=coordinate 1 based) or tabix indexed table (ending with .gz) (TAB format is required)
|  -j        |      Query (file or single positions as chr21:123456)
|  -C        |      Sequence name column [1]
|  -S        |      Start column [2]
|  -E        |      End column; can be identical to '-S' [2]
|  -p        |      Print position header (like a fasta header >chr21:123456)
|  -n        |      Print "Not found"
|  -s        |      Print simple statistics on standard error
|  -o        |      Save found/not found positions on file
|  -h        |      Print this help

Example:
```bash
 SearchInTable.py -i mytable.gz -j mypositions -s -o results
```
### selectPositions.py
selectPositions.py can filter an output REDItool table according to different criteria.

Options:
|   option   |   description  |
| :--------: | :------------- |
| -i         |     Table file from REDItools
| -d         |     Base distribution column for DNA-Seq (-1: no DNA-Seq) [-1]
| -c         |     Coverage RNA-Seq [5]
| -C         |     Coverage DNA-Seq [5]
| -v         |     Bases supporting RNA-Seq variation [1]
| -V         |     Bases supporting DNA-Seq variation [0]
| -s         |     Substitutions to select in RNA-Seq (separated by comma AG,CT) [all]
| -f         |     Frequency of variation in RNA-Seq [0.1]
| -F         |     Frequency of non-variation in DNA-Seq [0.95]
| -e         |     Exclude multiple substitutions in RNA-Seq
| -r         |     Exclude invariant sites in RNA-Seq
| -R         |     Exclude variant sites in DNA-Seq #
| -u         |     Use only positions supported by DNA-Seq
| -o         |     Save selected positions on outTable_533864766
| -h         |     Print this help

### SortTable.py (new in version 1.0.3)
SortTable.py is a facility to sort an input table according to genomic region and coordinates.
Input table may be TAB-delimited or use whatever delimiter. This script has been introduced for
users working on operating systems in which the sort program is not present. On very large input files
sorting time may increase respect to the unix sort program. On common GFFs including gene annotations
sorting time is lower or equal than the sort program.
Optionally, an input table can be outputted as TAB-delimited. Please tune the -b option according to
your memory capability. Default value should work well for many computers.

Options:
|   option   |   description  |
| :--------: | :------------- |
| -i	|	Table file
| -d	|	Delimiter character [\t] (default TAB)
| -s	|	Sequence name column [1]
| -c	|	Start column [4]
| -e	|	End column (can be identical to -c) [5]
| -m	|	Skip lines starting with [#]
| -o	|	Sorted output file [sortedTable_%s]
| -O	|	Output as TAB-delimited
| -b	|	Buffer size (as number of lines) [32000]
| -t	|	Temporary directory to use (multiple -t may be used)
| -h	|	Print this help

Example for input table with space-based columns:
```bash
 SortTable.py -i mytable.txt -d ' ' -s 1 -c 2 -e 2 -O
```
### tableToTabix.py (new in version 1.0.3)
tableToTabix.py create a tabix table and by default sort the input table. Useful to generate tabix tables
from GFF file to prepare in advance annotation tables for REDItools. As a specific alternative, GFFtoTabix.py may be used. Option -b should be tuned according to
memory capabilities in case of sorting. Tabix by default compresses the input table and then creates indices. If
a copy of the sorted and uncompressed table is needed, -u option should be used.

Options:
|   option   |   description  |
| :--------: | :------------- |
|  -i	|	TAB-delimited file
|  -s	|	Sequence name column [1]
|  -c	|	Start column [4]
|  -e	|	End column (can be identical to -c) [5]
|  -m	|	Skip lines starting with [#]
|  -0	|	Zero based coordinates
|  -S	|	Do not sort input file (sort by default)
|  -b	|	Buffer size (as number of lines) [32000]
|  -t	|	Temporary directory to use (multiple -t may be used)
|  -u	|	Save an uncompressed GFF copy (add _copy suffix)
|  -h	|	Print this help

Example for GFF file:
```bash
 tableToTabix.py -i mytable.gff -u
```
### SortGFF.py (new in version 1.0.3)
SortGFF.py is a facility to sort only GFF files and works as SortTable.py. Option -b should be tuned
according to memory capabilities. Useful for users working on machines in which the sort program is not present.

Options:
|   option   |   description  |
| :--------: | :------------- |
|   -i	|	GFF file
|   -o	|	Sorted output file [GFF_sorted_%s]
|   -b	|	Buffer size (as number of lines) [32000]
|   -t	|	Temporary directory to use (multiple -t may be used)
|   -h	|	Print this help

Example:
```bash
 SortGFF.py -i mytable.gff -u
```
### GFFtoTabix.py (new in version 1.0.3)
GFFtoTabix.py is a script to convert GFF files to tabix files. Useful to generate GFF annotations for REDItools.
It works as tableToTabix.py but specific for GFF files. Use -u to store an uncompressed version of input GFF.
Sorting is activated by default.

Options:
|   option   |   description  |
| :--------: | :------------- |
|  -i	|	GFF file
|  -S	|	Do not sort GFF (sort by default)
|  -b	|	Buffer size (as number of lines) [32000]
|  -t	|	Temporary directory to use (multiple -t may be used)
|  -u	|	Save an uncompressed GFF copy (add _copy suffix)
|  -h	|	Print this help

Example:
```bash
 GFFtoTabix.py -i mytable.gff -u
```
### TableToGFF.py (new in version 1.0.3)
TableToGFF.py is a facility to convert an output REDItool table to GFF in order to be used as input in REDItoolDnaRna.py.
Optionally, it can sort and tabix the output GFF. Options -b and -T work if -s is in effect.

Options:
|   option   |   description  |
| :--------: | :------------- |
|  -i	|	Table file from REDItools
|  -s	|	Sort output GFF
|  -t	|	Tabix output GFF (requires Pysam module)
|  -b	|	Buffer size (as number of lines) [32000] (requires -s)
|  -T	|	Temporary directory (requires -s)
|  -o	|	Outfile [outTable_%s.gff]
|  -h	|	Print this help

Example:
```bash
 TableToGFF.py -i mytable -s -t -o myoutputTable.gff
```
## Usage Example
### Test dataset
The following example shows how to use REDItools to detect RNA editing sites in Alu regions. To this aim, please
download and extract testREDItools.tar.gz:
```bash
 gunzip testREDItools.tar.gz
 tar -xvf testREDItools.tar
 cd testREDItools
 ls
 dna.bam      reditool-output-sample.tar.gz  reference.fa.fai       refGene.sorted.gtf.gz.tbi  rmsk.gtf.gz.tbi  rna.bam.bai
 dna.bam.bai  reference.fa            refGene.sorted.gtf.gz  rmsk.gtf.gz                rna.bam
```
The testREDItools folder contains:
 * rna.bam: aligned RNA-Seq reads in BAM format (coordinate sorted by samtools sort)
 * rna.bam.bai: index file for rna.bam (indexed by samtools index)
 * dna.bam: aligned DNA-Seq reads (coordinate sorted by samtools sort)
 * dna.bam.bai: index file for dna.bam (indexed by samtools index)
 * reference.fa: reference genome in fasta format
 * reference.fa.fai: index file for reference.fa (indexed by samtools faidx)
 * rmsk.gtf.gz: gzipped gtf file of repeated elements (sorted by sort -k1,1 -k2,4n)
 * rmsk.gtf.gz.tbi: tabix index file for rmsk.gtf.gz
 * refGene.sorted.gtf.gz: gzipped gtf file of RefSeq genes (sorted by sort -k1,1 -k2,4n)
 * refGene.sorted.gtf.gz.tbi: tabix index file for refGene.sorted.gtf.gz
 * reditool-output-sample.tar.gz: folder including output examples

RNA and DNA reads were extracted from the region chr21:47721047-47743785.
Reads were obtained according to Ramaswami et al. paper. RNA reads mapped by BWA where kindly provided by
Gokul Ramaswami.
If REDItools have been correctly installed in your path, you can call all REDItoolDnaRNA.py options:
```bash
 REDItoolDnaRna.py -h
```
If everything is ok, REDItoolDnaRna.py can be launched by:
```bash
 [epicardi@srv00 sample]$ REDItoolDnaRna.py -i rna.bam -j dna.bam -f reference.fa -o reditool-test -c 10,1 -Q 33,64 -q 25,25 -m 20,20 -s 2 -g 1 -u -a 6-0 -v 2 -n0.0 -N0.0 -V
 Script time --> START: 06/02/2013 20:44:48
 Analysis ID: 891206177
 Analysis on 1 regions.
 Started analysis on region: chr21
 Job completed for region: chr21
 Merging Tables.
 Results saved on reditool-test/DnaRna_891206177/outTable_891206177
 Script time --> END: 06/02/2013 20:45:13
```
In this case we are requiring to extract RNA and DNA positions with a minimal coverage of 10 for DNA and 1 for RNA, a minimum quality score of 25 for both,
a minimum mapping quality of 20 for both. Please consider that quality scores in DNA are in Sanger format, so 33 is used as offset to calculate the phred score
per site genomic. On the contrary, quality scores for RNA reads are in the Illumina 1.3+ format and, thus, the value 64 has to be used through the option -Q.
Since RNA reads are strand oriented and the second read of the pair maintains the RNA orientation, we ask REDItoolDnaRna.py to infer the strand by -s 2 (second in pair good for orientation),
and assign to each position the overrepresented strand (-g 1 option). In addition, we remove first 6 nucleotides from each read (-a 6-0) and require sites supported by at least 2
variant bases (-v 2) without taking into account the frequency of variation in both RNA (-n 0.0) and DNA (-N 0.0). Finally, we exclude positions not supported by DNA-Seq reads (-V).
At the end of the run, REDItoolDnaRna.py will create the folder reditool-test/DnaRna_891206177/ and the output table outTable_891206177. The folder will contain also the parameters.txt
file summarizing all parameter values used for the run.
Now enter into the output folder:
```bash
 cd reditool-test/DnaRna_891206177/
 ls
 outTable_891206177  parameters.txt
```
and run the accessory script selectPositions.py to filter out other positions:
```bash
 selectPositions.py -i outTable_891206177 -d 12 -c 2 -C 10 -v 2 -V 0 -f 0.1 -F 1.0 -e -u -o candidates.txt
 Script time --> START: 06/02/2013 20:46:20
 Reading table...
 Total lines: 17879
 Filtered in lines: 10
 Selected lines saved on candidates.txt
 Script time --> END: 06/02/2013 20:46:21
```
In this case we are excluding positions showing an editing frequency lower that 0.1 and supported by heterozigous DNA sites.
Option -u is used to collect only positions supported by DNA-Seq. Option -e is used to remove RNA positions showing multiple substitutions.
According to these criteria only 10 positions will be retained:

 more candidates.txt
 Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency
 chr21   47739131        A       0       20      34.60   [18, 0, 2, 0]   AG      0.10    30      30.40   [30, 0, 0, 0]   -       0.00
 chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00
 chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00
 chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00
 chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00
 chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00
 chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00
 chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00
 chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00
 chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00

Now we can annotate these positions using the information in the Repeat Mask database to look at Alu sites:
```bash
 AnnotateTable.py -a ../../rmsk.gtf.gz -i candidates.txt -u -c 1,2,3 -n RepMask -o candidates.rmsk.txt
 Script time --> START: 06/02/2013 20:48:29
 Table saved on candidates.rmsk.txt
 Script time --> END: 06/02/2013 20:48:29
```
 more candidates.rmsk.txt
 Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency      RepMask_feat    RepMask_gid     RepMask_tid
 chr21   47739131        A       0       20      34.60   [18, 0, 2, 0]   AG      0.10    30      30.40   [30, 0, 0, 0]   -       0.00    -       -       -
 chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
 chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00    SINE    AluSz   Alu-SINE
 chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE
 chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE

To remove positions not annotated in SINE regions we use the FilterTable.py script:
```bash
 FilterTable.py -i candidates.rmsk.txt -f ../../rmsk.gtf.gz -F SINE -E -o candidates.rmsk.alu.txt -p
 Script time --> START: 06/02/2013 20:50:25
 Reading Table file...
 All positions: 10
 Positions filtered in: 9
 Positions filtered out: 1
 Script time --> END: 06/02/2013 20:50:25
```
Finally we can add gene annotations using RefSeq database:
```bash
 AnnotateTable.py -a /home/epicardi/annotation/hg19/annotation/tabixAnn/refseq/refGene.sorted.gtf.gz -i candidates.rmsk.alu.txt -u -c 1,2 -n RefSeq -o candidates.rmsk.alu.ann.txt
 Script time --> START: 06/02/2013 20:52:48
 Table saved on candidates.rmsk.alu.ann.txt
 Script time --> END: 06/02/2013 20:52:48
```
 more candidates.rmsk.alu.ann.txt
 Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency      RepMask_feat    RepMask_gid     RepMask_tid     RefSeq_feat     RefSeq_gid
 chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
 chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00    SINE    AluSz   Alu-SINE     intron  C21orf58
 chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE     intron  C21orf58
 chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE     intron  C21orf58

## Contact
* Ernesto Picardi: ernesto.picardi at gmail dot com
