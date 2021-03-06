# REDItools: python scripts for RNA editing detection by RNA-Seq data

## Introduction
REDItools are python scripts developed with the aim to study RNA editing at genomic scale
by next generation sequencing data. RNA editing is a post-transcriptional phenomenon
involving the insertion/deletion or substitution of specific bases in precise RNA localizations.
In human, RNA editing occurs by deamination of cytosine to uridine (C-to-U) or mostly by the
adenosine to inosine (A-to-I) conversion through ADAR enzymes. A-to-I substitutions may have
profound functional consequences and have been linked to a variety of human diseases including
neurological and neurodegenerative disorders or cancer. Next generation sequencing technologies
offer the unique opportunity to investigate in depth RNA editing even though no dedicated
software has been released up to now.

REDItools are simple python scripts conceived to facilitate the investigation of RNA editing
at large-scale and devoted to research groups that would to explore such phenomenon in own
data but don’t have sufficient bioinformatics skills.
They work on main operating systems (although unix/linux-based OS are preferred), can handle reads from whatever
platform in the standard BAM format and implement a variety of filters.

## About this repo
Updated  REDItools to use a more modern version of [pysam](https://github.com/pysam-developers/pysam) (v0.15.0 is tested).

## Installation
```bash
git clone git://github.com/YongLiedu/reditools-updated
cd reditools-updated
python setup.py build
sudo python setup.py install
```

## Usage
See [documentation](https://github.com/YongLiedu/reditools-updated/blob/master/REDItools_documentation.md) or visit [here](http://srv00.recas.ba.infn.it/reditools/) for details.
