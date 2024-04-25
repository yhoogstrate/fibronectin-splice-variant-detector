# fibronectin-splice-var-determiner #

This application investigates/counts 2 splice isoforms in the gene FN1
in RNA-seq data: 7B89 and 11A12.

Using RNA-seq, variants can be determined based on reads that:

 - (1) fall exactly over the splice junction (use `-s` / `--spliced-only`)
 - (2) are 'spanning' the splice junction and are thus mapped perfectly within exons 1 and 8

**The only thing that you need are properly aligned RNA-seq BAM files and this tool.**

This package is derived from the EGFRvIII determiner application:

[10.1093/neuonc/noab231](https://doi.org/10.1093/neuonc/noab231)


## How to install and test ##

As provided in `./quick_test.sh`:

```
#!/bin/bash

if [ ! -d ".venv" ] 
then
    virtualenv -p python3 .venv
fi

if [ ! -r "tmp/test_001.bam.bai" ] || [ ! -r "tmp/test_001.bam" ] 
then
    samtools view -bS tests/data/test_001.sam > tmp/test_001.bam
    samtools index tmp/test_001.bam
fi

source .venv/bin/activate
pip install -e .


mkdir -p tmp/

python bin/fibronectin-splice-var-determiner -r hg38 -s -n tmp/test_001.bam > tmp/test_001_sj_01.txt
python bin/fibronectin-splice-var-determiner -r hg38 -s    tmp/test_001.bam > tmp/test_001_sj_02.txt

python bin/fibronectin-splice-var-determiner -r hg38 -n    tmp/test_001.bam > tmp/test_001_cov_01.txt
python bin/fibronectin-splice-var-determiner -r hg38       tmp/test_001.bam > tmp/test_001_cov_02.txt

```


## What is does ##

We designed a small python tool for estimnating the read counts and/or extracting
the actual read names that allows to further analysed the sequencing data.


```
$ fibronectin-splice-var-determiner -r hg38 -s tmp/test_001.bam
```

Will result in a text file like this:

| sample | 7B89: wt reads [24 -> 26] | 7B89: sv reads [24 -> 25] | 7B89: sv reads [25 -> 26] | 7B89: sv reads [24 -> 25 -> 26] | 7B89: discrepant reads | 11A12: wt reads [31 -> 33] | 11A12: sv reads [31 -> 32] | 11A12: sv reads [31 -> 32] | 11A12: sv reads [31 -> 32 -> 33] | 11A12: discrepant reads |
|--------|---------------------------|---------------------------|---------------------------|---------------------------------|------------------------|----------------------------|----------------------------|----------------------------|----------------------------------|-------------------------|
| tmp/test_001.bam | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |



### Genomic reference ###

The genomic reference (hg19/hg38) can be changed using the `-r` or
`--reference-build` argument. Genomic references starting with '>1' rather
than 'chr1' will be automatically resolved.


## LICENSE ##

This is FREE software without liability or warranty, following the GPL-3
software license.


