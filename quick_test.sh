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

