#!/bin/bash


if [ ! -d ".venv" ] 
then
    virtualenv -p python3 .venv
fi

if [ ! -r "tests/data/test_001.bam.bai" ] 
then
    samtools view -bS tests/data/test_001.sam > tests/data/test_001.bam
    samtools index tests/data/test_001.bam
fi

source .venv/bin/activate
pip install -e .


mkdir -p tmp/

python bin/fibronectin-splice-var-determiner -r hg38 -s -n tests/data/test_001.bam > tmp/test_001_sj_01.txt
python bin/fibronectin-splice-var-determiner -r hg38 -s    tests/data/test_001.bam > tmp/test_001_sj_02.txt

python bin/fibronectin-splice-var-determiner -r hg38 -n    tests/data/test_001.bam > tmp/test_001_cov_01.txt
python bin/fibronectin-splice-var-determiner -r hg38       tests/data/test_001.bam > tmp/test_001_cov_02.txt

