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

python bin/fibronectin-splice-var-determiner -r hg38 -s -n tests/data/test_001.bam

