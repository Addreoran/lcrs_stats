#!/bin/bash
mkdir data
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cp uniprot_sprot.fasta.gz  ./data/swiss.fasta.gz
rm uniprot_sprot.fasta.gz 
gunzip ./data/swiss.fasta.gz
python3 get_lcrs_sequences.py
python3 likelihood.py
