#!/bin/sh

##  software: wigToBigWig
## version: 
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
chmod +x wigToBigWig

## software: bigWigAverageOverBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed

## software: fetchChromSizes
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x fetchChromSizes

## reference genome (b37):
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai -O hs37d5.fa.fai

##repeat regions:
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

export PATH=$PATH:$PWD

wget https://bismap.hoffmanlab.org/raw/hg19.umap.tar.gz  
tar -zxvf hg19.umap.tar.gz  
cd hg19  
fetchChromSizes hg19 > hg19.chrom.sizes  
wigToBigWig <(zcat k24.umap.wg.gz) hg19.chrom.sizes k24.umap.wg.bw  
