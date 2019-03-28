#!/bin/sh

##  software: wigToBigWig
## version: 
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
chmod +x wigToBigWig

## software: bigWigAverageOverBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed

## software: samtools
## VERSION: 1.9
## TYPE: file format converter
## SOURCE_URL: https://github.com/samtools/samtools
#wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
#tar -xjf samtools-1.9.tar.bz2
#cd samtools-1.9
#make
#cd ..
#ln -s samtools-1.9 samtools
