#!/bin/bash

# the first and only argument passed to this script is the trimmed reads file
# give the output file a clean name
trimmed_reads=${1}
output_name="${trimmed_reads%_trimmed.fastq.gz}_assembly"

# unzip and untar SPAdes
gzip -d SPAdes-3.15.4-Linux.tar.gz
tar -xf SPAdes-3.15.4-Linux.tar
PATH=$PATH:$PWD/SPAdes-3.15.4-Linux/bin/

# may need to modify the output file names here (and below) according to your input files
spades.py -s $trimmed_reads -o $output_name

tar -cvzf ${output_name}.tar.gz $output_name

rm -r $output_name
rm $trimmed_reads
rm -r SPAdes-3.15.4-Linux
rm SPAdes-3.15.4-Linux.tar
