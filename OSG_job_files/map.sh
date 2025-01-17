#!/bin/bash

# rename input variables
trimmed_reads=${1}
ref_genome=${2}

# derive output file name
output_name="${1%_trimmed.fastq.gz}_out"

# might need to change the breseq version
tar -xzf breseq-0.37.1-Linux-x86_64.tar.gz

PATH=$PATH:$PWD/breseq-0.37.1-Linux-x86_64/bin/

# run breseq
breseq -r $ref_genome $trimmed_reads -o $output_name -j 8

# tar and gzip output directory and remove unecessary files
tar -cvzf ${output_name}.tar.gz $output_name
rm -r $output_name
rm -r $trimmed_reads