#!/bin/sh

#André Soares, 24Nov2017
#trimming primer sequences from
#cutadapt is here: https://github.com/marcelm/cutadapt

for read1 in ls SGG_16S_data/*R1_001.fastq.gz;
  do
    #i.e. this for loop reads through all read 1 files and then uses that path to look for read 2 files
    # echo $read1
    read2=$(echo $read1| sed 's/R1_001.fastq/R2_001.fastq/')
    #modified $read1 and $read2 variables to change the name of the outputs to _trim.fastq
    read1_trim=$(echo $read1| sed 's/R1_001.fastq/R1_001_trim.fastq/')
    read2_trim=$(echo $read2| sed 's/R2_001.fastq/R2_001_trim.fastq/')
    # echo $read2
    # -a is for the F primer and -A is the R primer
    # reason they are repeated is because https://github.com/marcelm/cutadapt/issues/229
    #got to write >99.1% of reads on most samples with this, after trimming
    cutadapt -a CCTACGGGNGGCWGCAG...CCTACGGGNGGCWGCAG \
    -A GACTACHVGGGTATCTAATCC...GACTACHVGGGTATCTAATCC \
    $read1 $read2 -o $read1_trim -p $read2_trim

done

mv SGG_16S_data/*trim.fastq.gz SGG_16S_data_trimmed
