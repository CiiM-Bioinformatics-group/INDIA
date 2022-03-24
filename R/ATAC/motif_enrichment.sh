#!/bin/bash

#$ -N findMotif
#$ -l arch=linux-x64
#$ -b n
#$ -i /dev/null
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -q all.q
export PATH=/vol/projects/CIIM/resources/tools/homer/bin/:$PATH

findMotifsGenome.pl /vol/projects/CIIM/INDIA/analysis/ATAC/output/DE_peaks.bed hg38 homer_enr/ -bg /vol/projects/CIIM/INDIA/analysis/ATAC/output/background_peaks.bed -size 200
