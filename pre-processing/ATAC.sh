#!/bin/bash

export NXF_OPTS='-Xms1g -Xmx4g'
export NEXTFLOW=/vol/projects/CIIM/resources/tools/nextflow

${NEXTFLOW} pull nf-core/atacseq

${NEXTFLOW} run nf-core/atacseq  -r 1.2.1 \
    -c '/home/mzoodsma/nextflow.config' \
    --input 'design.csv' \
    -profile 'singularity' \
    --genome 'GRCh38' \
    --igenomes_base '/vol/projects/CIIM/refs/igenomes' \
    -resume \
    --skip_merge_replicates
