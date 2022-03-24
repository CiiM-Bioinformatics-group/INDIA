#!/bin/bash

export NXF_OPTS='-Xms1g -Xmx4g'
export NEXTFLOW=/vol/projects/CIIM/resources/tools/nextflow

${NEXTFLOW} pull nf-core/rnaseq

${NEXTFLOW} run nf-core/rnaseq -r 2.0 \
  -c '/home/mzoodsma/nextflow.config' \
  --input 'design.csv' \
  -profile 'singularity' \
  --genome 'GRCh38' \
  --igenomes_base '/vol/projects/CIIM/refs/igenomes' \
  --skip_bigwig \
  --skip_stringtie \
  -resume
