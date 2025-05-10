#!/bin/bash
#SBATCH --job-name=counting
#SBATCH --mem=32G
#SBATCH -c 1

### SCRIPT DETERMINES NUMBER OF READS MAPPING TO GENE FEATURES IN BAM FILE(S)
### SOFTWARE USED: HTseq

## JOB LOG HEADER
  perl -E 'say"="x80'; echo "JOB STARTED: $(date)"; echo "NODE: $(hostname)"; echo "SCRIPT ${0}:"; cat $0; perl -E 'say"="x80'

## PRIMARY VARIABLES - usually need to change
  INDIR="/home/assmith/project_1/outfiles/b-mapping/"                  # INPUT DIRECTORY
  OUTDIR="/home/assmith/project_1/outfiles/c-counting"                # OUTPUT DIRECTORY

## SECONDARY VARIABLES - rarely need to change                     # GTF FILE CONTAINING GENE ANNOTATIONS
  GTF="/work/binf-courses/binf694/project_1/reference/human-STAR/GRCh38.p14-ensembl/Homo_sapiens.GRCh38.113.gene_name.gtf"
  BAM_EXT="Aligned.sortedByCoord.out.bam"                          # BAM FILE EXTENSION - REMOVED TO NAME OUTPUT FILES

## COMMAND(s) TO RUN
 # TIP: SCRIPT SHOWS HOW TO USE A "FOR" LOOP TO RUN COMMANDS ON MULTIPLE FILES

  mkdir -p ${OUTDIR}
  cd ${OUTDIR}

  # Next line says FOR every file (BAM_FILE) in INDIR that ends in BAM_EXT, repeat commands between "do" and "done"
  for BAM_FILE in $( ls -1 ${INDIR}/*.${BAM_EXT} )
  do
    echo "PROCESSING FILE: ${BAM_FILE}"
    PREFIX=$(basename -s ".${BAM_EXT}" ${BAM_FILE})     # FINDS UNIQUE PART OF FILENAME (REMOVES PATH & EXTENSION)

    # RUN _HT-Seq_ TO COUNT THE NUMBER OF READS IN EACH GENE FEATURE
    htseq-count \
      --mode intersection-nonempty \
      --stranded no \
      --format bam \
      --type exon \
      --idattr gene_name \
      ${BAM_FILE} \
      ${GTF} \
      > ${OUTDIR}/${PREFIX}.count.txt
  done

#JOB LOG FOOTER
  perl -E 'say"="x80'; echo "JOB COMPLETED: $(date)"; perl -E 'say"="x80'
