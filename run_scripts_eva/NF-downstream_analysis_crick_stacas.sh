#!/bin/bash
#SBATCH --job-name=10x-NPB
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=eva.hamrud@crick.ac.uk

export TERM=xterm
export NXF_VER=20.07.1
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/working/hamrude/NF_singularity
export NXF_HOME=/flask/scratch/briscoej/hamrude/10x_neural_plate_border_EH/NF-downstream_analysis
export NXF_WORK=work/

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

nextflow run ./NF-downstream_analysis/main.nf \
--input ./NF-downstream_analysis/samplesheet.csv \
--outdir ./output/NF-downstream_analysis_stacas \
--debug \
--integration STACAS \
-profile crick_eva \
-resume