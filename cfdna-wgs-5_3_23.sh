#! /bin/bash

module load snakemake
module load singularity
cd /data/panalc/cfdna-wgs

# Necessary to run conda snakemake command in shell script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

# NOTE: This is explicitly for the purpose of computing some missing data w/ an updated library sheet (cfdna_wgs-5_3_23.tsv)
# In particular, we were missing lib287 and lib297 from delfi data.
snakemake -s workflow/mpnst_5_3_23.smk --configfile config/mpnst.yaml --profile profile/

## START SCRIPT:
# sbatch --mem=1g --cpus-per-task=2 --time=12:00:00 --mail-type=BEGIN,TIME_LIMIT_90,END /data/panalc/cfdna-wgs/cfdna-wgs-5_3_23.sh
