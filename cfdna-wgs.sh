#! /bin/bash

module load snakemake
module load singularity
cd /data/panalc/cfdna-wgs

# Necessary to run conda snakemake command in shell script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

##################################################################
###                          TESTING                           ###
##################################################################
## dry-run 
# snakemake -s workflow/int_test.smk --configfile config/int_test_biowulf.yaml --profile profile/ -np
## actual run
# snakemake -s workflow/int_test.smk --configfile config/int_test_biowulf.yaml --profile profile/

##################################################################
###                        MPNST dataset                       ###
##################################################################
## dry-run 
# snakemake -s workflow/mpnst.smk --configfile config/mpnst.yaml --profile profile/ -np
## update timestamps
# snakemake -s workflow/mpnst.smk --configfile config/mpnst.yaml --profile profile/ --touch
## actual run
snakemake -s workflow/mpnst.smk --configfile config/mpnst.yaml --profile profile/

## START SCRIPT:
# sbatch --mem=1g --cpus-per-task=2 --time=24:00:00 --mail-type=BEGIN,TIME_LIMIT_90,END /data/panalc/cfdna-wgs/cfdna-wgs.sh
