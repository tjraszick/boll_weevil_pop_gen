#BSUB -L /bin/bash
#BSUB -J Raszick_SchistoRAD_dDocent_raw
#BSUB -n 40
#BSUB -R "select[mem2tb]"
#BSUB -R "span[ptile=40]"
#BSUB -R "rusage[mem=49750]"
#BSUB -M 49750
#BSUB -q xlarge
#BSUB -W 144:00
#BSUB -o stdout.%J
#BSUB -e stderr.%J

module load Westmere

module load dDocent/2.6.0-intel-2017A-Java-1.8.0

dDocent tjraszick_dDocent_config
