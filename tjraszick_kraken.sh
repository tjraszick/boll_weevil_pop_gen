#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J tjraszick_kraken_pe    # job name
#BSUB -n 4                      # assigns 4 cores for execution
#BSUB -R "span[ptile=4]"        # assigns 4 cores per node
#BSUB -R "rusage[mem=248750]"   # reserves 248750MB memory per core
#BSUB -M 248750                 # sets to 248750MB per process enforceable memory limit. (M * n)
#BSUB -q xlarge
#BSUB -W 96:00                  # sets to 96 hour the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load Westmere
module load Kraken/1.1-GCCcore-6.3.0-Perl-5.24.0

<<README
    - Kraken manual: https://ccb.jhu.edu/software/kraken/MANUAL.html
README

################################################################################
# TODO Edit these variables as needed:
export KRAKEN_DEFAULT_DB='bacteria'
export KRAKEN_DB_PATH='/scratch/datasets/kraken'
export KRAKEN_NUM_THREADS=10

################################################################################
# you only need to run --preload for the first command if you have multiple samples
kraken --fastq-input --fastq-output --gzip-compressed --unclassified-out unc_Chihuahua2014_14AA --classified-out class_Chihuahua2014_14AA --out-fmt paired --output Chihuahua2014_14AA_kraken.out --preload --paired Chihuahua2014_14AA_R1.fastq.gz Chihuahua2014_14AA_R2.fastq.gz
tamulauncher tjraszick_kraken_commands.txt

# annotate classified hits, currently comment disabled
# kraken-translate kraken.out > sequences.labels

# show report of percentage of hits to database entries, currently comment disabled
# kraken-report kraken.out > report.out

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/wiki/index.php/HPRC:AckUs
    - Kraken:
        Wood DE, Salzberg SL: Kraken: ultrafast metagenomic sequence classification
        using exact alignments. Genome Biology 2014, 15:R46.
CITATION