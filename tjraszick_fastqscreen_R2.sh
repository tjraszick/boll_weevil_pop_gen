#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J tjraszick_BW_fastqscreen            # job name
#BSUB -n 20                     # assigns 20 core for execution
#BSUB -R "span[ptile=20]"       # assigns 20 core per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2700MB process enforceable memory limit. (M * n)
#BSUB -W 72:00                  # sets to 92 hour the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid


module load FastQScreen/0.12.1-GCCcore-6.3.0-Perl-5.24.0

<<README
    - FastQ Screen homepage: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
README

################################################################################
# TODO Edit these variables as needed:
threads=20                       # make sure this is <= your BSUB -n value

aligner='bwa'                   # bwa, bowtie, bowite2
out_directory="filter"

################################################################################
# full path to aligners now required; uncomment DBs below as needed
# make sure the --filter option matches the number of DATABASEs below
# if you have two DATABASE entries then the value for --filter must have 2 characters like --filter 00
echo "
BWA $EBROOTBWA/bin/bwa

DATABASE mito_genome  ./mito_genome/Circularized_assembly_1_AgrandisMT.fasta

" > dbs.conf

tamulauncher tjraszick_fastqscreen_R2_commands.txt
