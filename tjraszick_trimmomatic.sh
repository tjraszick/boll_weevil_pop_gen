#BSUB -L /bin/bash
#BSUB -J tjraszick_trimmomatic
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=12300]"    # reserves 12300MB memory per core
#BSUB -M 12300                  # sets to 12300MB per process enforceable memory limit. (M * n)
#BSUB -W 72:00

module load Trimmomatic/0.38-Java-1.8.0

tamulauncher tjraszick_trimmomatic_commands.txt
