#!/bin/bash
#SBATCH --job-name=stringTie
#SBATCH -o slurmjob-%A
#SBATCH --mem=64G 
#SBATCH -p normal #### normal debug fast long smp gdec
#SBATCH --cpus-per-task=16 #nb thread
#SBATCH --mail-user=philippe.leroy.2@inrae.fr
#SBATCH --mail-type=ALL 
__author__='Emilie Gerard-Marchant'
__email__='emilie.gerard-marchant@inrae.fr'
##################################################################################
#### Modify Wed Apr 28 12:12:56 CEST 2021 - P. Leroy
#### Last modification Tue May 18 14:14:52 CEST 2021 - P. Leroy
##################################################################################
if [ "$#" -ne 2 ]; then
	echo
    echo "***** Illegal number of parameters for 02_assembly_stringTie.sh *****"
    echo "Needs: PathFolderBamFile; PathFolderOutput "
	echo
	echo "usage: sbatch 02_assembly_stringTie.sh PathFolderBamFile PathFolderOutput"
	echo "      1. PathFolderBamFile: complete path where are located the aligned hisat2 bam files"
	echo "      2. PathFolderOutput: complete path for stringtie gtf output result files"
else
	echo
	echo "usage: sbatch 02_assembly_stringTie.sh PathFolderBamFile PathFolderOutput"
	echo "      1. PathFolderBamFile: complete path where are located the aligned hisat2 bam files"
	echo "      2. PathFolderOutput: complete path for stringtie gtf output result files"
fi
DATA_DIR=$1
OUTPUT=$2
echo
echo '***********************************'
echo 'Transcripts assembly with StringTie'
echo '***********************************'
echo 'Folder where are located aligned hisat2 bam files: '${DATA_DIR}
echo 'Folder for stringtie output result files         : '${OUTPUT}
set -x # print cmd before execution
set -o errexit 
set -o nounset 
IFS=$'\n\t'
### modules
module purge
module load stringtie/2.0.3
echo
echo 'Set up directories ...'
SCRATCHDIR=/storage/scratch/"$USER"/"$SLURM_JOB_ID"
mkdir -p -m 700 ${SCRATCHDIR} 
mkdir -p ${OUTPUT}
cd ${SCRATCHDIR} 
echo
echo '**********************'
echo 'stringtie assembly ...'
echo '**********************'
echo "start job ${SLURM_JOB_NAME} with jobid ${SLURM_JOBID} on ${HOSTNAME} at `date`"
echo
#################################################################################
FILE=($(ls ${DATA_DIR}/*.sorted.bam))
for BAM in "${FILE[@]}";
do
	BAMFILE=$(basename ${BAM})
	NAME=${BAMFILE%.*}
	EXT=${BAMFILE##*.}
	if [ ${EXT} = "bam" ]; then
		GTF=${NAME}.gtf
		echo '***************** treatment of ' ${GTF}' *****************'
		stringtie --fr -p 16 ${BAM} -o ${GTF}
			#--fr : presume librairy on forward strand 
			#-p : nb de threads pour l'assemblage
			#-o : GTF ouput folder files
	fi
done
##################################################################################
mv  ${SCRATCHDIR} ${OUTPUT}
module purge
echo
echo '02_assembly_stringTie.sh - Job finished at: '`date`
#### end of bash script