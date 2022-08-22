#!/bin/bash
#SBATCH --job-name=mergeST
#SBATCH -o slurmjob-%A
#SBATCH --mem=64G
#SBATCH -p fast #### normal debug fast long smp gdec
#SBATCH --cpus-per-task=16 #nb thread
#SBATCH --mail-user=philippe.leroy.2@inrae.fr
#SBATCH --mail-type=ALL 
__author__='Emilie Gerard-Marchant'
__email__='emilie.gerard-marchant@inrae.fr'
##################################################################################
#### Modify Wed Apr 28 16:34:38 CEST 2021 - P. Leroy
#### Last modification Tue May 18 14:14:52 CEST 2021 - P. Leroy
##################################################################################
if [ "$#" -ne 2 ]; then
	echo
    echo "***** Illegal number of parameters for 03_mergeStringTie.sh *****"
    echo "Needs: PathFolderGTFfileList; PathFolderOutput "
	echo
	echo "usage: sbatch 03_mergeStringTie.sh PathFolderGTFfileList PathFolderOutput"
	echo "      1. PathFolderGTFfileList: complete path of the file (.txt) which contains the full path of each mRNAseq sample stringtie gtf file"
	echo "      2. PathFolderOutput: complete path for final stringtie merge gtf output file"
else
	echo
	echo "usage: sbatch 03_mergeStringTie.sh PathFolderGTFfileList PathFolderOutput"
	echo "      1. PathFolderGTFfileList: complete path of the file (.txt) which contains the full path of each mRNAseq sample stringtie gtf file"
	echo "      2. PathFolderOutput: complete path for stringtie output result files"
fi
GTF=$1
OUTPUT=$2
GTFlast="Renan_mRNAseqHisatStringTieMerge.gtf"
echo
echo '***********************************'
echo 'Transcripts merge with StringTie   '
echo '***********************************'
echo 'File with gtf stringtie path list                     : '${GTF}
echo 'Folder for latest gtf stringtie output result files   : '${OUTPUT}
#set -x #mode debug on
set -o errexit 
set -o nounset 
IFS=$'\n\t'
### modules
module purge
module load stringtie/2.0.3
echo "start job ${SLURM_JOB_NAME} with jobid ${SLURM_JOBID} on ${HOSTNAME} at `date`"
echo 'Set up directories ...'
SCRATCHDIR=/storage/scratch/"$USER"/"$SLURM_JOB_ID"
mkdir -p -m 700 ${SCRATCHDIR}
mkdir -p ${OUTPUT}
cd ${SCRATCHDIR}
echo
echo '**********************'
echo 'stringtie merge...    '
echo '**********************'
echo
#################################################################################
	CMD="stringtie --merge -p 16 -o ${GTFlast} ${GTF}"
	echo $CMD
	stringtie --merge -p 16 -o ${GTFlast} ${GTF}
		#-p : nb de threads pour l'assemblage
		#-o : nom fichier GTF d'output
##################################################################################
mv  ${SCRATCHDIR} ${OUTPUT}
module purge
echo
echo '03_mergeStringTie.sh - Job finished at: '`date`
##### end of bash script