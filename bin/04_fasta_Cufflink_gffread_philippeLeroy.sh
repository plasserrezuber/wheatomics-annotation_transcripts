#!/bin/bash
#SBATCH --job-name=gffread
#SBATCH -o slurmjob-%A
#SBATCH --mem=4G # 64
#SBATCH -p fast #### normal debug fast long smp gdec
#SBATCH --cpus-per-task=1 #nb thread
#SBATCH --mail-user=philippe.leroy.2@inrae.fr
#SBATCH --mail-type=ALL 
__author__='Emilie Gerard-Marchant'
__email__='emilie.gerard-marchant@inrae.fr'
##################################################################################
#### Modify Mon May  3 10:49:09 CEST 2021 - P. Leroy
#### Last modification Tue May 25 23:20:11 CEST 2021 - P. Leroy
##################################################################################
if [ "$#" -ne 4 ]; then
	echo
    echo "***** Illegal number of parameters for gffreadFasta.sh *****"
    echo "Needs: PathGTFfile; PathFastaRefSeqGenome; PathFolderOutput; FastaName "
	echo
	echo "usage: sbatch gffreadFasta.sh PathGTFfile PathFastaRefSeqGenome"
	echo "      1. PathGTFfile: complete path of the stringtie gtf list file"
	echo "      2. PathFastaRefSeqGenome: complete path for reference genome fasta file (folder should contains the .fasta.fai index)"
	echo "      3. PathFolderOutput: complete path for final fasta output file"
	echo "      4. FastaName: Name of the fasta file"
else
	echo
	echo "usage: sbatch agffreadFasta.sh PathGTFfile PathFastaRefSeqGenome"
	echo "      1. PathGTFfile: complete path of the stringtie gtf list file"
	echo "      2. PathFastaRefSeqGenome: complete path for reference genome fasta file (folder should contains the .fasta.fai index)"
	echo "      3. PathFolderOutput: complete path for final fasta output file"
	echo "      4. FastaName: Name of the fasta file with .fasta extension"
fi
GTF=$1
REF=$2
OUTPUT=$3
FASTA=$4
echo
echo '**************************************'
echo 'Transcripts extraction with gffread   '
echo '**************************************'
echo 'Stringtie final gtf file                      : '${GTF}
echo 'Fasta file of the reference genome            : '${REF}
echo 'Folder for the transcripts fasta result files : '${OUTPUT}
echo 'Name of the fasta file                        : '${FASTA}
#set -x #mode debug on
set -o errexit 
set -o nounset 
IFS=$'\n\t'
### modules
module purge
module load cufflinks/2.2.1 gcc/4.8.4 samtools/1.3
echo "Set up directories ..."
SCRATCHDIR=/storage/scratch/"$USER"/"$SLURM_JOB_ID"
mkdir -p -m 700 ${SCRATCHDIR} 
mkdir -p ${OUTPUT}
cd ${SCRATCHDIR}
echo
echo '*************************'
echo 'cufflinks gffread ...    '
echo '*************************'
echo
echo "start job ${SLURM_JOB_NAME} with jobid ${SLURM_JOBID} on ${HOSTNAME} at `date`"
echo
#################################################################################
	CMD="gffread ${GTF} -w ${FASTA} -g ${REF}"
	echo $CMD
	echo '****************  make fai index with samtools ****************'
	samtools faidx ${REF}
	echo '**************** make gffreads ****************'
	gffread ${GTF} -w ${FASTA} -g ${REF}
		#-w : output = transcripts fasta file
		#-g : input = reference genome
##################################################################################
mv  ${SCRATCHDIR} ${OUTPUT}
module purge
echo
echo '04_fasta_Cufflink_gffread.sh - Finished job : '`date`
##### end of bash script
