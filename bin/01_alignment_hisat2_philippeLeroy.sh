#!/bin/bash
#SBATCH --job-name=HT2map
#SBATCH -o slurmjob-%A
#SBATCH --mem=512G
#SBATCH -p smp #### normal debug fast long *smp* gdec
#SBATCH --cpus-per-task=16 #nb thread
#SBATCH --mail-user=philippe.leroy.2@inrae.fr
#SBATCH --mail-type=ALL 
__author__='Emilie Gerard-Marchant'
__email__='emilie.gerard-marchant@inrae.fr'
##################################################################################
#### Refactoring Wed Apr 28 09:58:55 CEST 2021 - P. Leroy
#### Last modification Wed May 19 12:08:44 CEST 2021 - P. Leroy
##################################################################################
if [ "$#" -ne 3 ]; then
	echo
    echo "***** Illegal number of parameters for 01_alignment_hisat2.sh *****"
    echo "Needs: PathFolderforFastqFiles; PathForRefSeqGenomeIdx; PathFolderOutput "
	echo
	echo "usage: sbatch 01_alignment_hisat2.sh PathForRefSeqGenomeIdx PathFolderOutput"
	echo "      1. PathFolderforFastqFiles: complete path where are located the mRNASeq fastq files"
	echo "      2. PathForRefSeqGenomeIdx: complete path where is located the hisat2 reference genome index"
	echo "      3. PathFolderOutput: complete path for hisat2 bam output result files"
else
	echo
	echo "usage: sbatch 01_alignment_hisat2.sh PathForRefSeqGenomeIdx PathFolderOutput"
	echo "      1. PathFolderforFastqFiles: complete path where are located the mRNASeq fastq files"
	echo "      2. PathForRefSeqGenomeIdx: complete path where is located the hisat2 reference genome index"
	echo "      3. PathFolderOutput: complete path for hisat2 bam output result files"
fi
DATA_DIR=$1
DATABANK=$2
OUTPUT=$3
echo
echo '******************************************'
echo 'Mapping on reference genome with StringTie'
echo '******************************************'
echo 'Folder where are located mRNAseq fastq files    : '${DATA_DIR}
echo 'Folder where are located reference genome index : '${DATABANK}
echo 'Folder for hisat2 output result files           : '${OUTPUT}
### errors management
### set -x #mode debug on
set -o errexit #pour que le script s'arrete en cas d erreur ignoree
set -o nounset #force l initialisation des variables
IFS=$'\n\t'
### modules management 
### samtools needs gcc
module purge
module load HISAT2/2.0.5 gcc/8.1.0 samtools/1.9
echo 'Set up directories ...'
### working directory analysis in scratch
SCRATCHDIR=/storage/scratch/"$USER"/"$SLURM_JOB_ID"
mkdir -p ${OUTPUT}
mkdir -p -m 700 ${SCRATCHDIR} 
	### -p : make the poarent folders if necessary
	### -m : right acesses u+rwx
cd ${SCRATCHDIR}
### Folder name contening the Renan mRNA .gz files : 28 samples : 16 grain / 2 root / 6 leaves / 4 steam
SAMPLES=(AACL AACM AACN AACO AACP AACQ AACR AACS AACT AACU AACV AACW AADE AADF AADU AADV AAHD AAHE AAHN AAHO AAIH AAII AAJB AAJC AAIR AAIS AAHX AAHY)
echo "================================================================================================="
echo "fastq samples: "${SAMPLES}
echo "================================================================================================="
echo
echo '*******************'
echo 'HISAT2 alignment...'
echo '*******************'
echo "start job ${SLURM_JOB_NAME} with jobid ${SLURM_JOBID} on ${HOSTNAME} at `date`"
echo
echo "summary_colDesc: Sample_Name | total_raw_reads | Properly_paired | Mate_mapped_and_singletons | Percentage"
echo
########################################################################################################
for sample in ${SAMPLES[@]};
do
	FILE=($(ls "$DATA_DIR"/"$sample"/RunsSolexa/*_JARVIS_*/*.fastq.gz))
	SHORTNAME=$(echo $(basename $FILE .fastq.gz) | cut -d'_' -f1-3,5-)
	VERYSHORTNAME=$(echo $(basename $FILE .fastq.gz) | cut -d'_' -f2)
	echo '***************** treatment of ' ${VERYSHORTNAME}' *****************'
	#####################################################################
	##### hisat2 alignment of raw reads against genome ##################
	#####################################################################
	SAMFILE=${SHORTNAME}.mapped.sam
	hisat2 -p 16 -x ${DATABANK} -1 ${FILE[0]} -2 ${FILE[1]} -S ${SAMFILE}
	wait
	####################################################################
	##### Make bam file and remove sam file            ##################
	#####################################################################
	BAMFILE=${SHORTNAME}.mapped.sorted.bam
	samtools sort -@ 16 ${SAMFILE} | samtools view -@ 16 -bS -o ${BAMFILE}
	wait
	######### Make bam file index #######################################
	samtools index -c -@ 16 ${BAMFILE}
	rm -rf ${SAMFILE}
	wait
	#####################################################################
	##### recover samtool flagstat for raw reads ########################
	#####################################################################
	RAWSTATS=${SHORTNAME}.mapped.sorted.raw_stats
	samtools flagstat -@16 ${BAMFILE} > ${RAWSTATS}
	ReadTotal=`grep "paired in sequencing" ${RAWSTATS} | cut -d '+' -f 1`
	wait
	PROPERLYPAIRED=`grep "properly paired" ${RAWSTATS} | cut -d '+' -f 1`
	wait
	echo "***********************************************************"
	echo "***********************************************************"
	echo "Number of raw reads: "${ReadTotal}
	echo "Number of raw reads properly paired : "${PROPERLYPAIRED}
	#####################################################################
	##### Filter bam file for mate mapped and singletons  ###############
	#####################################################################
	FILTERSTATS=${SHORTNAME}.mapped.sorted.filter_stats
	samtools view -@16 -F2308 -bq20 ${BAMFILE} | samtools flagstat -@16 - > ${FILTERSTATS}
	wait
	MateMapped=`grep "with itself and mate mapped" ${FILTERSTATS} | cut -d '+' -f 1`
	wait
	echo "After filtering: number of reads mapped with itself and mate mapped: "${MateMapped}
	Singleton=`grep "singletons" ${FILTERSTATS} | cut -d '+' -f 1`
	wait
	SUM=$(($MateMapped+$Singleton))
	echo "After filtering: number of reads mapped as singletons: "${Singleton}
	echo "------------------------------------------------------------"
	printf "%s\t%.3f\n" "Pourcentage of mapped reads after filtering: " "$((10**3 * $SUM/$ReadTotal * 100))e-3"
	echo "------------------------------------------------------------"
	printf "%s\t%s\t%s\t%s\t%s\t%.3f\n" "summary:" "$VERYSHORTNAME" "$ReadTotal" "$PROPERLYPAIRED" "$SUM" "$((10**3 * $SUM/$ReadTotal * 100))e-3"
done
##############################################################################################################################################
### move output data
mv  ${SCRATCHDIR} ${OUTPUT}
echo
echo '01_alignment_hisat2.sh - job finished at: '`date`
########### end of bash script