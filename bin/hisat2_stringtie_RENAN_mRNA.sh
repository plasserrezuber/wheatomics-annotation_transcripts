#!/bin/bash
#SBATCH --job-name=hst2Re
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -p gdec
#SBATCH --cpus-per-task=16
#SBATCH --array=0-27

### modules management 
module load HISAT2/2.0.5 gcc/8.1.0 samtools/1.9 stringtie/2.0.3

### errors management
### set -x #mode debug on
set -o errexit #pour que le script s'arrete en cas d erreur ignoree
set -o nounset #force l initialisation des variables
IFS=$'\n\t'

##################################################################################
#set up directories
DATA_DIR='/storage/groups/gdec/shared/triticum_aestivum/wheatomics/rnaseq/renan'
DATABANK='/storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/assembly/pseudo_v2/TaeRenan_refseq_v2.0.hisat2idx'
OUTPUT='/home/palasser/results/hisat2/TaeRenan_v2.0_mRNA'

### working directory analysis in scratch
SCRATCHDIR=/storage/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p -m 700 ${SCRATCHDIR}
cd ${SCRATCHDIR}
mkdir -p ${OUTPUT}

##################################################################################
### Folder name contening the Renan mRNA .gz files : 28 samples : 16 grain / 2 root / 6 leaves / 4 steam
SAMPLES=(AACL AACM AACN AACO AACP AACQ AACR AACS AACT AACU AACV AACW AADE AADF AADU AADV AAHD AAHE AAHN AAHO AAIH AAII AAJB AAJC AAIR AAIS AAHX AAHY)

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
FILE=($(\ls -1 $DATA_DIR"/"$sample"/RunsSolexa/"*_JARVIS_*/*".fastq.gz"))
SHORTNAME=$(basename ${FILE[0]} |cut -d'_' -f1-3,5)
########################################################################################################

echo '*****************************************************************************'
echo "start job  ${SLURM_JOBID} on ${HOSTNAME} at `date`"
echo '***************** hisat2 alignment of sample' ${SHORTNAME}' *****************'
#####################################################################
##### hisat2 alignment of raw reads against genome ##################
#####################################################################
hisat2 -p 16 -x ${DATABANK} -1 ${FILE[0]} -2 ${FILE[1]} |samtools sort -@ 16 --output-fmt BAM -o ${SHORTNAME}_sorted.bam -

######### Make bam file index #######################################
samtools index -c -@ ${SLURM_CPUS_PER_TASK} ${SHORTNAME}_sorted.bam

#######################################
##### Filter bam file   ###############
#######################################
samtools view -@ ${SLURM_CPUS_PER_TASK} -F2308 -b -q60 ${SHORTNAME}_sorted.bam > ${SHORTNAME}_sorted_q60.bam

#############################################################
##### samtool flagstat for RAWreads #########################
#############################################################
samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${SHORTNAME}_sorted.bam > ${SHORTNAME}_sorted_stats.txt

###########################################
##### samtool Flagstat mapQ 60  ###########
###########################################
samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${SHORTNAME}_sorted_q60.bam > ${SHORTNAME}_sorted_q60_stats.txt

echo 'hisat2 alignment - job finished at: '`date`

#######################################
#####       StringTie   ###############
#######################################
echo '***************** stringtie treatment of sample ' ${SHORTNAME}' *****************'
##option -p (threading) fait que le gtf n'est pas trie
##option -c 4 coverage 4 min per base pour considerer la presence d'un transcrit
stringtie -p ${SLURM_CPUS_PER_TASK} -c 4 ${SHORTNAME}_sorted_q60.bam -o ${SHORTNAME}_q60.gtf

echo 'stringTie assembly - Job finished at: '`date`
##############################################################################################################################################
### move output data
mv ${SCRATCHDIR}/${SHORTNAME}* ${OUTPUT}/ && rm -Rf ${SCRATCHDIR}


#######################################
##### StringTie merge on all GTF ######
#######################################
##sbatch -p fast -c 32 --mem=32G --wrap="ml stringtie/2.0.3; stringtie --merge -p 32 -c 2 -o /home/palasser/results/hisat2/TaeRenan_v2.0_mRNA/TaeRenan_v2.0_wheatomics_28mRNA_stringtie.gtf /home/palasser/results/hisat2/TaeRenan_v2.0_mRNA/CGC_*_q60.gtf"
