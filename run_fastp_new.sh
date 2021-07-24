#!/bin/bash

#SBATCH --account p31288
#SBATCH --partition short
#SBATCH --nodes=1
#SABTCH --ntasks-per-node=12
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb
#SBATCH --job-name="run_fastp_new.sh"

# load modules
module purge all
module load anaconda3
module load fastqc
module load multiqc

# Set working directory
cd $SLURM_SUBMIT_DIR

# Timestamp logs
timestamp=$(date +%d%m%Y_%H%M)


###########################
# PART 1: Download + Path #
###########################

numberSRR="SRX5509390"
readsPath="/home/cmw3681/hpylori_scripts/hpylori_reads"
workspaceDir="${SLURM_SUBMIT_DIR}/${numberSRR}_${timestamp}"
if [ -d $workspaceDir ]
then
	echo "Directory ${workspaceDir} already exists"
else
	mkdir ${workspaceDir}
	echo "Directory ${workspaceDir} has been made"
fi

rawReads="${workspaceDir}/raw_reads"
trimDir="${workspaceDir}/trimmed_fastp"
fastqcDir="${trimDir}/fastqc"
metaphlanDir="${workspaceDir}/metaphlan"

mkdir ${trimDir} ${fastqcDir} ${metaphlanDir} ${rawReads}

############################
### PART 1: FASTP, QUAST ###
############################

# Bring raw reads into directory

cd ${readsPath}

# FASTP analysis
source /software/anaconda2/etc/profile.d/conda.sh
conda activate sequence_run

for i in $(ls *_1.fastq.gz) # Iterate through R1
do
	i_sub=$(echo $i | cut -c1-10) # SRR accession number
	for j in $(ls *_2.fastq.gz) # Iterate R2
	do
		j_sub=$(echo $j | cut -c1-10)
		if [[ $i_sub == $j_sub ]] # Match reads via accession
		then
			fastp -i ${i} -I ${j} --out1 ${i_sub}_fastp_out.R1.fq.gz --out2 ${j_sub}_fastp_out.R2.fq.gz --detect_adapter_for_pe --thread 16 --length_required 50

		fi

	done

done

conda deactivate

# To move FASTP outputs into new directory
mv *fastp_out* ${trimDir}/
mv fastp.html fastp.json ${trimDir}/

cd ${trimDir}/

# FASTQC analysis
fastqc -t 12 *fq.gz

# To move FastQC output into new directory
mv *fastqc.html ${fastqcDir}
mv *fastqc.zip ${fastqcDir}

# MultiQC analysis
cd ${fastqcDir}
multiqc .
cd $SLURM_SUBMIT_DIR


#######################
## PART 2: MetaPhlAn ##
#######################
