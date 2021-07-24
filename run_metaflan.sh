#!/bin/bash

#SBATCH --account p31288
#SBATCH --partition short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb
#SBATCH --job-name="run_metaflan"

# Load modules
module purge all
module load anaconda3
module load bowtie2

# Activate env with metaphlan
source /software/anaconda2/etc/profile.d/conda.sh
conda activate metaphlan-new-3.0
module load bowtie2

# Set working directory
cd $SLURM_SUBMIT_DIR

#########################
## METAPHLAN PROFILING ##
#########################

# Metaphlan profiling with trimmed reads
cd SRX5509390_21072021_0958/trimmed_fastp # Direct path, don't have workspaceDir variable 

for i in $(ls *_fastp_out.R1.fq.gz) # Iterate R1
do
	i_num=$(echo $i | cut -c1-10) # SRR accession number for R1
	for j in $(ls *_fastp_out.R2.fq.gz) # Iterate R2
	do
		j_num=$(echo $j | cut -c1-10)
		if [[ $i_num == $j_num ]] # Match paired reads
		then
			metaphlan --bowtie2db ~/dependecies/mpa_db ${i},${j} --bowtie2out ${i_num}.bowtie2.bz2 --input_type fastq > ${i_num}_profile.txt

		fi

	done

done

# Move metaphlan output into new directory
mv *profile.txt ../metaphlan
mv *bowtie2* ../metaphlan

cd ../metaphlan

# Generate single tab-delimited table from files
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt

