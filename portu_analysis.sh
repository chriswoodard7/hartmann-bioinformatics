#!/bin/bash

#SBATCH --account p31288
#SBATCH --partition normal
#SBATCH --nodes=1
#SABTCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10gb
#SBATCH --job-name="portu_analysis"

# This script downloads runs from the BioProject PRJNA532515, trims the paired reads, performs quality control,
# and profiles each read based on a MetaPhlAn database. Tables are then generated that can be used for data vis.

# load modules
module purge all
module load sratoolkit
module load anaconda3
module load fastqc
module load multiqc
module load bowtie2

# Working directory is location of submitted script

###########################
###### PART 1: Path #######
###########################

timestamp=$(date +%d%m%Y_%H%M)

numberSRR="PRJNA532515"
workspaceDir="/projects/p31288/chris/portu_lira/${numberSRR}_${timestamp}"
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

###########################
#### PART 2: Download #####
###########################

cd ${rawReads}

fasterq-dump --split-files SRR8944124
fasterq-dump --split-files SRR8944125
fasterq-dump --split-files SRR8944126
fasterq-dump --split-files SRR8944127
fasterq-dump --split-files SRR8944128
fasterq-dump --split-files SRR8944129

cd ${workspaceDir}

gzip -r "raw_reads"

############################
### PART 3: FASTP, QUAST ###
############################

cd ${rawReads}
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

cd ${workspaceDir}

#######################
## PART 2: MetaPhlAn ##
#######################

# Activate env with metaphlan
source /software/anaconda2/etc/profile.d/conda.sh
conda activate metaphlan-new-3.0
module load bowtie2

# Metaphlan profiling with trimmed reads
cd ${trimDir} 

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
mv *profile.txt ${metaphlanDir}
mv *bowtie2* ${metaphlanDir}

cd ${metaphlanDir}

# Generate single tab-delimited table from files
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt

# Species only table and heatmap
grep -E "s__|clade" merged_abundance_table.txt | sed 's/^.*s__//g'\ | cut -f1,3-8 | sed -e 's/clade_name/body_site/g' > merged_abundance_table_species.txt
hclust2.py -i merged_abundance_table_species.txt -o abundance_heatmap_species.png --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
