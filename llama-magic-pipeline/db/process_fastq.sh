#!/bin/bash
#SBATCH --array=0-1 ## 2 DBs to process 
#SBATCH --output=%A_%a.o
#SBATCH --error=%A_%a.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=480:00:00 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kmolloy@rockefeller.edu

module purge
module load fastqc
module load trimmomatic
module load abyss/2.1.1

DIR=$1 
HINGE_PRIMERS=$2

if [ $HINGE_PRIMERS == 0 ]; then
	echo "Normal primers"  
	TRIM_LEN_R1=250
	TRIM_LEN_R2=200
	ABYSS_M=50
else
	echo "Hinge primers"
	TRIM_LEN_R1=275
	TRIM_LEN_R2=250
    ABYSS_M=25
fi

TRIM_STRICT=0.8
ABYSS_P=0.95

#FILES=($(ls $DIR/*fastq.gz | perl -pe 's/^.+[^\-]+\-(\d+)\-.+\_R[12].+\.fastq\.gz/$1/g' | sort | uniq))
FILES=($(ls $DIR/*fastq.gz | perl -pe 's/^.+[^\-\_]+[\-\_](\d+[\-\_][BH]).+\_R[12].+\.fastq\.gz/$1/g' | sort | uniq))


# the file to run 
ID=${FILES[$SLURM_ARRAY_TASK_ID]}
INPUT_R1=$(ls $DIR/*$ID*_R1*.fastq.gz)
INPUT_R2=$(ls $DIR/*$ID*_R2*.fastq.gz)

echo "Running pipeline for the files: ${INPUT_R1} ${INPUT_R2}"

# create a dir for output
echo "Making dir for output: $DIR/$ID"
mkdir $DIR/$ID
cd $DIR/$ID

#fastqc and trimmomatic and fastqc and merge and fastqc
echo "Running FASTQC before TRIM."
fastqc -o . ${INPUT_R1}
fastqc -o . ${INPUT_R2}

echo "Running TRIM."
java -jar /gpfs/share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE ${INPUT_R1} ${ID}_R1_trimmed.fastq.gz MAXINFO:$TRIM_LEN_R1:$TRIM_STRICT
java -jar /gpfs/share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE ${INPUT_R2} ${ID}_R2_trimmed.fastq.gz MAXINFO:$TRIM_LEN_R2:$TRIM_STRICT

echo "Running FASTQC after TRIM."
fastqc -o . ${ID}_R1_trimmed.fastq.gz
fastqc -o . ${ID}_R2_trimmed.fastq.gz

echo "Run Abyss merge pairs."
abyss-mergepairs -o $DIR/$ID/$ID -p $ABYSS_P -m $ABYSS_M --standard-quality ${ID}_R1_trimmed.fastq.gz ${ID}_R2_trimmed.fastq.gz

echo "Running FASTQC after Merge."
fastqc -o . ${ID}_merged.fastq
fastqc -o . ${ID}_reads_1.fastq
fastqc -o . ${ID}_reads_2.fastq

#convert fastq to fasta
echo "Convert fastq to fasta"
mkdir $DIR/$ID/merged
mkdir $DIR/$ID/merged/fasta
mv $DIR/$ID/${ID}_merged.fastq $DIR/$ID/merged/${ID}_merged.fastq
perl $DIR/convert_fastq_to_fasta.pl $DIR/$ID/merged $DIR/$ID/merged/fasta  ### this one okay but maybe fix

#detect/separate primers sequences 
echo "Detect primers (to help with finding the correct reading frame)."
mkdir $DIR/$ID/merged/fasta/removed
if [ $HINGE_PRIMERS == 0 ]; then
	perl $DIR/remove_primers.pl $DIR/$ID/merged/fasta $DIR/$ID/merged/fasta/removed  ### need to fix the pl script and then fix here (only ...new_hinge.pl fixed so far)
else 
	perl $DIR/remove_primers_new_hinge.pl $DIR/$ID/merged/fasta/${ID}_merged.fastq-1.fasta $DIR/$ID/merged/fasta/removed/dna_removed_primers.fasta
fi

#translate to protein sequence
echo "Translate to protein sequence"
mkdir $DIR/$ID/merged/fasta/removed/translated
perl $DIR/dnatoprot_longest_nr_with_primers.pl $DIR/$ID/merged/fasta/removed/dna_removed_primers.fasta none $DIR/$ID/merged/fasta/removed/translated 1 $HINGE_PRIMERS 0 1

#digest in silico
echo "Digest in silico (trypsin)"
perl $DIR/digest_fasta.pl $DIR/$ID/merged/fasta/removed/translated/longest_nr.fasta $DIR/$ID/merged/fasta/removed/translated 1 0 0

echo "Digest in silico (chymotrypsin)"
perl $DIR/digest_fasta.pl $DIR/$ID/merged/fasta/removed/translated/longest_nr.fasta $DIR/$ID/merged/fasta/removed/translated 4 0 1

echo "Pipeline complete."


# NOTE: 'use tail' sequences with 'new primers' is not compatible; 'use tail' has not been used since 'new primers' introduced

# Additional steps after this script completes:
# --------------------------------------------

# create DIR structure for each set of fastq files

# MOVE longest_nr.fasta, all_predigested.fasta, all_predigested_chymotrypsin.fasta, protein_peptides.fasta and protein_peptides_chymotrypsin.fasta to pipeline DIR 
# change names of all_predigested... to longest_nr_predigested...
# change permissions of files

# copy/correct taxonomy.xml and taxonomy_chymotrypsin.xml files from a previous pipeline DIR

# copy/correct db_params.txt from a previous pipeline DIR








