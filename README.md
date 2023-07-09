# Llama_Nanobodies

This repository contains scripts for running the Llama Magic pipeline.  The text in bold gives examples of program execution on the command line.

The first step in the pipeline is to prepare the in silico translated protein sequence database from HT-sequencing (FASTQ files) of the VhH variable region of the nanobody sequence from Bone Marrow aspirates of the immunized animal.  

The bash script db/process_fastq.sh is an example of a script to run the database steps on an HPC cluster, but will need to be modified for your particular system.

(1) Perform QC, trimming and merge read pairs.  We use Trimmomatic for trimming Abyss mergepairs for merging the reads:

**java -jar trimmomatic-0.36.jar SE read1.fastq.gz read1_trimmed.fastq.gz MAXINFO:275:0.8**

**java -jar trimmomatic-0.36.jar SE read2.fastq.gz read2_trimmed.fastq.gz MAXINFO:250:0.8**

**abyss-mergepairs -o my_dir -p 0.95 -m 25 --standard-quality read1_trimmed.fastq.gz read2_trimmed.fastq.gz**

_We also recommend performing QC of the reads with fastqc before/after trimming and after merging_
 
The following steps are performed using a series of custom perl scripts:

(2) Convert merged FASTQ file to FASTA format:

**perl convert_fastq_to_fasta.pl my_dir my_output_dir**

_First parameter to the script is the directory where the fastq files are located.  Second parameter is the directory to place the fasta files.  The script will process all fastq files in the directory._

(3) Search for the PCR primers near the 5’ and 3’ ends of each merged sequence.  This helps to determine the correct ORF for in silico translation, since the correct reading frame from the start position of the primer is known, and their orientation determines forward or reverse strandedness.  This script adds the primer information to the header line of the FASTA file.

**perl remove_primers_new_hinge.pl dna.fasta dna_removed_primers.fasta** 

_First parameter to the script is the sequence database in FASTA format and the second parameter is the file name of the output file.  This script is specifically for the new "hinge" primers in use as of LM 2.0_ 

(4) Perform an in silico translation to protein sequence.  This script uses the output from step 3 to correctly perform the translation to protein sequence. 

**perl dnatoprot_longest_nr_with_primers.pl dna_removed_primers.fasta none output_dir 1 1 0 1**

_The parameters are as follows:_

- _dna_removed_primers.fasta: file to process (output of previous step)_
- _none: indicates no 2nd file to process (the script will optionally process 2 files and combine results)_
- _output_dir: the location for the translated protein sequence database file(s)_
- _1: 1/0 flag indicating if primers should be used in the translation (otherwise longest non-redundant ORF is chosen as the protein sequence)_
- _1: 1/0 flag indicating primer type (a value of 1 indicates the new "hinge" primers in use as of LM 2.0)_
- _0: 1/0 flag indicating if singleton reads (protein sequences with copy number == 1) should be filtered from the DB (experimental, recommended to keep as 0)_
- _1: 1/0 flag indicating if in the case of singleton reads filtering, to output both filtered and non-filtered version of the DB (keep as 1)_

(5) Perform an in silico digestion of protein sequences with trypsin/chymotrypsin.  One or both digestions can be performed.  Allowed missed cleavages = 1 (trypsin) and 4 (chymotrypsin).  This script produces a FASTA file where each unique (chymo-)tryptic peptide is listed on a separate line.  It also produces an index file for input to the peptide mapping step (step 8).

**perl digest_fasta.pl longest_nr.fasta output_dir 1 0 0**  
(trypsin)

**perl digest_fasta.pl longest_nr.fasta output_dir 4 0 1**  
(chymotrypsin)

_The parameters (for trypsin digest) are as follows:_

- _longest_nr.fasta: protein sequence database file to process (output of step 4)_
- _output_dir: the location for the digested peptides database file(s)_
- 1: number of missed cleavages allowed (1 for trypsin and 4 for chymotrysin)
- 0: 1/0 flag no longer used (keep at 0)
- 0: 1/0 flag indicating digestion with trypsin (0) or chymotrypsin (1)

<hr>

Once the protein sequence database has been prepared, it is ready to be searched using MS/MS results data.  This involves performing X!Tandem peptide identification, mapping of identified peptides back to the protein sequences, and CDR finding and scoring of sequences:

We recommend setting up this pipeline on an HPC cluster.  The scripts "llama_pipeline.pl", "run_LM_script_pipeline.pl" and "run_LM_script.sh" are used in our implementation on an HPC cluster and are provided as an example but wil need to be modified for your particular system.

(6) Run X!Tandem: input is the Mass Spec data (mgf file) and the pre-digested peptide database created above.  The input parameters for X!Tandem are specified in XML format.  Examples of the input XML files can be found in the repository.  The main input file (input1.xml) references "taxonomy.xml" or "taxonomy_chymotrypsin.xml" which point to the pre-digested database files.  It also references "default_input.xml" containing default input parameters for XTandem.  

**tandem.exe input1.xml**

(7) Parse X!Tandem output xml file to a tab delimited text file with a peptide/spectra match on each line.  The script will filter based on statistical confidence (expectation values) log(e) (e = Expectation).

**perl parse_xtandem_llama.pl results 0.1**

_First parameter to the script is the directory where the X!Tandem output xml file is located.  Second parameter is the log(e) Expectation cutoff value.  The script will process all xml files in the directory_

(8) Map peptides to protein sequences, identify CDR regions and score sequences based on coverage.

**perl map_peptides_to_proteins.pl "output1.xml.peptide_list.0.1.txt" "" "longest_nr.fasta" "output1.xml" "" "protein_peptides.fasta" "" "1" "1" "0" "1"**

_The parameters are as follows:_

- _"output1.xml.peptide_list.0.1.txt": parsed X!Tandem output file_
- _"": 2nd parsed X!Tandem output file (left blank in example, can be used if performing both trypsin and chymotrypsin digest and results are to be combined_
- _"longest_nr.fasta": protein sequence database file (full sequences, not pre-digested)_
- _"output1.xml": X!Tandem output file_
- _"": 2nd X!Tandem output file (left blank in example, can be used if performing both trypsin and chymotrypsin digests)_
- _"protein_peptides.fasta": sequence to peptide look up (index) file allows for faster run time (this is an output of step 5 above)_
- _"": 2nd look up file (i.e. for chymotrypsin)_
- _"1": 1/0 flag to show sequence "score" (do not change)_
- _"1": 1/0 flag to use primer information from in silico translation for CDR positioning (do not change)_
- _"0": 1/0 flag no longer used (keep at 0)_
- _"1": 1/0 flag indicating primer type (a value of 1 indicates the new "hinge" primers in use as of LM 2.0)_

<hr>

The output of this script is HTML file(s) that can be opened with the Llama Magic web viewer.  These HTML files display the protein sequences (ordered by score) annotated with each identified peptide and the CDR regions.  Please see the documentation at https://hub.docker.com/r/snkeegan/llama-magic for details on setting up the LM web viewer.  It uses a local web server to display the HTML-formatted results and allows for a detailed view of each peptide identification.  

It is also possible to simply open the HTML files in any web browser but without the LM viewer application some functions will not be accessible.
