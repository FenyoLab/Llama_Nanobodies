#!c:/perl/bin/perl.exe 

use strict;
#use warnings;
use POSIX;
#use JSON;

my $img_source_plus = '/llama-magic-html/plus.gif';
my $img_source_minus = '/llama-magic-html/minus.gif';
my $img_source_star = '/llama-magic-html/greyplus.gif';

my $filename="";
my $filename2="";
my $fileroot;
my $filedir;
my $db_filename="";
my $tandem_output_filename="";
my $tandem_output_filename2="";
my $PEPTIDE_FILTER_MIN = 1.2; 
my $MIN_SEQ_LENGTH = 107;
my $FILTER_FIX_SEQ_END=1; #when primers not used
my $NUM_AA_AFTER_QVT=6; #when primers not used
my $MIN_CDR3_COV_PERC = 25;  
my $MAX_PROTEINS_IN_FILE = 5000;
my $MAX_GROUPS_IN_FILE = 500;
my $MAX_GROUP_PROTEINS_IN_FILE = 10; #25;
my $RK_FIXATION = 0; # NOTE: RK_FIXATION SHOULD BE SET TO 0!! [B/C CHYMOTRYPSIN ]
my @TAIL_SEQUENCES = ("SEPKIPQPQPKPQ", "SAHHSEDPSSKCP", "SEPKTPKPQPQPQPQ", "SGTNEVCKCPKCPAPEL", "EPKIPQPQPKPQ", "AHHSEDPSSKCP", "EPKTPKPQPQPQPQ", "GTNEVCKCPKCPAPEL");
my $NEW_CDR3_METHOD = 1;
my $incompletes = 1;

my $check_for_VH = 0;

my @nanobody_group; #the nanobodys are grouped by CDR3 region, allowing one mismatch

sub log10_
{
	my $n = shift;
	return log($n)/log(10);
}

my $cdr3_not_found=0;

##### ONLY FOR MAYO (HRB) DATA 02_20_2018
my $ADD_MAQ_TO_START_CHANGE_QVT_CHECK = 0; #if 1, then changes in filter_input, and also in CDR3-finding
my $QVT_LEN_FROM_END=25; #25 #-1; #in filter_input: if this is -1, won't filter on QVT at all, if this is > 0, then will look for QVT this amount from end
#####

##
#changes for llama_magic web app:
#changed protein database to file name not directory (ARGV[1])
#put log file and output file in directory where tandem results txt file is located (ARGV[0]) rather than in current directory of script
#output html file name always includes a '1' even if there is only one output file (see function open_OUT_ALL)
#put link at bottom of output html file to go to the next output file (if it exists)
##
my $output_excel = 1; #if true, the program will out put LM resuts to an excel file
my $stand_alone = 0; #this var will determine how we examine the input file containing XTamdem results - is it directly from XTAndem - i.e. run
		     #through llama-magic OR was it manually processed and provided to us - in each case the format is different
my $show_score;
my $use_primers;
my $new_primers;
my $use_tail;
my $peptide_index_filename;
my $peptide_index_filename2;

my $create_cdr3_fasta_all_seq = 0;
my $create_cdr3_fasta_by_cov = 0;

my $filter_cdr3 = 0;
my $filter_cdr3_filename = "output.xml.peptide_list.0.1.cdr3_cov75.fasta";

if ($ARGV[0]=~/\w/) { $filename="$ARGV[0]"; } 
else { $filename="/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/Llama/new_test/output1.xml.peptide_list.0.1.txt"; }

if ($ARGV[1]=~/\w/) { $filename2="$ARGV[1]"; } #this is for the second XTandem result file (When combining results as in trypsin+chymotrypsin)
else { $filename2=""; }

if ($ARGV[2]=~/\w/) { $db_filename="$ARGV[2]"; }
else { $db_filename="/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/Llama/new_test/longest_nr.fasta"; }

if ($ARGV[3]=~/\w/) { $tandem_output_filename="$ARGV[3]"; }
else { $tandem_output_filename="/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/Llama/new_test/output1.xml"; }

if ($ARGV[4]=~/\w/) { $tandem_output_filename2="$ARGV[4]"; }
else { $tandem_output_filename2=""; }

if ($ARGV[5]=~/\w/) { $peptide_index_filename="$ARGV[5]"; } #if file exists, index will be used else index won't be used
else{ $peptide_index_filename = "/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/Llama/new_test/protein_peptides.fasta"; }

if ($ARGV[6]=~/\w/) { $peptide_index_filename2="$ARGV[6]"; } #if file exists, index will be used else index won't be used
else{ $peptide_index_filename2 = ""; }

if ($ARGV[7]=~/\w/) { $show_score=$ARGV[7]; }
else { $show_score = 1; }

if ($ARGV[8]=~/\w/) { $use_primers=$ARGV[8]; }
else { $use_primers = 1; }

if ($ARGV[9]=~/\w/) { $use_tail=$ARGV[9]; }
else { $use_tail = 0; } 

if($ARGV[10]=~/\w/) { $new_primers=$ARGV[10]; }
else { $new_primers=1; }

#primer params
my $P1_SEQ_LENGTH;
my $P2_LH_SEQ_LENGTH;
my $P2_SH_SEQ_LENGTH;
my $P2_SEQ_LENGTH;
if($new_primers)
{
	$P1_SEQ_LENGTH = 8;
	$P2_SH_SEQ_LENGTH = 6; 
	$P2_LH_SEQ_LENGTH = 8; 
	$P2_SEQ_LENGTH = $P2_LH_SEQ_LENGTH;
}
else
{
	$P1_SEQ_LENGTH = 10;
	$P2_SEQ_LENGTH = 6;
}

my $cluster = 0;
if ($cluster)
{
	$stand_alone = 0;
	$show_score = 1;
	$output_excel = 1;
	$MIN_CDR3_COV_PERC = 1;
}

$filedir = $filename;
$filedir =~ s/[\/\\][^\/\\]+$//;

$filename =~ /[\/\\]([^\/\\]+)$/;
$fileroot = $1;
$fileroot =~ s/(\.\w\w\w)$//;
$fileroot =~ s/^output1/output/;

#open log file...
if(!open(LOG,">$filedir/$fileroot.count_proteins.log.txt"))
{
	print "Could not open for writing: $filedir/$fileroot.count_proteins.log.txt ($!)\n";
	exit(1);
}

#autoflush log file, so we an check progress more easily
select(LOG);
$|++; # autoflush LOG
select(STDOUT);

#output the parameters to log file
print LOG "filename=$filename\n";
print LOG "filename2=$filename2\n";
print LOG "db_filename=$db_filename\n";
print LOG "tandem_output_filename=$tandem_output_filename\n";
print LOG "tandem_output_filename2=$tandem_output_filename2\n";
print LOG "peptide_index_filename=$peptide_index_filename\n";
print LOG "peptide_index_filename2=$peptide_index_filename2\n";
print LOG "use_tail=$use_tail\n";
print LOG "show_score=$show_score\n";
print LOG "use_primers=$use_primers\n";
print LOG "new_primers=$new_primers\n";

###########Reading in peptide data################################################################
my %peptides=(); #the peptides, $peptides{seq} = 1

print LOG "Reading in peptide data...\n";
#read in database search results - the peptides that were found, and the expectation values
read_in_xtandem_results($filename, $tandem_output_filename);
if ($filename2)
{
    read_in_xtandem_results($filename2, $tandem_output_filename2);
}

#filter input peptides by expectation score
#filter: -log(e) * peptides_count >= 1.2 (peptides_count = # times peptide found in the file, log(e) = lowest expect value for peptide)
foreach (keys %peptides)  # will pass filter if: -log(e) > 1.2 OR found > 1x in file
{
	if(abs($peptides{$_}[0] * $peptides{$_}[3]) < $PEPTIDE_FILTER_MIN)
	{#remove peptide from list
		delete $peptides{$_};
	}
}

###########Reading in protein data###############################################################
my %proteins=(); #the proteins, $proteins{name} = seq
my %protein_gene_counts;
my $count_total_seq = 0;

print LOG "Reading in protein database data...\n";
$count_total_seq=read_in_protein_db($db_filename);

################Read in protein to peptide index##################################################
#if file exists, index will be used else index won't be used
my %peptide_index = ();
my $use_index = read_in_peptide_index($peptide_index_filename);
if ($use_index && $filename2 ne "")
{
    $use_index = read_in_peptide_index($peptide_index_filename2);
}

######Read in CDR3's that should be filtered out due to overlay with other animal DBs#############
my %filter_cdr3 = ();
if ($filter_cdr3)
{
	if(!open(FILTER_IN,"$filter_cdr3_filename"))
	{
		print LOG "Error: Could not open '$filter_cdr3_filename' ($!)\n";
		print LOG "No CDR3 filtering will not be used.\n";
		$filter_cdr3 = 0;
	}
	else
	{
		my $name;
		while ($name=<FILTER_IN>)
		{
			my $cur_cdr3 = <FILTER_IN>;
			
			$cur_cdr3 =~ s/\r\n$//;
			$cur_cdr3 =~ s/\n$//;
			#chomp($cur_cdr3);
			
			$filter_cdr3{$cur_cdr3} = 1;
		}
		
	}
}

################Find the 3 CDR regions and store their positions:##################################
my %CDR1; #$CDR1{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR2; #$CDR2{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR3; #$CDR3{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR3_seq;
my %CDR3_name;

#################Match peptides to proteins########################################################
my %peptide_proteins; #the proteins found (in the fasta dbs) that contain this peptide:
					  # peptide_proteins{pep-seq} = #prot-name1##prot-name2#...
my %peptide_proteins_lookup;
my %peptide_proteins_count; #the number of proteins found that match this peptide - protein_peptides_count{pep-seq} = n
my %protein_peptides; #the peptides found in this protein - protein_peptides{prot-name} - #pep-seq1##pep-seq2#...
my %protein_peptides_count; #the number of peptides found in this protein - protein_peptides_count{prot-name} = n
my %protein_peptides_pos;
my %protein_total_coverage_by_aa;
my %protein_cdr1_cover;
my %protein_cdr2_cover;
my %protein_cdr3_cover;
my %protein_total_cover;
my %protein_max_cover_tail;
my @protein_coverage_stats;

print LOG "Matching peptides to proteins...\n";

if ($create_cdr3_fasta_all_seq)
{
	if(!open(CDR3_FASTA,">$filedir/$fileroot.cdr3.fasta"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.cdr3.fasta ($!)\n";
	}
}

if ($create_cdr3_fasta_by_cov)
{
	if(!open(CDR3_FASTA_25, ">$filedir/$fileroot.cdr3_cov25.fasta"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.cdr3_cov25.fasta ($!)\n";
	}
	if(!open(CDR3_FASTA_50, ">$filedir/$fileroot.cdr3_cov50.fasta"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.cdr3_cov50.fasta ($!)\n";
	}
	if(!open(CDR3_FASTA_75, ">$filedir/$fileroot.cdr3_cov75.fasta"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.cdr3_cov75.fasta ($!)\n";
	}
	
}
	
my $count_proteins = 0; my $stat_i = 0;
my $start = time;
my %aa_pos_37; my %aa_pos_44; my %aa_pos_45; my %aa_pos_47; my $count_cdr1_diff_pos=0;
foreach my $name (keys %proteins)
{
	#CDR finding
	my $cdr1_l_pos; my $cdr1_r_pos;	my $cdr2_l_pos; my $cdr2_r_pos;	my $cdr3_l_pos; my $cdr3_r_pos;
	
	find_cdr1($proteins{$name}, $cdr1_l_pos, $cdr1_r_pos);
	find_cdr2($proteins{$name}, $cdr2_l_pos, $cdr2_r_pos);
	if ($NEW_CDR3_METHOD)
	{
		find_cdr3_new($proteins{$name}, $cdr2_r_pos, $cdr3_l_pos, $cdr3_r_pos, $name);
	}
	else { find_cdr3($proteins{$name}, $cdr3_l_pos, $cdr3_r_pos); }
	
	if ($check_for_VH)
	{
		#get stats for AAs in positions 37,44,45,47 from Llama VH/VHH paper
		#note in our numbering, the positions are 2 further back but also we index by zero so only add 1
		my $change_pos = $cdr1_r_pos - 36;
		my $aa = substr($proteins{$name}, 38+$change_pos, 1);
		if (defined $aa_pos_37{$aa}) { $aa_pos_37{$aa} += 1; }
		else { $aa_pos_37{$aa} = 1; }
		
		$aa = substr($proteins{$name}, 45+$change_pos, 1);
		if (defined $aa_pos_44{$aa}) { $aa_pos_44{$aa} += 1; }
		else { $aa_pos_44{$aa} = 1; }
		
		$aa = substr($proteins{$name}, 46+$change_pos, 1);
		if (defined $aa_pos_45{$aa}) { $aa_pos_45{$aa} += 1; }
		else { $aa_pos_45{$aa} = 1; }
		
		$aa = substr($proteins{$name}, 48+$change_pos, 1);
		if (defined $aa_pos_47{$aa}) { $aa_pos_47{$aa} += 1; }
		else { $aa_pos_47{$aa} = 1; }
		
		if ($cdr1_r_pos != 36)
		{
			print LOG "NOTE! CDR1 does not end at expected position: CDR1-end-pos = $cdr1_r_pos, Sequence = $proteins{$name}\n";
			$count_cdr1_diff_pos++;
		}
	}
	
	my $peptides_to_map = {};
	my %indexed_peptides = ();
	if ($use_index)
	{
		if (defined $peptide_index{$name})
		{ #peptides from this protein match peptides from xtandem results, they will be mapped
			$peptides_to_map = \%{$peptide_index{$name}};
		}
		else
		{#no peptides from this protein match
			$peptides_to_map = {};
		}
	}
	else # without index, must try to map all the XTandem matched peptides to the protein
	{ $peptides_to_map = \%peptides; }
	
	my $seq = $proteins{$name};
	
#	if ($FILTER_FIX_SEQ_END && !$use_primers) # ** do this even if $use_primers? or if $new_primers?
#	{#remove X's that were put on the sequence to standardize length and correctly identify the CDR regions
#		$seq =~ s/X+$//;
#	}
	$seq =~ s/X+$//;
	$seq =~ s/^X+//;
	
	my $orig_seq_len = length($seq);
	my $total_seq_len = $orig_seq_len;
	
	my %protein_coverage = ();
	my $protein_peptides='';
	my $protein_peptides_count=0;
	my $protein_peptides_pos='';
	my $max_cover_tail = ''; 
	if (scalar keys % { $peptides_to_map } > 0)
	{
		#mapping peptides to proteins
		$seq=~tr/L/I/;
		
		if ($use_tail)
		{
			my $max_cover = 0;
			my $max_cover_protein_coverage_ref;
			foreach my $tail (@TAIL_SEQUENCES)
			{
				$tail=~tr/L/I/;
				
				my %protein_coverage_tail=();
				my $protein_peptides_tail;
				my $protein_peptides_count_tail;
				my $protein_peptides_pos_tail;
				map_peptides($seq.$tail, \%protein_coverage_tail, $protein_peptides_tail, $protein_peptides_count_tail,
							 $protein_peptides_pos_tail, $peptides_to_map);
				
				if ((scalar keys %protein_coverage_tail) > $max_cover)
				{
					$max_cover = scalar keys %protein_coverage_tail;
					$max_cover_tail = $tail;
					$max_cover_protein_coverage_ref = \%protein_coverage_tail;
					
					$protein_peptides = $protein_peptides_tail;
					$protein_peptides_count = $protein_peptides_count_tail;
					$protein_peptides_pos = $protein_peptides_pos_tail;
				}
			}
			if ($max_cover == 0) { $max_cover_tail = $TAIL_SEQUENCES[0]; }
			$total_seq_len += length($max_cover_tail);
			$protein_peptides{$name} = $protein_peptides;
			$protein_peptides_count{$name} = $protein_peptides_count;
			$protein_peptides_pos{$name} = $protein_peptides_pos;
			foreach (keys %$max_cover_protein_coverage_ref)
			{
				$protein_coverage{$_} = 1;
			}
		}
		else
		{
			map_peptides($seq, \%protein_coverage, $protein_peptides, $protein_peptides_count, $protein_peptides_pos, $peptides_to_map);
			$protein_peptides{$name} = $protein_peptides;
			$protein_peptides_count{$name} = $protein_peptides_count;
			$protein_peptides_pos{$name} = $protein_peptides_pos;
		}	
	}
	else
	{
		if ($use_tail) { $max_cover_tail = $TAIL_SEQUENCES[0]; }
		$protein_peptides{$name} = '';
		$protein_peptides_count{$name} = 0;
		$protein_peptides_pos{$name} = '';
	}
	
	#calculate coverage percent for this protein (for sorting and display) ( including tail)
	my $cover_n = 0; my $cover_cdr = 0; my $cover_cdr2 = 0; my $cover_cdr1 = 0; my $cover_cdr3 = 0;
	foreach (keys %protein_coverage)
	{
		if($cdr1_l_pos <= $_ && $cdr1_r_pos >= $_) { $cover_cdr++; $cover_cdr1++; }
		if($cdr2_l_pos <= $_ && $cdr2_r_pos >= $_) { $cover_cdr++; $cover_cdr2++; }
		if($cdr3_l_pos <= $_ && $cdr3_r_pos >= $_) { $cover_cdr++; $cover_cdr3++; }
		$cover_n++;
	}
	
	my $temp = $protein_peptides{$name};
	while($temp=~s/^#([^#]+)#//)
	{#read in peptides for this protein
		my $pep = $1;
		$peptide_proteins{$pep}.="$name<br>";
		$peptide_proteins_lookup{$pep}{$name} = 1;
		$peptide_proteins_count{$pep}++;
	}
	
	#total 3 cdrs length
	my $cdr_len = ($cdr1_r_pos - $cdr1_l_pos + 1) + ($cdr2_r_pos - $cdr2_l_pos + 1) + ($cdr3_r_pos - $cdr3_l_pos + 1);
	my $cdr3_len = $cdr3_r_pos - $cdr3_l_pos + 1;
	my $cdr3_cover_perc = ($cover_cdr3 / $cdr3_len) * 100;
	my $cdr3_seq = substr($proteins{$name}, $cdr3_l_pos, $cdr3_r_pos-$cdr3_l_pos+1);
	
	if ($cdr3_cover_perc >= $MIN_CDR3_COV_PERC && $cover_n > 0) 
	{
		#create the sorting array:
		$protein_coverage_stats[$stat_i][0] = sprintf("%.1f", $cdr3_cover_perc);
		$protein_coverage_stats[$stat_i][1] = $cdr3_len;
		$protein_coverage_stats[$stat_i][2] = $protein_gene_counts{$name};
		$protein_coverage_stats[$stat_i][3] = sprintf("%.1f", ($cover_cdr / $cdr_len) * 100);
		$protein_coverage_stats[$stat_i][4] = $cdr_len;
		$protein_coverage_stats[$stat_i][5] = sprintf("%.1f", ($cover_n / $total_seq_len) * 100);
		$protein_coverage_stats[$stat_i][6] = $total_seq_len;
		$protein_coverage_stats[$stat_i][7] = $name;
		
		my $s = 8*($cover_cdr3 / $cdr3_len) +
				2*($cover_cdr2 / ($cdr2_r_pos - $cdr2_l_pos + 1)) +
				2*($cover_cdr1 / ($cdr1_r_pos - $cdr1_l_pos + 1)) +
				2*($cover_n / $total_seq_len) +
				$cdr3_len/15 +
				($cdr2_r_pos - $cdr2_l_pos + 1)/10 +
				($cdr1_r_pos - $cdr1_l_pos + 1)/9;
		
		$protein_coverage_stats[$stat_i][8] = sprintf("%.1f", $s);
		$stat_i++;
		
		#store more info about the protein for display:
		foreach my $key (keys %protein_coverage)
		{
			$protein_total_coverage_by_aa{$name}{$key} = 1;
		}
		$protein_cdr1_cover{$name} = $cover_cdr1;
		$protein_cdr2_cover{$name} = $cover_cdr2;
		$protein_cdr3_cover{$name} = $cover_cdr3;
		$protein_total_cover{$name} = $cover_n;
		$protein_max_cover_tail{$name} = $max_cover_tail;
		
		$CDR1{$name} = [$cdr1_l_pos, $cdr1_r_pos];
		$CDR2{$name} = [$cdr2_l_pos, $cdr2_r_pos];
		$CDR3{$name} = [$cdr3_l_pos, $cdr3_r_pos];
		$CDR3_seq{$name} = $cdr3_seq;
	}
	if ($create_cdr3_fasta_all_seq)
	{
		if($cdr3_len > 2)
		{ print CDR3_FASTA ">$name $proteins{$name}\n$cdr3_seq\n"; }
	}
	if ($create_cdr3_fasta_by_cov)
	{
		if ($cdr3_cover_perc >= 25)
		{
			if($cdr3_len > 2)
			{ print CDR3_FASTA_25 ">$cdr3_seq\n$cdr3_seq\n"; }
		}
		if ($cdr3_cover_perc >= 50)
		{
			if($cdr3_len > 2)
			{ print CDR3_FASTA_50 ">$cdr3_seq\n$cdr3_seq\n"; }
		}
		if ($cdr3_cover_perc >= 75)
		{
			if($cdr3_len > 2)
			{ print CDR3_FASTA_75 ">$cdr3_seq\n$cdr3_seq\n"; }
		}
	}
	
	if($count_proteins % 50000 == 0) { print LOG qq!$count_proteins ($count_total_seq)\n!; }
	if($count_proteins % 5000 == 0 && $count_proteins > 0)
	{
		my $duration = time - $start;
		print qq!$count_proteins ($count_total_seq) ($duration s)\n!;
		
		$start = time;
	}
	$count_proteins++;
}

if ($check_for_VH)
{
	print LOG "VH / VHH ANALYSIS:\n";
	print LOG "Position 37: \n";
	my $total = 0;
	foreach my $aa (keys %aa_pos_37)
	{
		my $perc = sprintf("%.2f", ($aa_pos_37{$aa}/$count_total_seq)*100);
		print LOG "$aa: $aa_pos_37{$aa}/$count_total_seq ($perc)\n";
		$total += $aa_pos_37{$aa};
	}
	print LOG "Total $total\n";
	print "\n";
	
	print LOG "Position 44: \n";
	$total = 0;
	foreach my $aa (keys %aa_pos_44)
	{
		my $perc = sprintf("%.2f", ($aa_pos_44{$aa}/$count_total_seq)*100);
		print LOG "$aa: $aa_pos_44{$aa}/$count_total_seq ($perc)\n";
		$total += $aa_pos_44{$aa};
	}
	print LOG "Total $total\n";
	print "\n";
	
	print LOG "Position 45: \n";
	$total = 0;
	foreach my $aa (keys %aa_pos_45)
	{
		my $perc = sprintf("%.2f", ($aa_pos_45{$aa}/$count_total_seq)*100);
		print LOG "$aa: $aa_pos_45{$aa}/$count_total_seq ($perc)\n";
		$total += $aa_pos_45{$aa};
	}
	print LOG "Total $total\n";
	print "\n";
	
	print LOG "Position 47: \n";
	$total = 0;
	foreach my $aa (keys %aa_pos_47)
	{
		my $perc = sprintf("%.2f", ($aa_pos_47{$aa}/$count_total_seq)*100);
		print LOG "$aa: $aa_pos_47{$aa}/$count_total_seq ($perc)\n";
		$total += $aa_pos_47{$aa};
	}
	print LOG "Total $total\n";
	print "\n";
	
	print LOG "Number sequences where CDR1 does not end in 'usual' position: $count_cdr1_diff_pos.\n";
}

if ($create_cdr3_fasta_all_seq) { close(CDR3_FASTA); }
if ($create_cdr3_fasta_by_cov)
{
	close(CDR3_FASTA_25);
	close(CDR3_FASTA_50);
	close(CDR3_FASTA_75);
}

#my $num_CDR3s = scalar keys %CDR3_name;
#print LOG "Number of CDR3s (exact matching): $num_CDR3s\n";
#print LOG "Done!\n";

#save JSON
#my $utf8_encoded_json_text = encode_json \%CDR3_name;
#if (open (JSON_OUT,">C:/temp/cdr3_repeats.json"))
#{
#	print JSON_OUT $utf8_encoded_json_text;
#	close(JSON_OUT);
#}

############Sort proteins based on CDR coverage, length, etc.########################################
print LOG "Sorting the proteins for ouput...\n";

my @proteins_to_sort;
my $proteins_to_sort_count=0;
my @proteins_sorted;

if ($show_score)
{
	@proteins_sorted = sort { $b->[8] <=> $a->[8] || $b->[2] <=> $a->[2] } @protein_coverage_stats;
}
else
{
	@proteins_sorted = sort
	{ $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] || $b->[3] <=> $a->[3] || $b->[4] <=> $a->[4] || $b->[5] <=> $a->[5] || $b->[6] <=> $a->[6] }
	@protein_coverage_stats;
}

$proteins_to_sort_count = $#proteins_sorted + 1;

print LOG "Proteins sorted! ($proteins_to_sort_count)\n";

############Group the proteins by CDR3 region w/ one AA difference allowed, combine groups###########
print LOG "Grouping proteins on CDR3 region...\n";

my %lowest_matching_i;
for(my $i = 0; $i <= $#proteins_sorted; $i++)
{
	if($i % 100 == 0) { print LOG "$i...\n"; }
	if(defined ${$nanobody_group[$i]}{-1}) { next; }
	my %matching_is; my @matching_js;
	my $cdr3_seq_i = $CDR3_seq{$proteins_sorted[$i][7]};
	for(my $j = $i+1; $j <= $#proteins_sorted; $j++)
	{
		if(defined ${$nanobody_group[$j]}{-1}) { next; }
		my $cdr3_seq_j = $CDR3_seq{$proteins_sorted[$j][7]};
		if(compare_cdrs($cdr3_seq_i, $cdr3_seq_j))
		{#if cdrs match...
			if(defined $lowest_matching_i{$j})
			{#if protein j already matches a previous protein, record this protein's position 
				$matching_is{$lowest_matching_i{$j}} = 1;
			}
			push @matching_js, $j; #add protein j to list of matches to protein i
		}
	}
	if($#matching_js >= 0)
	{#if there are proteins that matched to protein i
		if(keys %matching_is)
		{#if there are other proteins that matched to the matches of protein i
			
			#get lowest group and add matching js
			my @sorted_matching_is = sort {$a <=> $b} keys %matching_is;
			my $lowest_i = $sorted_matching_is[0]; #the lowest matching protein group is the group to put all matches into
			
			foreach my $cur_j (@matching_js)
			{#for all matching j's to protein i, put them in the lowest i group
				${$nanobody_group[$lowest_i]}{$cur_j} = 1; #mark j as being in lowest matching protein group
				$lowest_matching_i{$cur_j} = $lowest_i; #set lowest matching i for this j
				${$nanobody_group[$cur_j]}{-1} = 1; #protein j is in another group, x out its entry
			}
			
			foreach my $cur_i (@sorted_matching_is)
			{#for each of the other proteins found that matched proteins matching protein i, we need to put them in the lowest i group
				if($cur_i != $lowest_i)
				{
					foreach my $cur_j (keys %{$nanobody_group[$cur_i]})
					{
						${$nanobody_group[$lowest_i]}{$cur_j} = 1;
						$lowest_matching_i{$cur_j} = $lowest_i;
					}
					${$nanobody_group[$cur_i]}{-1} = 1;
				}
			}
			
			#set the group for protein i as well
			${$nanobody_group[$lowest_i]}{$i} = 1;
			$lowest_matching_i{$i} = $lowest_i;
		}
		else
		{#its just new proteins that match, no groups
			foreach my $cur_j (@matching_js)
			{
				${$nanobody_group[$i]}{$cur_j} = 1; #mark protein j as being in protein i's group
				$lowest_matching_i{$cur_j} = $i; #set lowest matching i to the current i for protein j
				${$nanobody_group[$cur_j]}{-1} = 1; #protein j is in another group, x out its entry
			}
		}
	}
	else
	{
		$nanobody_group[$i] = ();
	}
}

#printing Stats and verifying groups...
my $num = $#nanobody_group+1;
print LOG "nanobody_group has $num members.\n";
$num = $#proteins_sorted+1;
print LOG "proteins_sorted has $num members.\n";
print LOG "cdr3_not_found=$cdr3_not_found\n";


#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	print LOG "$i: ";
#	if(defined ${$nanobody_group[$i]}{-1}) { print LOG "x"; }
#	else { print LOG join ', ', sort { $a <=> $b } keys %{$nanobody_group[$i]}; }
#	print LOG "\n";
#}

#do some checking
#print "Verifying groups...\n";
#my @checking_group;
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	if(!${$nanobody_group[$i]}{-1})
#	{
#		$checking_group[$i]++;
#		foreach (keys %{$nanobody_group[$i]})
#		{
#			$checking_group[$_]++;
#		}
#		
#	}
#}
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	if($checking_group[$i] != 1) { print "position $i of checking_group is checking_group[$i].\n"}
#}

if ($filter_cdr3)
{
	#remove groups that contain CDR3s that we want to filter
	for(my $i = 0; $i <= $#nanobody_group; $i++)
	{
		if(! (defined ${$nanobody_group[$i]}{-1})) 
		{
			#group head
			my $name = $proteins_sorted[$i][7];
			my $cdr3_seq = substr($proteins{$name}, $CDR3{$name}[0], $CDR3{$name}[1]-$CDR3{$name}[0]+1);
			if (defined $filter_cdr3{$cdr3_seq})
			{
				${$nanobody_group[$i]}{-1} = 1;
			}
			
			#group members
			if(! (defined ${$nanobody_group[$i]}{-1})) 
			{
				my @current_group = keys %{$nanobody_group[$i]};
				for(my $k = 0; $k <= $#current_group; $k++)
				{
					my $j = $current_group[$k];
					$name = $proteins_sorted[$j][7];
					$cdr3_seq = substr($proteins{$name}, $CDR3{$name}[0], $CDR3{$name}[1]-$CDR3{$name}[0]+1);
					if (defined $filter_cdr3{$cdr3_seq})
					{#remove entire group of which this sequence is a member
						${$nanobody_group[$i]}{-1} = 1;
						last;
					}
				}
				
			}
			
		}
	}
}

###### Outputting results to a txt file ############### <- only outputs data for > min coverage of CDR3 region - output all data above
if ($output_excel)
{
	#open log file...
	if(!open(STAT,">$filedir/$fileroot.count_proteins.group-stat.csv"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.count_proteins.group-stat.csv ($!)\n";
	}
	else
	{
		print STAT "RANK,ID,GROUP_HEAD,SCORE,CDR3_PERC,CDR3_LEN,HT_SEQ_COUNT,CDR123_PERC,CDR123_LEN,SEQ_PERC,SEQ_LEN,CDR1_PERC,CDR1_LEN,CDR2_PERC,CDR2_LEN,CDR3_SEQ,SEQ\n";
		my $head_count = 1;
		for(my $i = 0; $i <= $#nanobody_group; $i++)
		{
			
			if(! (defined ${$nanobody_group[$i]}{-1})) 
			{
				my $name = $proteins_sorted[$i][7];
				#print out stats for head of group
				
				my $cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
				my $cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
				my $cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
				my $cdr3_seq = substr($proteins{$name}, $CDR3{$name}[0], $CDR3{$name}[1]-$CDR3{$name}[0]+1);
				
				print STAT "$head_count,$name,1,$proteins_sorted[$i][8],$proteins_sorted[$i][0],$proteins_sorted[$i][1],$proteins_sorted[$i][2],$proteins_sorted[$i][3],$proteins_sorted[$i][4],$proteins_sorted[$i][5],$proteins_sorted[$i][6],$cdr1_perc,$cdr1_len,$cdr2_perc,$cdr2_len,$cdr3_seq,$proteins{$name}\n";
				
				#print out stats for group members
				my @current_group = sort { $a <=> $b } keys %{$nanobody_group[$i]};
				for(my $k = 0; $k <= $#current_group; $k++)
				{
					my $j = $current_group[$k];
					$name = $proteins_sorted[$j][7];
					
					###
					$cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
					$cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
					$cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
					$cdr3_seq = substr($proteins{$name}, $CDR3{$name}[0], $CDR3{$name}[1]-$CDR3{$name}[0]+1);
					print STAT "$head_count,$name,0,$proteins_sorted[$j][8],$proteins_sorted[$j][0],$proteins_sorted[$j][1],$proteins_sorted[$j][2],$proteins_sorted[$j][3],$proteins_sorted[$j][4],$proteins_sorted[$j][5],$proteins_sorted[$i][6],$cdr1_perc,$cdr1_len,$cdr2_perc,$cdr2_len,$cdr3_seq,$proteins{$name}\n";
				}
				$head_count += 1;
			}
		}
		close(STAT);
	}
}

##############Output results to the html file###################################################################################

my %protein_peptides_unique_count;
my %peptides_done;
my $proteins_unique_count=0;

print LOG "Outputting results to the file ($proteins_to_sort_count proteins with > $MIN_CDR3_COV_PERC% CDR3 coverage)...\n";

open_OUT_ALL(1);

if($proteins_to_sort_count == 0)
{
	print OUT_ALL "No matches found.\n";
	close_OUT_ALL(0, 1); 
	close(LOG);
	exit(0);
}
my @names_in_group = (); my %peptide_match_count = (); 
my $group_div = 0; my @current_group; my $group_div_head; my $indent = ""; my $make_new_file = 0; my $file_num = 1; my $num_group_members;
for(my $group_count = -1, my $cur_group_i = 0, my $real_group_count = 1, my $total_proteins_in_file = 0, my $real_group_count_on_page = 1;
	$group_count <= $#nanobody_group;
	$cur_group_i++, $total_proteins_in_file++)
{
	my $proteins_sorted_count; my $group_head;
	
	if($total_proteins_in_file > $MAX_PROTEINS_IN_FILE)
	{ $total_proteins_in_file = 0; $make_new_file = 1; }
	
	if($real_group_count_on_page > $MAX_GROUPS_IN_FILE)
	{ $real_group_count_on_page = 0; $total_proteins_in_file = 0; $make_new_file = 1; }
	
	if($cur_group_i > $#current_group || ($cur_group_i+1) >= $MAX_GROUP_PROTEINS_IN_FILE) #$cur_group_i starts at -1
	{#go to the next group
		if($cur_group_i <= $#current_group) { print OUT_ALL '          ...<br>'; } #there's more in the group, show ...
		$group_count++;
		
		if($group_count > $#nanobody_group) { last; } #this isn't a new group, it's the end!
		#if($group_count > 30) { last; } 
		
		if(defined ${$nanobody_group[$group_count]}{-1})
		{#next group doesn't exist (it's part of another group), skip
			@current_group = (); $cur_group_i = 0;
			$total_proteins_in_file--;
			next;
		}
		else
		{
			$proteins_sorted_count = $group_count;
			
			#@current_group = sort { $a <=> $b } keys %{$nanobody_group[$group_count]};
			
			####### sort group members by HT-Seq count ######
			my @unsorted_current_group = keys %{$nanobody_group[$group_count]};
			my @ht_seq_count_in_group = ();
			for(my $k=0; $k <= $#unsorted_current_group; $k++)
			{
				push(@ht_seq_count_in_group, [$proteins_sorted[$unsorted_current_group[$k]][2], $unsorted_current_group[$k]])
			}
			my @ht_seq_count_in_group_sorted = sort { $b->[0] <=> $a->[0] } @ht_seq_count_in_group;
			@current_group = ();
			for(my $k=0; $k <= $#ht_seq_count_in_group_sorted; $k++)
			{
				push(@current_group, $ht_seq_count_in_group_sorted[$k][1]);
			}
			#################################################
			
			
			$num_group_members = $#current_group+2; #add in an extra 1 for the group head
			$group_head = 1; $cur_group_i = -1;
			
			#make a list of all the names in this group, for counting up # group-proteins that match a peptide
			%peptide_match_count = ();
			@names_in_group = ();
			push(@names_in_group, $proteins_sorted[$group_count][7]); #group head
			for(my $k = 0; $k <= $#current_group; $k++)
			{
				push(@names_in_group, $proteins_sorted[$current_group[$k]][7]);
			}	
		}
		
		if($make_new_file)
		{
			close_OUT_ALL($group_div, ++$file_num);
			open_OUT_ALL($file_num);
			$make_new_file = 0;
		}
	}
	else
	{
		$proteins_sorted_count = $current_group[$cur_group_i];
		$group_head = 0;
	}
	
	#read in the name of the next protein in the sorted array
	my @protein_peptides; 
	my $num_pp = 0;
	my $name;
	$name = $proteins_sorted[$proteins_sorted_count][7];
	
	$protein_peptides_unique_count{$name}=0;
	my $temp=$protein_peptides{$name};
	my $temp1 = $protein_peptides_pos{$name};
	while($temp=~s/^#([^#]+)#//)
	{#read in peptides for this protein
		my $pep=$1;
		$protein_peptides[$num_pp][0] = $pep;
		$temp1 =~ s/^#([^#]+)#//;
		$protein_peptides[$num_pp][1] = $1; #the starting position of the peptide in the protein
		$num_pp++
	}
	#sort peptides based on starting position in the protein
	my @sorted_protein_peptides = sort { $a->[1] <=> $b->[1] } @protein_peptides;
	
	#print out the protein with stats
	#format: CDR1: X% (X/X); CDR2: X% (X/X); CDR3: X% (X/X); combined CDR: X% (X/X); overall: X% (X/X)
	my $cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
	my $cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
	my $cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
	my $total_cdr_cover = $protein_cdr1_cover{$name} + $protein_cdr2_cover{$name} + $protein_cdr3_cover{$name};
	if($group_head)
	{
		if($group_div) { print OUT_ALL qq!</div>!; $group_div = 0; }
		print OUT_ALL "<br><b title='rank of the group'>$real_group_count</b>";
		
		#add spaces after rank so all results line up
		my $num_sp = 3 - floor(log10_($real_group_count)+1);
		for(my $s = 0; $s < $num_sp; $s++) { print OUT_ALL ' ';  }
	
		$real_group_count++;
		$real_group_count_on_page++;
		if($#current_group >= 0)
		{
			$group_div_head = $proteins_sorted_count;	
			print OUT_ALL qq!<img src="$img_source_plus" title="click to expand the group for viewing more sequences" onclick="ec('d_g_$proteins_sorted_count', 'i_g_$proteins_sorted_count')" id="i_g_$proteins_sorted_count" style="cursor:hand;" alt="+" />!;
		}
		else
		{
			print OUT_ALL qq!<img src="$img_source_star" style="cursor:hand;" alt="*" disabled="true" title="no more sequences in this group"/>!;
		}
	}
	elsif($cur_group_i == 0)
	{#its the first sequence in group, print out div before it
		print OUT_ALL qq!<div id="d_g_$group_div_head" style="display:none">!;
		$group_div = 1;
	}
	
	if ($group_head)
	{
		if ($num_group_members > $MAX_GROUP_PROTEINS_IN_FILE) { print OUT_ALL "<b title='number of sequences in this group (first $MAX_GROUP_PROTEINS_IN_FILE are shown)'>($num_group_members)</b>"; }
		else { print OUT_ALL "<b title='number of sequences in this group'>($num_group_members)</b>"; }
		
		#add spaces after num group members so all results line up
		my $num_sp = 3 - floor(log10_($num_group_members)+1);
		for(my $s = 0; $s < $num_sp; $s++) { print OUT_ALL ' '; }
	}
	else
	{#print out indentation so group members line up with head of group
		print OUT_ALL '          ';
	}
	print OUT_ALL qq!<img src="$img_source_plus" title="click to expand the sequence for peptide mapping information" onclick="ec('d_$proteins_sorted_count', 'i_$proteins_sorted_count')" id="i_$proteins_sorted_count" style="cursor:hand;" alt="+" />!;
	if ($show_score)
	{
		print OUT_ALL qq!<b title="sequence name">$name</b> <b title="SCORE=8*CDR3-COV + 2*CDR2-COV + 2*CDR1-COV + 2*SEQ-COV + CDR3-LEN/15 + CDR2-LEN/10 + CDR1-LEN/9">Score: $proteins_sorted[$proteins_sorted_count][8];</b> <b title="CDR1 coverage">CDR1: $cdr1_perc\% ($protein_cdr1_cover{$name}/$cdr1_len);</b> <b title="CDR2 coverage">CDR2: $cdr2_perc\% ($protein_cdr2_cover{$name}/$cdr2_len);</b> <b title="CDR3 coverage">CDR3: $proteins_sorted[$proteins_sorted_count][0]\% ($protein_cdr3_cover{$name}/$proteins_sorted[$proteins_sorted_count][1]);</b> <b title="CDR 1, 2 and 3 combined coverage">combined CDR: $proteins_sorted[$proteins_sorted_count][3]\% ($total_cdr_cover/$proteins_sorted[$proteins_sorted_count][4]);</b> <b title="overall sequence coverage">overall: $proteins_sorted[$proteins_sorted_count][5]\% ($protein_total_cover{$name}/$proteins_sorted[$proteins_sorted_count][6]);</b> <b title="number of HT-DNA sequencing reads that produced this sequence">HT-seq count: $proteins_sorted[$proteins_sorted_count][2]</b>\n             !;
	}
	else
	{
		print OUT_ALL qq!<b title="sequence name">$name</b> <b title="CDR1 coverage">CDR1: $cdr1_perc\% ($protein_cdr1_cover{$name}/$cdr1_len);</b> <b title="CDR2 coverage">CDR2: $cdr2_perc\% ($protein_cdr2_cover{$name}/$cdr2_len);</b> <b title="CDR3 coverage">CDR3: $proteins_sorted[$proteins_sorted_count][0]\% ($protein_cdr3_cover{$name}/$proteins_sorted[$proteins_sorted_count][1]);</b> <b title="CDR 1, 2 and 3 combined coverage">combined CDR: $proteins_sorted[$proteins_sorted_count][3]\% ($total_cdr_cover/$proteins_sorted[$proteins_sorted_count][4]);</b> <b title="overall sequence coverage">overall: $proteins_sorted[$proteins_sorted_count][5]\% ($protein_total_cover{$name}/$proteins_sorted[$proteins_sorted_count][6]);</b> <b title="number of HT-DNA sequencing reads that produced this sequence">HT-seq count: $proteins_sorted[$proteins_sorted_count][2]</b>\n             !;
	}
	
	#underline AA's that correspond to identified peptides:
	#show the 3 CDR regions in red
	my $cov_by_aa_ref = $protein_total_coverage_by_aa{$name};
	my $sequence = $proteins{$name};
	my $seq_len = length($sequence);
	
	if ($FILTER_FIX_SEQ_END && !$use_primers)
	{ $sequence =~ s/X+$//; }
	
	if ($use_tail)
	{ $sequence .= $protein_max_cover_tail{$name}; }
	
	my @seq = split('', $sequence);
	my $u_on = 0;
	my $cdr_i = 1;
	
	for(my $i = 0; $i <= $#seq; $i++)
	{
		my $pre = ""; my $post = "";
		
		#bold primers if using
		if ($use_primers && ($i == 0 || $i == ($seq_len-$P2_SEQ_LENGTH))) { $pre .= "<b>"; }
		if ($use_primers && ($i == ($seq_len-1) || $i == ($P1_SEQ_LENGTH-1))) { $post .= "</b>"; }
		
		if(($i == $CDR1{$name}[0]) || ($i == $CDR2{$name}[0]) || ($i == $CDR3{$name}[0]))
		{ $pre .= "<font title=\"CDR $cdr_i\" color=\"\#FF0000\">"; $cdr_i++;}
		if(($i == $CDR1{$name}[1]) || ($i == $CDR2{$name}[1]) || ($i == $CDR3{$name}[1]))
		{ $post .= "</font>"; }
		
		if(defined $$cov_by_aa_ref{$i}) { if(!$u_on) { $pre .= "<u title='peptide coverage'>"; $u_on = 1; } }
		elsif($u_on) { $pre .= "</u>"; $u_on = 0; }
		
		print OUT_ALL $pre.$seq[$i].$post;	
	}
	if($u_on) { print OUT_ALL "</u>"; }	
	print OUT_ALL "\n"; 

	#print out all sorted peptides	
	print OUT_ALL qq!<div id="d_$proteins_sorted_count" style="display:none">!; 
	for(my $i = 0; $i <= $#sorted_protein_peptides; $i++) 
	{
		my $expect = $peptides{$sorted_protein_peptides[$i][0]}[0];
		my $orig = $peptides{$sorted_protein_peptides[$i][0]}[1]; 
		
		if($i > 0) { print OUT_ALL "\n"; }
		printf OUT_ALL "   <b title='the best score (log(e)) for this peptide'>%-9s</b> ", $expect;
		
		for(my $s = 0; $s < $sorted_protein_peptides[$i][1]; $s++) { print OUT_ALL " "; }
		
		my $source = $peptides{$sorted_protein_peptides[$i][0]}[2];
		my $num_found = $peptides{$sorted_protein_peptides[$i][0]}[3];
		my $protein_uid = $peptides{$sorted_protein_peptides[$i][0]}[4]; 
		my $domain_id = $peptides{$sorted_protein_peptides[$i][0]}[5]; 
		my $tandem_xml_file = $peptides{$sorted_protein_peptides[$i][0]}[6];
		my $num_proteins = $peptide_proteins_count{$sorted_protein_peptides[$i][0]};
		my $protein_names = $peptide_proteins{$sorted_protein_peptides[$i][0]}; #could make a hash out of this
		my $all_protein_names = 1;
		
		#for each peptide, find out how many proteins from this group it matches
		my $num_group_matches = 0;
		if (defined $peptide_match_count{$sorted_protein_peptides[$i][0]})
		{
			$num_group_matches = $peptide_match_count{$sorted_protein_peptides[$i][0]};
		}
		else
		{
			for(my $k = 0; $k <= $#names_in_group; $k++)
			{
				if ($peptide_proteins_lookup{$sorted_protein_peptides[$i][0]}{$names_in_group[$k]}) 
				#if (index($protein_names, $names[$k] . '<br>') != -1)
				{
					$num_group_matches++;
				}
			}
			$peptide_match_count{$sorted_protein_peptides[$i][0]} = $num_group_matches;	
		}
		
		#score:
		my $score = sprintf("%.2f", $num_proteins/$num_group_members);
		my $score2;
		if ($num_group_matches/$num_proteins == 1) {
			$score2 = "100%";
		}
		else
		{
			$score2 = sprintf("%.2f", $num_group_matches/$num_proteins * 100);
			$score2 .= "%";
		}
		
		print OUT_ALL qq!<a title="click here to view best matching spectrum" href="/llama-magic-cgi/peptide.pl?path=$tandem_xml_file&uid=$protein_uid&id=$domain_id&label=$orig">$orig</a>!;
		print OUT_ALL " <b title='uniqueness score (# matching group members/# total matching sequences)'>($score2) ($num_group_matches / $num_proteins)</b>  <b title='spectra count'>($num_found)</b>";
		print OUT_ALL qq!<img src="$img_source_plus" title="click to expand for spectra information" onclick="ec('d_$proteins_sorted_count-$i', 'i_$proteins_sorted_count-$i')" id="i_$proteins_sorted_count-$i" style="cursor:hand;" alt="+" />!;
		print OUT_ALL qq!<div id="d_$proteins_sorted_count-$i" style="display:none">$source</div>!;
		if($num_proteins < 100)
		{
			print OUT_ALL "  <b title='number of sequences matched for this peptide across entire database'>($num_proteins)</b>";
			print OUT_ALL qq!<img src="$img_source_plus" title="click to expand for all matched sequences" onclick="ec('d_$proteins_sorted_count-$i-2', 'i_$proteins_sorted_count-$i-2')" id="i_$proteins_sorted_count-$i-2" style="cursor:hand;" alt="+" />!;
			print OUT_ALL qq!<div id="d_$proteins_sorted_count-$i-2" style="display:none">$protein_names</div>!;
		}
		else
		{
			print OUT_ALL "  <b title='number of sequences matched for this peptide across entire database'>($num_proteins)</b>";
			print OUT_ALL qq!<img src="$img_source_star" style="cursor:hand;" alt="*" disabled="true" title="100 or more sequences matched this peptide"/>!;
		}	
	}
	print OUT_ALL qq!</div>!; 
	
}

close_OUT_ALL($group_div, 1); 

close(LOG);


##############subroutines##########################################################################################
sub filter_input
{#trim until the first M and check if length is >= 107
 #also check if 'QVT' is is at the end (last 10 AA's of sequence)
	my $seq = $_[0];
	my $m_i = index($seq, 'M');
	if($m_i >= 0) { $_[0] = substr($seq, $m_i); } 
	else { return 0; }
	
	if(length($_[0]) >= $MIN_SEQ_LENGTH)
	{
		#length okay, 'M' found, check 'QVT' in the last 10 characters
		
		if($ADD_MAQ_TO_START_CHANGE_QVT_CHECK)
		{
			#DON'T CHECK AT ALL
			if($QVT_LEN_FROM_END == -1) { return 1; }
			
			#CHECK OLD WAY, WITH LENGTH FROM END FOR QVT expanded
			my $end_seq = substr($_[0], length($_[0])-$QVT_LEN_FROM_END); 
			if($end_seq =~ /QVT/)
			{
				return 1;
			}
			return 0;	
		}
		else
		{
			if ($FILTER_FIX_SEQ_END)
			{
				#QVT is found, which is part of primer
				#in order for the cdr3-finding to work consistently for different length primers,
				#put a standard amount of AA's at the end, instead of allowing differing lengths after the QVT
				#put X's at the end so it will not match to any peptides
				my $pos = index($_[0], 'QVT', length($_[0])-10);
				if ($pos != -1)
				{
					my $clipped_seq = substr($_[0], 0, length($_[0])-(length($_[0])-$pos));
					my $end_part = substr($_[0], $pos);
					my $l = length($end_part);
					while ($l < ($NUM_AA_AFTER_QVT + 3)) # length(QVT) == 3
					{
						$end_part .= 'X';
						$l++;
					}
					$_[0] = $clipped_seq . $end_part;
					return 1;
				}
				#print "Filter: No QVT.\n";
				return 0;
			}
			else
			{
				#OLD WAY
				my $end_seq = substr($_[0], length($_[0])-10); 
				if($end_seq =~ /QVT/)
				{
					return 1;
				}
				return 0;	
			}
		}
	}
	else { return 0;  }
	
	return 0; 
}

sub find_cdr1
{#finds the cdr1 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	my $left_area_start=20; my $left_area_end=26;
	my $right_area_start=32; my $right_area_end=40;
	my $default_left_area=28; my $default_right_area=36;
	my $default_cdr_len = 8;
	
	if($new_primers)
	{#shift for new primers to the right
		$left_area_start+=9; $left_area_end+=9;
		$right_area_start+=9; $right_area_end+=9;
		$default_left_area+=9; $default_right_area+=9;
	}
	
	#STARTING POS. OF CDR1:
	my $left_area = substr($seq, $left_area_start, $left_area_end-$left_area_start+1); #look from pos. 20 - 26 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if(($la_i = index($left_area, 'SC')) < 0)
	{#didn't find 'SC', look for 'C'
		$la_i = index($left_area, 'C');
	} else { $la_i++; } #'C' is our marker, so advance past 'S'
	if($la_i >= 0) { $left_cdr = $la_i + $left_area_start + 5; } #CDR1 starts at 'C' + 5 (add 20 to put it back in the full sequence)
	
	#ENDING POS. OF CDR1:
	my $right_area = substr($seq, $right_area_start, $right_area_end-$right_area_start+1); #look from pos. 32 - 40 of seq (0-based)
	my $ra_i; my $right_cdr = -1;
	if($right_area =~ /(W.R)/)
	{#if we found 'WXR', find its index
		$ra_i = index($right_area, $1);
	}
	else { $ra_i = index($right_area, 'W'); } #didn't find 'WXR', look for 'W'
	if($ra_i >= 0) { $right_cdr = $ra_i + $right_area_start - 1; } #CDR1 ends at 'W' - 1 (add 32 to put it back in the full sequence)
	
	#check if st/end found and if not follow rules:
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = $default_left_area; $right_cdr = $default_right_area; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - $default_cdr_len; }
	elsif($right_cdr == -1) { $right_cdr = $left_cdr + $default_cdr_len; }
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	return 1;
}

sub find_cdr2
{#finds the cdr2 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	my $left_area_start=32; my $left_area_end=40;
	my $right_area_start=63; my $right_area_end=72;
	my $default_left_area=51; my $default_right_area=60;
	my $default_cdr_len = 9;
	
	if($new_primers)
	{#shift for new primers to the right
		$left_area_start+=9; $left_area_end+=9;
		$right_area_start+=9; $right_area_end+=9;
		$default_left_area+=9; $default_right_area+=9;
	}
	
	#STARTING POS. OF CDR2:
	my $left_area = substr($seq, $left_area_start, $left_area_end-$left_area_start+1); #look from pos. 32 - 40 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if($left_area =~ /(W.R)/)
	{#if we found 'WXR', find its index
		$la_i = index($left_area, $1);
	}
	else { $la_i = index($left_area, 'W'); } #didn't find 'WXR', look for 'W'
	if($la_i >= 0) { $left_cdr = $la_i + $left_area_start + 14; } #CDR2 starts at 'W' + 14 (add 32 to put it back in the full sequence)
	
	#ENDING POS. OF CDR2:
	my $right_area = substr($seq, $right_area_start, $right_area_end-$right_area_start+1); #look from pos. 63 - 72 of seq (0-based)  
	my $ra_i; my $right_cdr = -1;
	$ra_i = index($right_area, 'RF');
	if($ra_i >= 0) { $right_cdr = $ra_i + $right_area_start - 8; } #CDR2 ends at 'R' - 8 (add 63 to put it back in the full sequence)
	
	#check if st/end found and if not follow rules:
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = $default_left_area; $right_cdr = $default_right_area; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - $default_cdr_len; }
	elsif($right_cdr == -1) { $right_cdr = $left_cdr + $default_cdr_len; }
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	return 1;
}

sub find_cdr3_new
{
	my $seq = $_[0];
	my $name = $_[4];
	my $n = length($seq)-1;
	
	#below params only used as a fall back if the new cdr3 finding doesn't work
	my $left_area_start=92; my $left_area_end=102;
	my $right_area_start; my $right_area_end;
	my $default_left_area; my $default_right_area;
	my $default_cdr_len = 11;
	
	if ($use_primers)
	{#if primers were located, p2 is placed in sequence, but the random 12-mer portion following p2 is removed
	 #so CDR3-finding must be adjusted - reduced by 4 AA's (corresponds to random 12-mer of nucleotide seq)
		#below is amount from end
		$right_area_start = 10;
		$right_area_end = 0;
		
		$default_left_area = 17;
		$default_right_area = 6;
		if($new_primers)
		{
				# these 2 are from start of string
				$left_area_start+=9; $left_area_end+=9;

				# these are from end of string
				my $add_amt = 10;
				$right_area_start+=$add_amt;
				$right_area_end+=$add_amt;
				$default_left_area+=$add_amt;
				$default_right_area+=$add_amt;
		}
	}
	else
	{
		#below is amount from end
		$right_area_start = 14;
		$right_area_end = 4;
		
		$default_left_area = 21;
		$default_right_area = 10;
	}
	
	#STARTING POS. OF CDR3:
	my $end_cdr2 = $_[1];
	my $offset = $end_cdr2+1;
	my $left_area = substr($seq, $offset); #look from end of cdr2 until end of seq
	
	my $la_i=-1; my $left_cdr = -1;
	if($left_area =~ /(DTAVY.C)/)
	{
		$la_i = index($left_area, $1);
		$la_i += 6; #'C' is our marker, so advance past 'DXXXYX'
	}
	elsif($left_area =~ /(D.AVY.C)/)
	{
		$la_i = index($left_area, $1);
		$la_i += 6; #'C' is our marker, so advance past 'DXXXYX'
	}
	elsif($left_area =~ /(DT.VY.C)/)
	{
		$la_i = index($left_area, $1);
		$la_i += 6; #'C' is our marker, so advance past 'DXXXYX'
	}
	elsif($left_area =~ /(DTA.Y.C)/)
	{
		$la_i = index($left_area, $1);
		$la_i += 6; #'C' is our marker, so advance past 'DXXXYX'
	}
	elsif($left_area =~ /(D...Y.C)/)
	{
		$la_i = index($left_area, $1);
		$la_i += 6; #'C' is our marker, so advance past 'DXXXYX'
	}
	else
	{#incorrect CDR2 finding leading to incorrect CDR3 finding
		#fall back to old way
		$offset = $left_area_start;
		$left_area = substr($seq, $offset, $left_area_end-$left_area_start+1); #look from pos. 92 - 102 of seq (0-based)
		if($left_area =~ /(Y.C)/)
		{#if we found 'YXC', find its index
			$la_i = index($left_area, $1);
			$la_i += 2; #'C' is our marker, so advance past 'YX'
		}
		else { $la_i = index($left_area, 'C'); } #didn't find 'YXC', look for 'C'	
		
	}
	if($la_i >= 0) { $left_cdr = $la_i + $offset + 3; } #CDR3 starts at 'C' + 3 (add 92 to put it back in the full sequence)
	
	#ENDING POS. OF CDR3:
	if($left_cdr == -1) { $offset = $end_cdr2 + 1 } #look from end of cdr2 until end of seq
	else { $offset = $left_cdr+1; } # look from cdr3 start until end of seq
	my $right_area = substr($seq, $offset); 
	my $ra_i=-1; my $right_cdr = -1;
	my $subtract_amount = 4;
	if ($right_area =~ /(GTQVTV)/) 
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(.TQVTV)/)
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(G.QVTV)/)
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(GT.VTV)/)
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(GTQ.TV)/)
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(GTQV.V)/)
	{
		$ra_i = index($right_area, $1);
	}
	elsif($right_area =~ /(GTQVT.)/)
	{
		$ra_i = index($right_area, $1);
	}
	else
	{
		if($ADD_MAQ_TO_START_CHANGE_QVT_CHECK)
		{#don't want to look from end (old way), since there's different amounts of AA's at the end...
		 #try with only QVT
			if($right_area =~ /(..QVT.)/)
			{
				$ra_i = index($right_area, $1);
			}
			elsif($right_area =~ /(WGQ)/)
			{
				$ra_i = index($right_area, $1);
				$subtract_amount=1;
			}
			elsif($right_area =~ /(WG.)/)
			{
				$ra_i = index($right_area, $1);
				$subtract_amount=1;
			}
			elsif($right_area =~ /(.GQ)/)
			{
				$ra_i = index($right_area, $1);
				$subtract_amount=1;
			}
		}
		else
		{
			#fall back to old way
			$offset = $n-$right_area_start;
			$right_area = substr($seq, $offset, $right_area_start-$right_area_end+1); #look from pos. n-14 - n-4 of seq (n = last index of seq)
			$subtract_amount = 1;
			if(($ra_i = index($right_area, 'WGQ')) < 0)
			{#didn't find 'WGQ', look for 'WG'
				if(($ra_i = index($right_area, 'WG')) < 0)
				{#didn't find 'WG', look for 'W'
					if(($ra_i = index($right_area, 'W')) < 0)
					{#didn't find 'W', look for 'GQ'
						$ra_i = index($right_area, 'GQ');
						if($ra_i >= 0) { $subtract_amount = 2; } #if 'GQ' found, CDR3 ends at 'G' - 2 
					}
				}
			}
		}
	}
	
	if($ra_i >= 0) { $right_cdr = $ra_i + $offset - $subtract_amount; } 
	
	#check if st/end found and if not follow rules: 
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = $n-$default_left_area; $right_cdr = $n-$default_right_area; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - $default_cdr_len; }
	elsif($right_cdr == -1) { $right_cdr = ($left_cdr + $default_cdr_len) <= $n ? ($left_cdr + $default_cdr_len) : $n; }
	
	if($left_cdr > $right_cdr)
	{
		#print "Alert! left_cdr ($left_cdr) > right_cdr ($right_cdr) in find_cdr3: $seq\n";
		$left_cdr = $n-1;
		$right_cdr = $n;
	} 
	
	$_[2] = $left_cdr;
	$_[3] = $right_cdr;
	
	return 1;
}

sub find_cdr3
{#finds the cdr3 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	my $left_area_start=92; my $left_area_end=102;
	my $right_area_start; my $right_area_end;
	my $default_left_area; my $default_right_area;
	my $default_cdr_len = 11;
	if ($use_primers)
	{#if primers were located, p2 is placed in sequence, but the random 12-mer portion following p2 is removed
	 #so CDR3-finding must be adjusted - reduced by 4 AA's (corresponds to random 12-mer of nucleotide seq)
		#below is amount from end
		$right_area_start = 10;
		$right_area_end = 0;
		
		$default_left_area = 17;
		$default_right_area = 6;
	}
	else
	{
		#below is amount from end
		$right_area_start = 14;
		$right_area_end = 4;
		
		$default_left_area = 21;
		$default_right_area = 10;
	}
	
	#STARTING POS. OF CDR3:
	my $left_area = substr($seq, $left_area_start, $left_area_end-$left_area_start+1); #look from pos. 92 - 102 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if($left_area =~ /(Y.C)/)
	{#if we found 'YXC', find its index
		$la_i = index($left_area, $1);
		$la_i += 2; #'C' is our marker, so advance past 'YX'
	}
	else { $la_i = index($left_area, 'C'); } #didn't find 'YXC', look for 'C'
	if($la_i >= 0) { $left_cdr = $la_i + $left_area_start + 3; } #CDR3 starts at 'C' + 3 (add 92 to put it back in the full sequence)
	
	#ENDING POS. OF CDR3:
	my $n = length($seq)-1; my $n1 = $n-$right_area_start;
	my $subtract_amount = 1; 
	my $right_area = substr($seq, $n1, $right_area_start-$right_area_end+1); #look from pos. n-14 - n-4 of seq (n = last index of seq)
	my $ra_i; my $right_cdr = -1;
	if(($ra_i = index($right_area, 'WGQ')) < 0)
	{#didn't find 'WGQ', look for 'WG'
		if(($ra_i = index($right_area, 'WG')) < 0)
		{#didn't find 'WG', look for 'W'
			if(($ra_i = index($right_area, 'W')) < 0)
			{#didn't find 'W', look for 'GQ'
				$ra_i = index($right_area, 'GQ');
				if($ra_i >= 0) { $subtract_amount = 2; } #if 'GQ' found, CDR3 ends at 'G' - 2 
			}
		}
	} 
	if($ra_i >= 0) { $right_cdr = $ra_i + $n1 - $subtract_amount; } #CDR3 ends at 'W' - 1 (or 'G' - 2) (add n-14 to put it back in the full sequence)
	
	#check if st/end found and if not follow rules: 
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = $n-$default_left_area; $right_cdr = $n-$default_right_area; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - $default_cdr_len; }
	elsif($right_cdr == -1) { $right_cdr = ($left_cdr + $default_cdr_len) <= $n ? ($left_cdr + $default_cdr_len) : $n; }
	
	if($left_cdr > $right_cdr)
	{
		#print "Alert! left_cdr ($left_cdr) > right_cdr ($right_cdr) in find_cdr3: $seq\n";
		$left_cdr = $n-1;
		$right_cdr = $n;
	} 
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	
	return 1;
}

sub compare_cdrs
{#return true if the 2 sequences differ by 0 or 1 AA
	my @cdr1 = split '', $_[0];
	my @cdr2 = split '', $_[1];
	
	if(abs($#cdr2-$#cdr1) > 1)  { return 0; }
	
	my $num_diff = $#cdr2 > $#cdr1 ? ($#cdr2-$#cdr1) : 0;
	
	my $i = 0;
	foreach (@cdr1)
	{
		if(($i > $#cdr2) || ($cdr2[$i] ne $_)) { $num_diff++; }
		if($num_diff > 1) { return 0; }
		$i++;
	}
	
	#if($num_diff > 1) { return 0; }
	#else { return 1; }
	return 1;
}

sub consolidate_groups
{#adds group $j to group $i, if a matching sequence in both groups is found ($nanobody_group[$i] should be defined)
 #then consolidates all previous groups (back to $i) (in case adding group $j caused additional matches w/ previous groups)
	my $i = shift;
	my $j = shift;
	
	#if($nanobody_group[$i])
	if($j <= $i) { return 0; }
	
	my $consolidate = 0;
	if(!defined ${$nanobody_group[$j]}{-1} && !defined ${$nanobody_group[$i]}{-1}) #make sure these groups weren't removed already
	{
		if(!defined ${$nanobody_group[$i]}{$j})
		{
			foreach (keys %{$nanobody_group[$j]})
			{
				if(defined ${$nanobody_group[$i]}{$_})
				{ $consolidate = 1; last; }
			}
		}
		else { $consolidate = 1; }
	}
	
	if($consolidate) 
	{#put $j and all in its group into the $i group
	 #remove $j group
		${$nanobody_group[$i]}{$j} = 1;
		foreach (keys %{$nanobody_group[$j]})
		{
			${$nanobody_group[$i]}{$_} = 1;
		}
		%{$nanobody_group[$j]} = ();
		${$nanobody_group[$j]}{-1} = 1; #mark this group as removed (different than a sequence that matches no other sequences)
		
		
	}
	consolidate_groups($i, $j-1);
}

sub close_OUT_ALL
{
	my $group_div = shift;
	my $next_num = shift;
	
	if($group_div) { print OUT_ALL qq!</div>!; }
	print OUT_ALL "</pre>";
	
	if ($next_num > 1)
	{
		my $cur_num = $next_num-1;
		#put link to go to next page
		print OUT_ALL qq(<br><br>page $cur_num | <a href="$fileroot.$next_num.cdr_coverage.html">next page</a>);
	}
	
	close(OUT_ALL);
}

sub open_OUT_ALL
{
	my $num = shift;
	
	#if($num == 1) { open(OUT_ALL,">$filedir/$fileroot.cdr_coverage.html") || die "Could not open for writing: $filedir/$fileroot.cdr_coverage.html\n"; }
	#else { open(OUT_ALL,">$filedir/$fileroot.$num.cdr_coverage.html") || die "Could not open for writing: $filedir/$fileroot.$num.cdr_coverage.html\n"; }
	if(!open(OUT_ALL,">$filedir/$fileroot.$num.cdr_coverage.html"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.$num.cdr_coverage.html ($!)\n";
		close(LOG);
		exit(1);
	}
	
	print OUT_ALL <<JSCRIPT;
	<SCRIPT LANGUAGE="JavaScript">
	// Row Hide function.
	function ec(tId, clickIcon)
	{
		dstyle = document.getElementById(tId).style.display;
		if (dstyle == "none")
		{
			document.getElementById(tId).style.display = "";
			document.getElementById(clickIcon).src = "$img_source_minus";
			document.getElementById(clickIcon).alt = "-";
		}
		else
		{
			document.getElementById(tId).style.display = "none";
			document.getElementById(clickIcon).src = "$img_source_plus";
			document.getElementById(clickIcon).alt = "+";
		}
	}
	</SCRIPT>
JSCRIPT
	print OUT_ALL "<pre>\n";
}

sub strip_primer_info
{
	my $name = shift;
	my $passed = 0;
	if($name =~ s/(_fr[012])_fwd_(\w*)_rev_(\w*)$//)
	{
		my $fr = $1;
		my $p1 = $2;
		my $p2 = $3;
		my $x1 = '';
		my $x2 = '';
		if ($p1 ne '' && $p2 ne '')
		{
			$passed = 1;
			while ($P1_SEQ_LENGTH-(length($p1)+length($x1)) > 0) { $x1 .= 'X'; }
			while ($P2_SEQ_LENGTH-(length($p2)+length($x2)) > 0) { $x2 .= 'X'; }
		}
		
		return ($passed, $name.$fr, $x1, $x2);
	}
	else { return ($passed, $name, '', ''); }
	
}


sub map_peptides
{
	my $seq = $_[0];
	my $protein_cov_ref = $_[1];
	
	$_[2] = '';
	$_[3] = 0;
	$_[4] = '';
	my $peptides_ref = $_[5];
	
	my $pep_count = 0;
	foreach my $pep (keys %{$peptides_ref}) #%peptides)
	{
		#$pep_count++;
		#if ($pep_count > 40) {
		#	last;
		#}
		
		my $pos = 0; 
		while(($pos = index($seq, $pep, $pos)) >= 0)
		{
			my $char = '';
			if($pos > 0) { $char = substr($seq, $pos-1, 1); }
			if(!$RK_FIXATION || ($char eq 'K' || $char eq 'R'))
			#$RK_FIXATION: not allowing peptides at beginning of sequence, MUST have R/K before peptide (simplification)
			{
				$_[2] .= "#$pep#";
				$_[3]++;
				$_[4] .= "#$pos#";
				
				#tally peptide coverage for this protein
				my $pep_len = length($pep);
				for(my $i = $pos; $i < ($pos+$pep_len); $i++) { $$protein_cov_ref{$i} = 1; }
			}
			$pos = $pos + length($pep);
			
		}
	}
}

sub read_in_xtandem_results
{
	my $input_filename = shift;
	my $tandem_filename = shift;
	if(!open(IN,"$input_filename"))
	{
		print LOG "Error: Could not open '$input_filename' ($!)\n";
		close(LOG);
		exit(1);
	}
	
	my $cur_count_pep;
	my $cur_count_pep_unique;
	my $line=<IN>;
	while($line=<IN>)
	{
		$line =~ s/\r\n$//;
		$line =~ s/\n$//;
		#chomp($line);
		
		if ((!$stand_alone && $line=~/^([A-Za-z]+)\t([0-9\-\+edED\.]*)\t([^\t]+)\t([^\t]+)\t([^\t]*)/) ||
			#use this if we are running through llama magic
			($stand_alone && $line=~/^([A-Za-z]+)\t([0-9\-\+edED\.]*)/))
			#use this if we don't have protein_uid and domain_id, and spectrum fields
		{
			my $original_pep=$1;
			my $expect=$2;
			
			my $protein_uid;
			my $domain_id;
			my $spectrum;
			
			if(defined $3) { $protein_uid = $3; }
			else { $protein_uid = ""; }
			if(defined $4) { $domain_id = $4; }
			else { $domain_id = ""; }
			if(defined $5) { $spectrum = $5; }
			else { $spectrum = ""; }
		
			my $pep = $original_pep;
			$pep=~tr/L/I/;
			if (!(defined $peptides{"\U$pep"})) 
			{#if peptide found more than once, only the lowest expectation value 
			 #and the associated original sequence is stored
				$cur_count_pep_unique++; 
				if($spectrum) { $peptides{"\U$pep"}=[$expect,$original_pep,"$spectrum ($expect)", 1, $protein_uid, $domain_id, $tandem_filename]; }
				else { $peptides{"\U$pep"}=[$expect,$original_pep,"", 1, $protein_uid, $domain_id, $tandem_filename]; }
				
			}
			else
			{
				my $prev_expect = $peptides{"\U$pep"}[0];
				if($prev_expect > $expect)
				{
					$peptides{"\U$pep"}[0] = $expect;
					$peptides{"\U$pep"}[4] = $protein_uid;
					$peptides{"\U$pep"}[5] = $domain_id;
				}
				if($spectrum)
				{
					if($peptides{"\U$pep"}[2]) { $peptides{"\U$pep"}[2] .= "<br><br>$spectrum ($expect)"; }
					else { $peptides{"\U$pep"}[2] = "$spectrum ($expect)"; }
				}
				
				$peptides{"\U$pep"}[3]++;
			}
			$cur_count_pep++;
		} 
		else { if ($line=~/\w/) { print LOG qq!Error: $line\n!; } }	
	}
	close(IN);
	print LOG qq!$input_filename\n$cur_count_pep peptides ($cur_count_pep_unique)\n!;
}

sub read_in_protein_db
{
	my $protein_db_filename = shift;

	if (!open (IN,"$protein_db_filename"))
	{
		print LOG "Error: Could not open '$protein_db_filename' ($!)\n";
		close(LOG);
		exit(1);
	}
	
	print LOG "$protein_db_filename\n";
	
	my $name="";
	my $gene_count;
	my $sequence="";
	my $count_seq=0;
	my $all_length_count_seq=0;
	my $line;
	while ($line=<IN>)
	{
		my $name_; my $gene_count_ = 0;
		$line =~ s/\r\n$//;
		$line =~ s/\n$//;
		#chomp($line);
		
		#print "LINE= '$line'\n";
		
		if ($line=~/^>(\S+) \+ (\d+) other/)
		{
			$name_=$1;
			$gene_count_ = $2 + 1;
			
		}
		elsif($line=~/^>(\S+)/)
		{
			$name_=$1;
			$gene_count_ = 1;
		}
		
		if($gene_count_ > 0)
		{#its a name line
			#print "$name $sequence\n";
			if ($name=~/\w/ and $sequence=~/\w/)
			{
				my $protein_passed = 0; my $x1; my $x2;
				if ($use_primers)
				{#strip primers from name, if both were found, then protein passes filter
					
					my $primers_found = 0;
					($primers_found, $name,$x1,$x2) = strip_primer_info($name);
					if ($primers_found) #only use seq len as a filter, since both primers found
					{
						#if primers were cut off, add Xs at beg/end (makes all primers same length)
						$sequence = $x1 . $sequence . $x2;
						if(length($sequence) > $MIN_SEQ_LENGTH) { $protein_passed = 1; }
						
					}
					
				}
				else
				{
					if($ADD_MAQ_TO_START_CHANGE_QVT_CHECK)
					{
						$sequence = 'MAQ' . $sequence;
					}
					#check length from first M is atleast 107 and trim sequence, also check for QVT
					if(filter_input($sequence)) { $protein_passed = 1; }
					else { print "protein failed.\n"; }
				}
				
				if($protein_passed)
				{
					$proteins{$name}=$sequence;
					$protein_gene_counts{$name} = $gene_count;
					$count_seq++;
					$all_length_count_seq++;
				}
				else { $all_length_count_seq++; }
			}
			$name=$name_;
			$gene_count = $gene_count_;
			$sequence="";
		}
		else
		{#its a sequence line
			$sequence.="\U$line";
		}
		#if($count_seq > 10000) { last; }
	}	
	if ($name=~/\w/ and $sequence=~/\w/)
	{
		my $protein_passed = 0; my $x1; my $x2;
		if ($use_primers)
		{#strip primers from name, is both were found, then protein passes filter
			my $primers_found = 0;
			($primers_found, $name,$x1,$x2) = strip_primer_info($name);
			if ($primers_found) #only use seq len as a filter, since both primers found
			{
				#if primers were cut off, add Xs at beg/end (makes all primers same length)
				$sequence = $x1 . $sequence . $x2;
				if(length($sequence) > $MIN_SEQ_LENGTH) { $protein_passed = 1; }
				
			}
		}
		else
		{
			if($ADD_MAQ_TO_START_CHANGE_QVT_CHECK)
			{
				$sequence = 'MAQ' . $sequence;
			}
			#check length from first M is atleast 107 and trim sequence, also check for QVT
			if(filter_input($sequence)) { $protein_passed = 1; }  
		}
		 
		if($protein_passed)
		{
			$proteins{$name}=$sequence;
			$protein_gene_counts{$name} = $gene_count;
			$count_seq++;
			$all_length_count_seq++;
		}
		else 
		{ $all_length_count_seq++; }
	}
	close(IN);
	print LOG "Found $count_seq out of $all_length_count_seq sequences with minimum length.\n";
	return $count_seq;
}

sub read_in_peptide_index
{
	my $index_filename = shift;
	if ($index_filename)
	{
        if(!open(IND_IN,"$index_filename"))
		{
			print LOG "Error: Could not open '$index_filename' ($!)\n";
			print LOG "Index will not be used.\n";
			return 0;
		}
		else
		{
			#read in each name and set of peptides
			#for each set of peptides, look in XTandem found peptides and if present add to index for that protein
			my $name;
			while ($name=<IND_IN>)
			{
				$name =~ s/^>//;
				
				$name =~ s/\r\n$//;
				$name =~ s/\n$//;
				#chomp($name);
				
				my $peps_list = <IND_IN>;
				while($peps_list=~s/^#([^#]+)#//)
				{#read in peptides for this protein
					my $cur_pep = $1;
					$cur_pep=~tr/L/I/;
					if (defined $peptides{$cur_pep})
					{
						$peptide_index{$name}{$cur_pep} = 1;
					}
					
				}
			}
			return 1;
		}
    }
	else
	{
		return 0;
	}
}

