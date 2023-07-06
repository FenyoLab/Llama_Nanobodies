#!/usr/local/bin/perl 

use warnings;
use strict;

my %primers = ('p2'=>'TGGAGACGGTGACCTGGGT',
               'p1'=>'ATGGCTGAGGTGCAGCTGGTGGAGTCTGG');

my %primers_rc = ('p2'=>'ACCCAGGTCACCGTCTCCA',  
                  'p1'=>'CCAGACTCCACCAGCTGCACCTCAGCCAT');
               

my %p1_alt_nuc = (6 => 'C', 8 => 'T');
my %p1_alt_nuc_rc = (20 => 'A', 22 => 'G');

my $max_p1_allowed_mismatches = 8;
my $max_p2_allowed_mismatches = 4;

my $num_sequences_with_primers = 0;
my $num_sequences_with_primers1 = 0;
my $num_sequences_with_primers2 = 0;
my $num_sequences = 0;

my $dir="";
my $out_dir = "";

if ($ARGV[0] && $ARGV[0]=~/\w/) { $dir="$ARGV[0]";}
else { $dir = ""; }

if ($ARGV[1] && $ARGV[1]=~/\w/) { $out_dir="$ARGV[1]";}
else { $out_dir = ""; }

print "Executing remove_primers.pl\n";

if (!opendir(DIR,"$dir")) { print "Error reading $dir ($!)\n"; exit(2); }

my @allfiles=readdir DIR;
closedir DIR;
my $count_fasta = 0;

my $outfile = "$out_dir/dna_removed_primers.fasta";
if(!open(OUT, ">$outfile")) { print "Error: Could not open $outfile for writing.\n"; exit(1); }
		
foreach my $filename (@allfiles)
{#for each fasta file  or fastq file in the directory:
	my $ftype;
	if ($filename =~ /\.fas?t?a?$/i) { $ftype = 'FASTA'; }
	elsif($filename =~ /\.fastq$/i || $filename =~ /\.fq$/i) { $ftype = 'FASTQ'; }
	else { next; }
	
	if (!open (IN,"$dir/$filename")) { print "Error opening $dir/$filename ($!).\n"; next; }
	print "Opened $dir/$filename\n";
	$count_fasta++;
	
	my $name="";
	my $description="";
	my $sequence="";
	my $line="";
	
	while ($line=<IN>)
	{
		chomp($line);
		if ($line=~/^[>@](\S+)\s?(.*)$/)
		{
			my $name_=$1; 
			my $description_=$2;
			$name_ =~ s/[\:\-\\\/]/_/g;
			if ($name=~/\w/ and $sequence=~/\w/)
			{#the entire sequence has been read in
				#look for primers in the sequence, record p-start-type, p-start-seq, p-end-type-p-end-seq
				my @ret_vals = remove_primers($sequence);
				
				#save the new sequence to a file, in the description indicate primers removed
				print OUT ">$name";
				if ($ret_vals[0]) { print OUT " $ret_vals[3]=$ret_vals[5]"; }
				else { print OUT " X"; }
				if ($ret_vals[1]) { print OUT " $ret_vals[4]=$ret_vals[6]"; }
				else { print OUT " X"; }
				print OUT "\n$ret_vals[2]\n";
				
				if ($ret_vals[0] && $ret_vals[1])
				{ $num_sequences_with_primers++; }
				else
				{
					if ($ret_vals[0]) { $num_sequences_with_primers1++; }
					if ($ret_vals[1]) { $num_sequences_with_primers2++; }
				}
				$num_sequences++;
				
				#if ($ret_vals[3] eq 'p1' and length($ret_vals[5]) < 21) #< length($primers{'p1'}))
				#{
				#	print "r1/p1: $name\n";
				#}
				#if ($ret_vals[4] eq 'p2' and length($ret_vals[6]) < 15) #< length($primers{'p2'}))
				#{
				#	print "r2/p2: $name\n";
				#}
				
			
			}
			$name=$name_;
			$description=$description_;
			$sequence="";
		}
		else
		{
			$sequence .= "\U$line";
			if ($ftype eq 'FASTQ')
			{#skip next 2 lines
				$line=<IN>;
				$line=<IN>;
			}
			
		}
		#if ($num_sequences > 100000)
		#{
		#	last;
		#}
		
	}
	if ($name=~/\w/ and $sequence=~/\w/)
	{#do the same as above for the last sequence in the file
		my @ret_vals = remove_primers($sequence);
		
		#save the new sequence to a file, in the description indicate primers removed
		print OUT ">$name";
		if ($ret_vals[0]) { print OUT " $ret_vals[3]=$ret_vals[5]"; }
		else { print OUT " X"; }
		if ($ret_vals[1]) { print OUT " $ret_vals[4]=$ret_vals[6]"; }
		else { print OUT " X"; }
		print OUT "\n$ret_vals[2]\n";
		
		if ($ret_vals[0] && $ret_vals[1])
		{ $num_sequences_with_primers++; }
		else
		{
			if ($ret_vals[0]) { $num_sequences_with_primers1++; }
			if ($ret_vals[1]) { $num_sequences_with_primers2++; }
		}
		$num_sequences++;
			
	}
	close(IN);
}
print "Primers removed finished: $count_fasta files processed, $num_sequences sequences processed, $num_sequences_with_primers sequences found with both primers.\n";
print "$num_sequences_with_primers1 sequences found with only R1 primer, $num_sequences_with_primers2 found with only R2 primer.\n";
if ($count_fasta == 0)
{
	print "Warning: No fasta files processed ($!)!\n";
	exit(1);
}
close(OUT);
print "Exiting program...\n";
exit(0);

sub remove_primers
{
	my $seq = shift;
	my $seq_len = length($seq);
	my $fixed_seq = '';
	#look for p1 adaptor at position 0+12 of sequence (all 4 versions) - allow 8 mismatches
	#look for p2 adaptor at position 0+13 of sequence - allow 4 mismatches
	#see which one matches the best, and remove based on the best match
	my $r1_max_match_ratio = -1;
	my $r1_max_match_mismatches = 0;
	my $r1_max_match_p_id = '';
	my $r1_max_match_p_seq = '';
	my $r1_max_match_clip_pos = -1;
	my $r2_max_match_ratio = -1;
	my $r2_max_match_mismatches = 0;
	my $r2_max_match_p_id = '';
	my $r2_max_match_p_seq = '';
	my $r2_max_match_clip_pos = -1;
	my $r1_primer_found = 0; my $r2_primer_found = 0;
	foreach my $p_id (keys %primers)
	{
		my @start_pos_range = (-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20); 
		foreach my $start_pos (@start_pos_range)
		{ #positions in the sequence to look for adaptor
			if ($p_id eq 'p2' && $start_pos < -1*$max_p2_allowed_mismatches) { next; }
			
			my $start_at; my $end_at;
			if($p_id eq 'p2')
			{
			    $start_at = $start_pos;
			    $end_at = $start_pos+19;
			    
			}
			else
			{
				$start_at = $start_pos;
				$end_at = $start_pos+29;
			}
			
			#compare
                        my $primer_seq = $primers{$p_id};        
			if ($start_at < 0)
			{
				$primer_seq = substr($primer_seq,-1*$start_pos);
				$start_at = 0;
			}
			
			my $compare_seq = substr($seq, $start_at, $end_at-$start_at);
                        my $num_matches = 0;
			for(my $j = 0; $j < length($compare_seq); $j++) 
			{
				my $next_nucl = substr($compare_seq,$j,1); 
				my $primer_nucl = substr($primer_seq,$j,1);
				if ($p_id eq 'p1' and defined $p1_alt_nuc{$j})
				{#check both possibilities
					if($next_nucl eq $primer_nucl or $next_nucl eq $p1_alt_nuc{$j})
					{
						$num_matches++;
					}
				}
				else
				{
					if($next_nucl eq $primer_nucl)
					{
						$num_matches++;
					}
				}
				
			}
                        if ($num_matches/length($primers{$p_id}) > $r1_max_match_ratio)
			{
				$r1_max_match_ratio = $num_matches/length($primers{$p_id});
				$r1_max_match_p_id = $p_id;
				$r1_max_match_p_seq = $compare_seq;
				$r1_max_match_mismatches = length($primers{$p_id})-$num_matches;
				$r1_max_match_clip_pos = $end_at;
			}
		}
	}
	
	if( ($r1_max_match_p_id eq 'p2' && ($r1_max_match_mismatches <= $max_p2_allowed_mismatches)) ||
            ($r1_max_match_p_id ne 'p2' and ($r1_max_match_mismatches <= $max_p1_allowed_mismatches)) )
	{ $r1_primer_found = 1; }
	
	#look for r2 primer:
	#look for rev. complement of p1 adaptor at position seq_len-12 of sequence (all 4 versions) - allow 8 mismatches
	#look for rev. complement of p2 adaptor at position seq_len-13 of sequence - allow 4 mismatches
	#see which one matches the best, and remove based on the best match
	
	foreach my $p_id (keys %primers_rc)
	{
		my @start_pos_range = (-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20); #
		foreach my $start_pos (@start_pos_range)
		{ #positions in the sequence to look for adaptor
			if ($p_id eq 'p2' && $start_pos < -1*$max_p2_allowed_mismatches) { next; }
			
			my $start_at; my $end_at;
			if($p_id eq 'p2')
			{
				$start_at = $seq_len-$start_pos-19; 
				$end_at = $seq_len-$start_pos;
			    
			}
			else
			{
				$start_at = $seq_len-$start_pos-29; 
				$end_at = $seq_len-$start_pos;
			}
			
			#compare
			my $primer_seq = $primers_rc{$p_id};
			if ($end_at > $seq_len)
			{
				$primer_seq = substr($primer_seq,0,$start_pos);
				$end_at = $seq_len;
			}
			
			my $compare_seq = substr($seq,$start_at,$end_at-$start_at);
			my $num_matches = 0;
			for(my $j = 0; $j < length($compare_seq); $j++) 
			{
				my $next_nucl = substr($compare_seq,$j,1); 
				my $primer_nucl = substr($primer_seq,$j,1); 
				
				if ($p_id eq 'p1' and defined $p1_alt_nuc_rc{$j})
				{#check both possibilities
					if($next_nucl eq $primer_nucl or $next_nucl eq $p1_alt_nuc_rc{$j})
					{
						$num_matches++;
					}
				}
				else
				{
					if($next_nucl eq $primer_nucl)
					{
						$num_matches++;
					}
				}
			}
			if ($num_matches/length($primers_rc{$p_id}) > $r2_max_match_ratio)
			{
				$r2_max_match_ratio = $num_matches/length($primers_rc{$p_id});
				$r2_max_match_p_id = $p_id;
				$r2_max_match_p_seq = $compare_seq;
				$r2_max_match_mismatches = length($primers_rc{$p_id})-$num_matches;
				$r2_max_match_clip_pos = -1*($seq_len-$start_at);
			}   	
		}
	}
	if( ($r2_max_match_p_id eq 'p2' && ($r2_max_match_mismatches <= $max_p2_allowed_mismatches)) ||
	    ($r2_max_match_p_id ne 'p2' && ($r2_max_match_mismatches <= $max_p1_allowed_mismatches)) )
	{ $r2_primer_found = 1; }
	
	#output fixed sequence to new fasta file
	if($r2_primer_found)
	{
		$fixed_seq = substr($seq,0,$r2_max_match_clip_pos);
	}
	else { $fixed_seq = $seq; }
	
	if ($r1_primer_found)
	{
		$fixed_seq = substr($fixed_seq,$r1_max_match_clip_pos);
	}

	return ($r1_primer_found, $r2_primer_found, $fixed_seq, $r1_max_match_p_id, $r2_max_match_p_id, $r1_max_match_p_seq, $r2_max_match_p_seq);
}