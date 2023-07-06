#!/usr/local/bin/perl 
#
#require "./masses_and_fragments.pl";

# reads each file in the directory - either the first argument on the command line, or the current directory
#for each fasta file in the directory, the name/sequence is read in and then the protein is digested with 
#trypsin to get the resulting peptides - 
#a single output file is created - "all_predigested.fasta" that contains, 
#for the description line: the peptide sequence, followed by a number representing the number of proteins that 
#  resulted in that peptide when digested with trypsin
#for the sequence line: the peptide sequence

#use warnings;
use strict;
use File::Basename;

my $filename;
my $out_dir;
my $incompletes;
my $use_tail;
my $use_chymotrypsin;

if ($ARGV[0] && $ARGV[0]=~/\w/) { $filename=$ARGV[0];}  
else { $filename = ""; }

if ($ARGV[1] && $ARGV[1]=~/\w/) { $out_dir=$ARGV[1];}  
else { $out_dir = ""; }

if ($ARGV[2] && $ARGV[2]=~/\w/) { $incompletes=$ARGV[2];}
else { $incompletes=1; }    # Trypsin:1, Chymotrypsin:4

if ($ARGV[3] && $ARGV[3]=~/\w/) { $use_tail=$ARGV[3];}
else { $use_tail=0; }

if ($ARGV[4] && $ARGV[4]=~/\w/) { $use_chymotrypsin=$ARGV[4];}
else { $use_chymotrypsin=0; }

my %protein_peptides = ();

my $total_pep_count=0;
my %PEP=();
my %PEP_proteins=();
my %PEP_proteins_count=();
my $files_count=0;
my $peptides_count=0;
my $proteins_count=0;
my $line="";
my @TAIL_SEQUENCES = ("SEPKIPQPQPKPQ", "SAHHSEDPSSKCP", "SEPKTPKPQPQPQPQ", "SGTNEVCKCPKCPAPEL", "EPKIPQPQPKPQ", "AHHSEDPSSKCP", "EPKTPKPQPQPQPQ", "GTNEVCKCPKCPAPEL");

my $file_suffix="";
if ($use_chymotrypsin)
{
    $file_suffix="_chymotrypsin";
}

if (!open (IN,"$filename"))
{ print "Error opening $filename.\n"; exit(0); }

my $out_pre = basename($filename);
$out_pre =~ s/\.[^\.]+$//;

if (!open (IND_OUT, ">$out_dir/$out_pre" . "_protein_peptides" . $file_suffix . ".fasta"))
{ print "Error creating $out_dir/$out_pre" . "_protein_peptides" . $file_suffix . ".fasta\n"; close(IN); exit(0); }

if (!open (OUT,">$out_dir/$out_pre" . "_predigested" . $file_suffix . ".fasta"))
{ print "Error creating $out_dir/$out_pre" . "_predigested" . $file_suffix . ".fasta\n"; close(IN); close(IND_OUT); exit(0); }

my $name="";
my $sequence="";
my $output_check=0;
while ($line=<IN>)
{
	chomp($line);
	if ($line =~ s/^>//)
	{#if the current line starts with a '>', then it is the description line, remove the '>' and input the 'name'
		
		my $name_=$line;
		if ($name_ =~ s/^\"//)
		{#if the line begins with quotes (")
		 #remove quotes at beginning and end
		 #replace any non-word character (not letters or numbers) with '_'
			$name_ =~ s/\"$//;
			$name_ =~ s/[^\w]/_/g; 
		}
		else
		{#line does not begin with quotes ("), remove anything after (and including) the first whitespace character 
			$name_ =~ s/\s.*$//;
		}
		if ($name =~ /\w/ and $sequence =~ /\w/)
		{#both name and sequence were inputted, add the protein to the count
			if ($output_check and $proteins_count < 5) { print "$name\n$sequence\n"; }
			
			#digest the protein with trypsin, the function fills the %PEP hash with the resulting peptides
			%PEP=();
			
			#add tail sequence to protein so that we don't miss c terminus peptides:
			if($use_tail)
			{
				my $new_sequence;
				foreach my $tail (@TAIL_SEQUENCES)
				{
					$new_sequence = $sequence . $tail;
					
					if ($use_chymotrypsin) {
						my $str = DigestChymoTrypsin($name,$new_sequence,$incompletes);
						if ($output_check and $proteins_count < 5) { print "$str\n"; }
					}
					else
					{
						my $str = DigestTrypsin($name,$new_sequence,$incompletes);
						if ($output_check and $proteins_count < 5) { print "$str\n"; }
					}
					
				}
			}
			else
			{
				if ($use_chymotrypsin)
				{
					my $str = DigestChymoTrypsin($name,$sequence,$incompletes);
					if ($output_check and $proteins_count < 5) { print "$str\n"; }
				}
				else
				{
					my $str = DigestTrypsin($name,$sequence,$incompletes);
					if ($output_check and $proteins_count < 5) { print "$str\n"; }
				}
			}
			
			if($proteins_count % 50000 == 0) { print "$proteins_count\n"; } #if($protein_count >= 10000) { last; } }
			#if ($proteins_count > 10000) { last; }
			
			#strip primer info from name before outputting!
			$name =~ s/(_fr[012])_fwd_(\w*)_rev_(\w*)$//;
			if($1)
			{ $name = $name . $1; }
			
			print IND_OUT ">$name\n";
			foreach my $peptide (keys %PEP)
			{#we want to keep a count of the number of proteins that contained a particular peptide when digested 
			 #with trypsin - a list of each 'name' of the proteins is stored in the hash %PEP_proteins and the number 
			 #of proteins for a given peptide is stored in the hash %PEP_proteins_count
				if (length($peptide)>6)
				{#disregard peptides that are too short
					#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
					#{
					#	$PEP_proteins_count{$peptide}++;
					#	$PEP_proteins{$peptide} .= qq!#$name#!;
					#}
					#else
					#{
					#	print "$name $peptide\n";
					#}
					$PEP_proteins_count{$peptide}++;
					print IND_OUT "#$peptide#";
				}
			}
			print IND_OUT "#\n";
			
			$proteins_count++;
			
		}
		$name=$name_;
		$sequence="";
	}
	else
	{#it is a line of the sequence
		$sequence.="$line";
	}
}	

#the last protein left in the file, do the same as above
if ($name =~ /\w/ and $sequence =~ /\w/)
{
	#digest the protein with trypsin, the function fills the %PEP hash with the resulting peptides
	%PEP=();
	
	#add tail sequence to protein so that we don't miss c terminus peptides:
	if($use_tail)
	{
		my $new_sequence;
		foreach my $tail (@TAIL_SEQUENCES)
		{
			$new_sequence = $sequence . $tail;
			
			if ($use_chymotrypsin)
			{
				DigestChymoTrypsin($name,$new_sequence,$incompletes);
			}
			else
			{
				DigestTrypsin($name,$new_sequence,$incompletes);
			}
		}
	}
	else
	{
		if ($use_chymotrypsin)
		{
			DigestChymoTrypsin($name,$sequence,$incompletes);
		}
		else
		{
			DigestTrypsin($name,$sequence,$incompletes);
		}
	}
	
	#strip primer info from name before outputting!
	$name =~ s/(_fr[012])_fwd_(\w*)_rev_(\w*)$//;
	if($1) { $name = $name . $1; }
	
	print IND_OUT ">$name\n";
	foreach my $peptide (keys %PEP)
	{
		#print qq!$peptide\n!;
		if (length($peptide)>6)
		{
			#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
			#{
			#	$PEP_proteins_count{$peptide}++;
			#	$PEP_proteins{$peptide} .= qq!#$name#!;
			#}
			$PEP_proteins_count{$peptide}++;
			print IND_OUT "#$peptide#";
		}
	}
	print IND_OUT "#\n";
	
	$proteins_count++;
}
close(IN);
close(IND_OUT);

print qq!$filename $proteins_count\n!;

foreach my $peptide (keys %PEP_proteins_count) #(keys %PEP_proteins)
{
    print OUT qq!>$peptide $PEP_proteins_count{$peptide}\n$peptide\n!;
}
close(OUT);
print("DONE\n");

sub DigestChymoTrypsin
{ # arg1 - name, arg2 = sequence, arg3 = # of incompletes
	# cleavage at c-terminal of FLYWM with 0-4 missed cleavages
	
	my $name = shift();
	my $seq = shift();
	my $incompletes = shift();
	
	my $pep_str="";
	
	my $temp=$seq;
	my @pep=();
	my @start=();
	my @end=();
	my $aa="";
	my $aa_="";
	my $i=0;

	for($i=0;$i<=$incompletes;$i++)
	{
		$start[$i]=0;
		$end[$i]=-1;
		#$pep[$i]="[";
	}
	
	my $aa_count=0;
	while ($temp =~ s/^\s*([A-Z\*])//)
	{
		$aa="\U$1";
		
		
		if ( (($aa_=~/F/ or $aa_=~/L/ or $aa_=~/Y/ or $aa_=~/W/ or $aa_=~/M/) and $aa!~/P/) or $aa_=~/\*/)
		{
			for($i=0;$i<=$incompletes;$i++)
			{
				if(defined $pep[$i] and (not ($PEP{"$pep[$i]"}==1)))
				{
					$PEP{"$pep[$i]"}=1;
					$pep_str.="\n$pep[$i]";
				}
				$pep[$i]=$pep[$i+1];
				$start[$i]=$start[$i+1];
				$end[$i]=$end[$i+1];
			}
			$start[$incompletes]=$aa_count;
			$end[$incompletes]=$aa_count-1;
		}
		for($i=0;$i<=$incompletes;$i++)
		{
			my $aa__=$aa;
			$aa__=~s/I/L/g;
			if ($aa__!~/\*/) { $pep[$i].=$aa__; }
			$end[$i]++;
		}
		$aa_=$aa;
		$aa_count++;
	}
	for($i=0;$i<=$incompletes;$i++)
	{
		if(defined $pep[$i] and (not ($PEP{"$pep[$i]"}==1)))
		{
			$PEP{"$pep[$i]"}=1;
			$pep_str.="\n$pep[$i]";
		}
		
	}
	
	return $pep_str;
}

sub DigestTrypsin
{#fills the %PEP hash with the peptides resulting from digesting the sequence with trypsin
 #arg1 - name, arg2 = sequence, arg3 = # of incompletes
	my $name = shift();
	my $seq = shift();
	my $incompletes = shift();
	
	my $pep_str="";
	
	my $temp=$seq;
	my @pep=();
	my @start=();
	my @end=();
	my $aa="";
	my $aa_="";
	my $i=0;

	for($i=0;$i<=$incompletes;$i++)
	{
		$start[$i]=0;
		$end[$i]=-1;
		#$pep[$i]="[";
	}
	my $aa_count=0;
	while ($temp =~ s/^\s*([A-Z\*])//)
	{
		$aa="\U$1";
		$aa=~s/I/L/g;
		if ( (($aa_=~/R/ or $aa_=~/K/) and $aa!~/P/) or $aa_=~/\*/)
		{
			for($i=0;$i<=$incompletes;$i++)
			{
				if(defined $pep[$i] and (not ($PEP{"$pep[$i]"}==1)))
				{
					$PEP{"$pep[$i]"}=1;
					$pep_str.="\n$pep[$i]";
				}
				$pep[$i]=$pep[$i+1];
				$start[$i]=$start[$i+1];
				$end[$i]=$end[$i+1];
			}
			$start[$incompletes]=$aa_count;
			$end[$incompletes]=$aa_count-1;
		}
		for($i=0;$i<=$incompletes;$i++)
		{
			if ($aa!~/\*/) { $pep[$i].=$aa; }
			$end[$i]++;
		}
		$aa_=$aa;
		$aa_count++;
	}
	for($i=0;$i<=$incompletes;$i++)
	{
		if(defined $pep[$i] and (not ($PEP{"$pep[$i]"}==1)))
		{
			$PEP{"$pep[$i]"}=1;
			$pep_str.="\n$pep[$i]";
		}
	}
	return $pep_str;
}


