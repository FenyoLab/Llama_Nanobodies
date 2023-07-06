#!/usr/local/bin/perl 

use warnings;
use strict;

my $dir="";
my $out_dir = "";
my $seqs_per_file = 100000;
my $split_fasta = 0;

if ($ARGV[0] && $ARGV[0]=~/\w/) { $dir="$ARGV[0]";} else
{ $dir=""; }
#"C:\\NCDIR\\Llama\\FASTQ\\Merged\\08_05_2014\\R1Trim50R2Trim50MinMatch50Identity0.95"; }
#{ $dir="C:\\NCDIR\\Llama\\FASTQ\\Merged\\10_07_2013\\4762_merged_iden0.9_notrim"; } #"C:\\NCDIR\\Llama\\FASTQ\\Reads\\08_05_2014"; }

#{ $dir="/ifs/data/proteomics/projects/NCDIR/Llama/data-and-results/db/MISEQ_DB_2014_08_05/ABySSMerge"; }
#{ $dir="C:\\NCDIR\\Llama\\FASTQ\\Merged\\08_05_2014"; }

if ($ARGV[1] && $ARGV[1]=~/\w/) { $out_dir="$ARGV[1]";} else
{ $out_dir=""; }
#"C:\\NCDIR\\Llama\\FASTQ\\Merged\\08_05_2014\\R1Trim50R2Trim50MinMatch50Identity0.95\\fasta"; }
#{ $out_dir="C:\\NCDIR\\Llama\\FASTQ\\Merged\\10_07_2013\\4762_merged_iden0.9_notrim"; } #"C:\\NCDIR\\Llama\\FASTQ\\Reads\\08_05_2014\\fasta"; }

#{ $out_dir="/ifs/home/snk218/llama/fasta-split"; }
#{ $out_dir="C:\\NCDIR\\Llama\\FASTQ\\Merged\\08_05_2014\\fasta-split"; } 

if (!opendir(DIR,"$dir")) { print "Error reading $dir ($!)\n"; exit(2); }

my @allfiles=readdir DIR;
closedir DIR;
my $count_fasta = 0;
foreach my $filename (@allfiles)
{#for each fastq file  
	my $ftype;
	if($filename =~ /\.fastq$/i || $filename =~ /\.fq$/i) { $ftype = 'FASTQ'; }
	else { next; }
	
	if (!open (IN,"$dir/$filename")) { print "Error opening $dir/$filename ($!).\n"; next; }
	print "Opened $dir/$filename\n";
	
	my $file_count = 1;
	my $out_count = 1;
	if (!open (OUT,">$out_dir/$filename-$file_count.fasta")) { print "Error opening $out_dir/$filename-$file_count.fasta ($!).\n"; next; }
	print "Opened $out_dir/$filename-$file_count.fasta\n";
	
	$count_fasta++;
	
	my $name="";
	my $description="";
	my $sequence="";
	my $line="";
	my $entry_count = 0;
	
	while ($line=<IN>)
	{
		chomp($line);
		if ($line=~/^@(\S+)\s?(.*)$/)
		{#if the current line is the description line - begins with '@', and then 
		 #one or more non-whitespace, followed by 0 or more whitespace char's, then anything until the end
			$entry_count++;
			#if ($entry_count % 1000 == 0) { print "$entry_count\n"; }
			my $name_=$1; 
			my $description_=$2;
			$name_ =~ s/[\:\-\\\/]/_/g;
			if ($name=~/\w/ and $sequence=~/\w/)
			{#the entire sequence has been read in, so do the conversion:
				
				print OUT qq!>$name\n$sequence\n!;
				if ($split_fasta && $out_count % $seqs_per_file == 0)
				{
					close(OUT);
					$file_count++;
					if (!open (OUT,">$out_dir/$filename-$file_count.fasta")) { print "Error opening $out_dir/$filename-$file_count.fasta ($!).\n"; last; }
					print "Opened $out_dir/$filename-$file_count.fasta\n";
					
				}
				$out_count++;
				
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
	}	
	if ($name=~/\w/ and $sequence=~/\w/)
	{#do the same as above for the last sequence in the file
				
		print OUT qq!>$name\n$sequence\n!;
	}
	close(IN);
	close(OUT);
}
print "Conversions finished: $count_fasta files processed.\n";



