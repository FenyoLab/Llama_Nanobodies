#!/usr/local/bin/perl 

use warnings;
use strict;

my $filename="";
my $filename2="";
my $out_dir = "";
my $use_primers=1;
my $new_primers=1;
my $remove_singletons=1;  # also outputs file with singletons removed (protein sequences that come from only 1 DNA sequence)
my $output_full_db=1;

if ($ARGV[0]=~/\w/) { $filename="$ARGV[0]"; }
else { $filename = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/code/NCDIR/Llama/dna_removed_primersL1.fasta"; }

if ($ARGV[1]=~/\w/) { $filename2="$ARGV[1]"; }
else { $filename2 = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/code/NCDIR/Llama/dna_removed_primersL2.fasta"; } 

if ($ARGV[2]=~/\w/) { $out_dir="$ARGV[2]"; }
else { $out_dir = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/code/NCDIR/Llama/"; } 

if ($ARGV[3]=~/\w/) { $use_primers=$ARGV[3]; }
else { $use_primers=1; }

if ($ARGV[4]=~/\w/) { $new_primers=$ARGV[4]; }
else { $new_primers=1; }

if ($ARGV[5]=~/\w/) { $remove_singletons=$ARGV[5]; }
else { $remove_singletons=1; }

if ($ARGV[6]=~/\w/) { $output_full_db=$ARGV[6]; }
else { $output_full_db=1; }

print "filename=$filename, filename2=$filename2, out_dir=$out_dir, use_primers=$use_primers, new_primers=$new_primers, remove_singletons=$remove_singletons, output_full_db=$output_full_db\n";
if ($filename2 eq "none")
{
    $filename2="";
}

my %redundancy_listing;
my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
		"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
		"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
		"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
		
		"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
		"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
		"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
		"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
		
		"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
		"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
		"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
		"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
		
		"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
		"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
		"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
		"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");


my $LENGTH_P1_PRIMER;
if ($new_primers)
{
    print "Using new primers\n";
    $LENGTH_P1_PRIMER = 23;
}
else
{
	$LENGTH_P1_PRIMER = 29;
	print "Using old primers\n";
}

print "Executing dnatoprot_longest_nr.pl\n";

my %longest_orf=();
my %filtered_longest_orf=();
my %orf=();
my %orf_trim=();

my $BOTH_PRIMERS = 1;
my $OUTPUT_FILTERED_SEQUENCES = 0; #set this to false when running through llama magic

if ($output_full_db)
{
    if (!open (OUT,">$out_dir/longest_nr.fasta")) { print "Error creating $out_dir/longest_nr.fasta ($!)\n"; exit(0); }
    print "Opened $out_dir/longest_nr.fasta\n";
}

if($remove_singletons)
{
    if (!open (OUT2,">$out_dir/longest_nr_no_singletons.fasta"))
    {
        print "Error creating $out_dir/longest_nr_no_singletons.fasta ($!)\n";
        if ($output_full_db) { close(OUT); }
        exit(0);
    }
    print "Opened $out_dir/longest_nr_no_singletons.fasta\n";
}

### Just testing opening of input files
### want to make sure they can both be opened before starting
if (!open (IN,"$filename"))
{
    print "Error opening $filename ($!).\n";
    if ($output_full_db) { close(OUT); }
    if ($remove_singletons) { close(OUT2); }  
}
else { close(IN); }
    
if ($filename2)
{
    if (!open (IN,"$filename2"))
    {
        print "Error opening $filename2 ($!).\n";
        if ($output_full_db) { close(OUT); }
        if ($remove_singletons) { close(OUT2); }
        exit(0);
    }
    else { close(IN); }
}

###################################################
print("Reading input...\n");
my $err = read_input($filename);
if (!$err && $filename2) { $err = read_input($filename2); }

if ($err)
{
    if ($output_full_db) { close(OUT); }
    if ($remove_singletons) { close(OUT2); }
    exit(0);
}

#OUT will contain the set of the unique, longest length reading frames found with a count in their description of how many dna sequences 
#resulted in this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
print("Writing output...\n");
my $total_count=0;
my $unique_count=0;
my $non_singletons_count=0;
my $size = keys %longest_orf;
print("Number of keys: $size.\n");
foreach my $seq (sort keys %longest_orf)
{
    my $temp=$longest_orf{$seq};
    if ($temp =~ s/^([^\,]+)\, //) #get characters up to next ', '
    {
        my $name=$1;
        my $count=0;
        
        #get fwd and rev primer - put them back onto sequence
        my $p1 = '';
        my $p2 = '';
        if ($use_primers)
        {
            if($name =~ /(_fr[012])_fwd_(\w*)_rev_(\w*)$/)
            {
                my $fr = $1;
                $p1 = $2;
                $p2 = $3;
            }
        }
        
        my $found_name = '';
        while($temp =~ s/^([^\,]+)\, //) #count how many dna sequences result in current $seq and print result
        {
            $count++;
            my $next_name = $1;
        }
        
        if ($count==0) 
        { 
            print OUT qq!>$name\n$p1$seq$p2\n!; 
            $total_count++;
        }
        else
        {
            if($output_full_db) { print OUT qq!>$name + $count other\n$p1$seq$p2\n!; }
            if($remove_singletons) { print OUT2 qq!>$name + $count other\n$p1$seq$p2\n!; $non_singletons_count++; }
            #print OUT qq!>$longest_orf{$seq}\n$p1$seq$p2\n!;
            $total_count += ($count+1);
        }
        $unique_count++;
        if ($unique_count % 100000 == 0)
        {
            print "$unique_count\n";
        }
    }
    else
    {
        print "Error in outputting to file.\n";
    }
}

if($output_full_db) { close(OUT); }
if($remove_singletons) { close(OUT2); }
print "Total non-singletons output: $non_singletons_count\n";
print "Total sequences output: $unique_count\n";
print "Total sequences output (including duplicates): $total_count\n";


if ($OUTPUT_FILTERED_SEQUENCES)
{
	if (open (OUT,">$out_dir/longest_nr.filtered.fasta"))
	{#this file will contain the set of the unique, longest length reading frames found with a count in their description of how many dna sequences 
	 #resulted in this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
		print "Success opening $out_dir/longest_nr.filtered.fasta for writing.\n";
		my $total_count=0;
		my $unique_count=0;
		foreach my $seq (sort keys %filtered_longest_orf)
		{
			my $temp=$filtered_longest_orf{$seq};
			if ($temp =~ s/^([^\,]+)\, //) #get characters up to next ', '
			{
				my $name=$1;
				my $count=0;
				
				#get fwd and rev primer - put them back onto sequence
				my $p1 = '';
				my $p2 = '';
				if ($use_primers)
				{
					if($name =~ /(_fr[012])_fwd_(\w*)_rev_(\w*)$/)
					{
						my $fr = $1;
						$p1 = $2;
						$p2 = $3;
					}
				}
				
				my $found_name = '';
				while($temp =~ s/^([^\,]+)\, //) #count how many dna sequences result in current $seq and print result
				{
					$count++;
					my $next_name = $1;
					
				}
				
				if ($count==0) 
				{ 
					print OUT qq!>$name\n$p1$seq$p2\n!; 
					$total_count++;
				}
				else
				{
					print OUT qq!>$name + $count other\n$p1$seq$p2\n!;
					#print OUT qq!>$filtered_longest_orf{$seq}\n$p1$seq$p2\n!;
					$total_count += ($count+1);
				}
				$unique_count++;
			}
			else
			{
				print "Error in outputting to file."
			}
		}
		close(OUT);
		print "Total sequences output: $unique_count\n";
		print "Total sequences output (including duplicates): $total_count\n";
	}
	else
	{
		print "Error creating $out_dir/longest_nr.filtered.fasta ($!)\n";
		exit(0);
	}
    my $filtered_sequences = scalar keys %filtered_longest_orf;
    print "Sequences filtered: $filtered_sequences\n";
}

sub read_input
{
    my $fname = shift;
    if (!open (IN,"$fname")) { print "Error opening $fname ($!).\n"; return 1; }
    print "Opened $fname\n";

    my $entry_count = 0;
    my $name="";
    my $description="";
    my $sequence="";
    my $line="";
    
    while ($line=<IN>)
    {
        chomp($line);
        if ($line=~/^[>@](\S+)\s?(.*)$/)
        {#if the current line is the description line - begins with '>', and then 
         #one or more non-whitespace, followed by 0 or more whitespace char's, then anything until the end
            
            my $name_=$1; 
            my $description_=$2;
            $name_ =~ s/[\:\-\\\/]/_/g;
            if ($name=~/\w/ and $sequence=~/\w/)
            {#the entire sequence has been read in, so do the conversion:
            
                $entry_count++;
                if ($entry_count % 100000 == 0)
                {
                    print "$entry_count\n";
                }
            
                my @ret_val = translate_dna($sequence, $description);
                my $sequence_filtered = 1;
                if ($ret_val[0] > 6)
                {#seq must be > 6 in length
                    if ($use_primers)
                    {
                        if($ret_val[1] !~ /\*/ && $ret_val[1] !~ /X$/)
                        {#no stop codons in primer-translated seq and no X at the end (frame shift from primer reading frame due to insertion/deletion)
                            if (!$BOTH_PRIMERS || ($ret_val[4] && $ret_val[5]))
                            {#if flag set, require both primers found to output sequence
                                $longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3]_fwd_$ret_val[4]_rev_$ret_val[5], ";
                                $sequence_filtered = 0;
                            } 
                        } 
                    }
                    else
                    {
                        $longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3], ";
                        $sequence_filtered = 0;
                    }
                }
                if ($sequence_filtered)
                {
                    if ($use_primers)
                    {
                        $filtered_longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3]_fwd_$ret_val[4]_rev_$ret_val[5], ";
                    }
                    else
                    {
                        $filtered_longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3], ";
                    }
                }
            }
            $name=$name_;
            $description=$description_;
            $sequence="";
        }
        else
        {
            $sequence .= "\U$line";
        }
    }	
    if ($name=~/\w/ and $sequence=~/\w/)
    {#do the same as above for the last sequence in the file
        $entry_count++;
        my @ret_val = translate_dna($sequence, $description);
        my $sequence_filtered = 1;
        if ($ret_val[0] > 6)
        {#seq must be > 6 in length
            if ($use_primers)
            {
                if($ret_val[1] !~ /\*/ && $ret_val[1] !~ /X$/)
                {#no stop codons in primer-translated seq and no X at the end (frame shift from primer reading frame due to insertion/deletion)
                    if (!$BOTH_PRIMERS || ($ret_val[4] && $ret_val[5]))
                    {#if flag set, require both primers found to output sequence
                        $longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3]_fwd_$ret_val[4]_rev_$ret_val[5], ";
                        $sequence_filtered = 0;
                    } 
                } 
            }
            else
            {
                $longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3], ";
                $sequence_filtered = 0;
            }
        }
        if ($sequence_filtered)
        {
            if ($use_primers)
            {
                $filtered_longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3]_fwd_$ret_val[4]_rev_$ret_val[5], ";
            }
            else
            {
                $filtered_longest_orf{$ret_val[1]} .= "$name\_$ret_val[2]\_fr$ret_val[3], ";
            }
        }
        
    }
    close(IN);
    return 0;
}

sub translate_dna
{
	my $sequence = shift;
	my $description = shift;
	my $size = length($sequence);
	
	#if primers info is present, read it in, if atleast one primer found, use the info
	#else, translate with longest orf
	if ($use_primers && $description !~ /^X X$/)
	{
		my $frame;
		my $direction;
		
		my $r1_primer="";
		my $r1_primer_seq="";
		my $r2_primer="";
		my $r2_primer_seq="";
		
		if($new_primers)
		{
			$description =~ /^([pX][12]?_?[SL]?H?)=?(\w*) ([pX][12]?_?[SL]?H?)=?(\w*)/;
			$r1_primer = $1;
			$r1_primer_seq = $2;
			$r2_primer = $3;
			$r2_primer_seq = $4;
		}
		else
		{
			$description =~ /^([pX]\S?)=?(\w*) ([pX]\S?)=?(\w*)/;
			$r1_primer = $1;
			$r1_primer_seq = $2;
			$r2_primer = $3;
			$r2_primer_seq = $4;
		}
		
		if($r1_primer eq 'X' or $r2_primer eq 'X')
		{
			my $hi = 1;
		}
		if ($r1_primer eq 'p1')
		{
			$direction = 'fwd';
			$frame = 1;
			
			#correct for unusual case where same primer found in both r1 and r2, by ignoring r2 primer
			if ($r2_primer eq 'p1')
			{
				$r2_primer_seq = '';
			}
			
		}
		elsif($r1_primer eq 'p2' || $r1_primer eq 'p2_LH' || $r1_primer eq 'p2_SH')
		{
			if ($r2_primer eq 'p1')
			{
				$direction = 'rev';
				$frame = 1;
			}
			else
			{ #r1 primer is p2 but p1 not found in r2
				$direction = 'rev';
				if ($r1_primer eq 'p2_SH')
				{
                    $frame = ($size % 3) + 1;
					if ($frame == 3) { $frame = 0; }
                }
                else
				{
					$frame = $size % 3;
				}
				
				#correct for unusual case where same primer found in both r1 and r2, by ignoring r2 primer
				if ($r2_primer eq 'p2' || $r2_primer eq 'p2_LH' || $r2_primer eq 'p2_SH')
				{
					$r2_primer_seq = '';
				}
			}
		}
		else
		{#r1 primer not found
			if ($r2_primer eq 'p1')
			{
				$direction = 'rev';
				$frame = 1;
			}
			else # $r2 primer is p2 but p1 not found in r1
			{
				$direction = 'fwd';
				if ($r2_primer eq 'p2_SH')
				{
                    $frame = ($size % 3) + 1;
					if ($frame == 3) { $frame = 0; }
                }
                else
				{
					$frame = $size % 3;
				}
			}
		}
		my $seq;
		if ($direction =~ /^fwd$/)
		{
			$seq=$sequence;
			if($r2_primer eq 'p2_SH')
			{ $seq = $seq . 'G'; }
		} 
		else 
		{
			$seq = reverse $sequence;
			$seq =~ tr/ATCG/TAGC/;
			if($r1_primer eq 'p2_SH')
			{ $seq = $seq . 'C'; }
		}
		
		my @ret_vals = translate_frame($seq,$frame,$size, 0);	
			
		if ($direction eq 'fwd')
		{
			#p2 translation
			if ($r2_primer_seq)
			{
				my @ret_vals;
				if ($r2_primer eq 'p2_SH')
				{
					@ret_vals = translate_frame($r2_primer_seq, 1, length($r2_primer_seq), 1);
                }
				else
				{# p2_LH or p2
					@ret_vals = translate_frame($r2_primer_seq, 0, length($r2_primer_seq), 1);
				}
                
				$r2_primer_seq = $ret_vals[1];
			}
			#p1 translation
			if ($r1_primer_seq)
			{
				my $r1_frame=0;
				if (length($r1_primer_seq) < $LENGTH_P1_PRIMER)
				{
					my $l_diff = $LENGTH_P1_PRIMER - length($r1_primer_seq);
					if (($l_diff % 3) == 2) { $r1_frame = 1; }
					if (($l_diff % 3) == 1) { $r1_frame = 2; }
					if (($l_diff % 3) == 0) { $r1_frame = 0; }
				}
				
				if ($new_primers) { $r1_primer_seq = $r1_primer_seq . 'T'; }
				else { $r1_primer_seq = $r1_primer_seq . 'G'; }
                
				my @ret_vals = translate_frame($r1_primer_seq, $r1_frame, length($r1_primer_seq), 1);
				$r1_primer_seq = $ret_vals[1];
			}
			
			return ($ret_vals[0], $ret_vals[1], $direction, $frame, $r1_primer_seq, $r2_primer_seq);
		}
		else
		{
			#p2 translation
			if ($r1_primer_seq)
			{
				$r1_primer_seq = reverse $r1_primer_seq;
				$r1_primer_seq =~ tr/ATCG/TAGC/;
				
				my @ret_vals;
				if ($r1_primer eq 'p2_SH')
				{
					@ret_vals = translate_frame($r1_primer_seq, 1, length($r1_primer_seq), 1);
				}
				else
				{
					@ret_vals = translate_frame($r1_primer_seq, 0, length($r1_primer_seq), 1);
				}
				$r1_primer_seq = $ret_vals[1]; 
			}
			
			#p1 translation
			if ($r2_primer_seq)
			{
				my $r1_frame=0;
				if (length($r2_primer_seq) < $LENGTH_P1_PRIMER)
				{
					my $l_diff = $LENGTH_P1_PRIMER - length($r2_primer_seq);
					if (($l_diff % 3) == 2) { $r1_frame = 1; }
					if (($l_diff % 3) == 1) { $r1_frame = 2; }
					if (($l_diff % 3) == 0) { $r1_frame = 0; }
				}
				
				if ($new_primers) { $r2_primer_seq = 'A' . $r2_primer_seq; }
				else { $r2_primer_seq = 'C' . $r2_primer_seq; }
				
				$r2_primer_seq = reverse $r2_primer_seq;
				$r2_primer_seq =~ tr/ATCG/TAGC/;
				
				my @ret_vals = translate_frame($r2_primer_seq, $r1_frame, length($r2_primer_seq), 1);
				$r2_primer_seq = $ret_vals[1];
			}
			
			return ($ret_vals[0], $ret_vals[1], $direction, $frame, $r2_primer_seq, $r1_primer_seq);
		}
	}
	else
	{
		my $longest_orf_direction="";
		my $longest_orf_frame="";
		my $longest_orf_length=0;
		my $longest_orf_seq="";
		foreach my $direction ("fwd","rev")
		{#read both forward and reverse
			my $seq="";
			if ($direction =~ /^fwd$/) { $seq=$sequence; } 
			else 
			{ 
				$seq = reverse $sequence;
				$seq =~ tr/ATCG/TAGC/;
			}
			for (my $k=0;$k<3;$k++)
			{#for each reading frame:
				my @ret_vals = translate_frame($seq, $k, $size, 1);
				if ($longest_orf_length<$ret_vals[0])
				{#store the info for the longest orf, will print only longest orf's to a separate file
					$longest_orf_direction=$direction;
					$longest_orf_frame=$k;
					$longest_orf_length=$ret_vals[0];
					$longest_orf_seq=$ret_vals[1];
				}
				
			}
		}
		return ($longest_orf_length, $longest_orf_seq, $longest_orf_direction, $longest_orf_frame, "", "");
	}
}

sub translate_frame
{
	my $seq = shift;
	my $k = shift; #the reading frame
	my $size = shift;
	my $split_orfs = shift;
	
	my $longest_orf_length=0;
	my $longest_orf_seq="";
	
	my $protein="";
	for(my $n=$k;$n<$size;$n=$n+3)
	{
		my $triplet = substr($seq, $n, 3);
		if ($mapping{$triplet} and $mapping{$triplet} =~ /[\w\*]/) { $protein .= $mapping{$triplet}; } # '*' is stop codon
		else { $protein.="X"; } # X is unknown, doesn't code for anything must be error in sequence
	}
	
	#remove X's at beginning and end
	$protein =~ s/X+$//;
	$protein =~ s/^X+//;
	
	if (!$split_orfs)
	{
		return (length($protein),$protein);
	}
	
	my $temp="$protein*";
	my $index=0;
	
	#find all the orf's - the sequence up to the next stop codon
	while ($temp =~ s/^([^\*]*)\*//) 
	{# $1 is 0 or more of any char's except '*', followed by '*', starting at the beginning
	 #so we get the AA sequence up to the next stop codon, and remove it
		my $orf_seq=$1;
		
		#remove X at beginning and end
		$orf_seq =~ s/X+$//;
		$orf_seq =~ s/^X+//;
		
		my $orf_length=length($orf_seq); 
		$index+=$orf_length+1;
		if ($longest_orf_length<$orf_length)
		{#store the info for the longest orf, will print only longest orf's to a separate file
			$longest_orf_length=$orf_length;
			$longest_orf_seq=$orf_seq;
		}
	}
	#count and return how many orf's ( > 1 means stop codon in the seq) - record this for seq'a with both primers found,e tc.
	return ($longest_orf_length,$longest_orf_seq);
}


