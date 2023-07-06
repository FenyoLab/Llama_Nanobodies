#!/usr/local/bin/perl 

use warnings;
use strict;

my %primers = ('p1'=>'GTCCTGGCTGCTCTTCTACAAGG',
               'p2_LH'=>'GTGGTTGTGGTTTTGGTGTCTTGGG',
               'p2_SH'=>'CTGGGGTCTTCGCTGTGGTGC');

my %primers_rc = ('p1'=>'CCTTGTAGAAGAGCAGCCAGGAC',  
                  'p2_LH'=>'CCCAAGACACCAAAACCACAACCAC',
                  'p2_SH'=>'GCACCACAGCGAAGACCCCAG');
               

my %p1_alt_nuc = (15 => 'T');
my %p1_alt_nuc_rc = (7 => 'A');

my $max_p1_allowed_mismatches = 8;
my $max_p2_LH_allowed_mismatches = 8;
my $max_p2_SH_allowed_mismatches = 8;

my $num_sequences_with_primers = 0;
my $num_sequences_with_primers1 = 0;
my $num_sequences_with_primers2 = 0;
my $num_sequences = 0;

my $filename="";
my $outfile = "";

if ($ARGV[0] && $ARGV[0]=~/\w/) { $filename="$ARGV[0]";}
else { $filename = ""; }

if ($ARGV[1] && $ARGV[1]=~/\w/) { $outfile="$ARGV[1]";}
else { $outfile = ""; }

print "Executing remove_primers.pl\n";

if(!open(OUT, ">$outfile")) { print "Error: Could not open $outfile for writing ($!).\n"; exit(0); }
print "Opened $outfile\n";

if (!open (IN,"$filename")) { print "Error opening $filename ($!).\n"; close(OUT); exit(0); }
print "Opened $filename\n";

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
            #look for primers in the sequence, record p-start-type, p-start-seq, p-end-type, p-end-seq
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
            
            if ($num_sequences % 10000 == 0) { print "$num_sequences\n"; }
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

print "Primers removed finished: $num_sequences sequences processed, $num_sequences_with_primers sequences found with both primers.\n";
print "$num_sequences_with_primers1 sequences found with only R1 primer, $num_sequences_with_primers2 found with only R2 primer.\n";

close(OUT);
close(IN);

sub remove_primers
{
	my $seq = shift;
	my $seq_len = length($seq);
	my $fixed_seq = '';
	#look for p1 at position 0+6 of sequence - allow 8 mismatches
	#look for p2-LH,p2-SH at position 0+6 of sequence - allow 8 mismatches for both
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
		my @start_pos_range = (-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10); 
		foreach my $start_pos (@start_pos_range)
		{ #positions in the sequence to look for primer
			
			my $start_at; my $end_at;
            if($p_id eq 'p2_LH')
			{
			    $start_at = $start_pos;
			    $end_at = $start_pos+25;
			    
			}
			elsif($p_id eq 'p2_SH')
			{
			    $start_at = $start_pos;
			    $end_at = $start_pos+21;
			    
			}
			else #p1
			{
				$start_at = $start_pos;
				$end_at = $start_pos+23;
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
	
	if( ($r1_max_match_p_id eq 'p2_LH' && ($r1_max_match_mismatches <= $max_p2_LH_allowed_mismatches)) ||
        ($r1_max_match_p_id eq 'p2_SH' && ($r1_max_match_mismatches <= $max_p2_SH_allowed_mismatches)) ||
        ($r1_max_match_p_id eq 'p1' && ($r1_max_match_mismatches <= $max_p1_allowed_mismatches)) )
	{ $r1_primer_found = 1; }
	
	#look for r2 primer:
	#look for rev. complement of p1 at position seq_len-6 of sequence - allow 8 mismatches
	#look for rev. complement of p2_LH/SH at position seq_len-6 of sequence - allow 8 mismatches for both
	#see which one matches the best, and remove based on the best match
	
	foreach my $p_id (keys %primers_rc)
	{
		my @start_pos_range = (-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10); 
		foreach my $start_pos (@start_pos_range)
		{ #positions in the sequence to look for adaptor
			
			my $start_at; my $end_at;
			if($p_id eq 'p2_LH')
			{
				$start_at = $seq_len-$start_pos-25; 
				$end_at = $seq_len-$start_pos;
			    
			}
			elsif($p_id eq 'p2_SH')
			{
				$start_at = $seq_len-$start_pos-21; 
				$end_at = $seq_len-$start_pos;
			}
            else #p1
            {
                $start_at = $seq_len-$start_pos-23; 
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
	if( ($r2_max_match_p_id eq 'p2_LH' && ($r2_max_match_mismatches <= $max_p2_LH_allowed_mismatches)) ||
        ($r2_max_match_p_id eq 'p2_SH' && ($r2_max_match_mismatches <= $max_p2_SH_allowed_mismatches)) ||
	    ($r2_max_match_p_id eq 'p1'    && ($r2_max_match_mismatches <= $max_p1_allowed_mismatches)) )
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