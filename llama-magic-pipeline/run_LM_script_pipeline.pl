#!/usr/local/bin/perl

#runs the scripts sent in by command line arguments for pipeline
#search_and_map_scripts (runs xtandem (command line/xml interface), parse_xtandem.pl, filter_tandem_results.pl, and map_peptides_to_proteins.pl)

# "search_and_map_scripts" "C:/Code/NCDIR/Llama/results/64/31" "C:/Code/NCDIR/Llama/results/64" "10" "0.4" "ppm" "Daltons" "/results/64/31" 1 1 1 0

use warnings;
use strict;

my $TEST=0;

my $expect_threshold = '0.1';

my $results_dir = $ARGV[0];
my $taxon_dir = $ARGV[1];
my $parent_err = $ARGV[2];
my $fragment_err = $ARGV[3];
my $parent_err_units = $ARGV[4];
my $frag_err_units = $ARGV[5];
my $parent_err2 = $ARGV[6];
my $fragment_err2 = $ARGV[7];
my $parent_err_units2 = $ARGV[8];
my $frag_err_units2 = $ARGV[9];
my $enzyme = $ARGV[10];
my $enzyme2 = $ARGV[11];
my $cgi_results_dir = $ARGV[12];
my $show_score = $ARGV[13];
my $use_primers = $ARGV[14];
my $add_tails = $ARGV[15];
my $use_constants = $ARGV[16];
my $new_primers = $ARGV[17];
my $mgf_name = $ARGV[18];
my $mgf_name2 = $ARGV[19];

#create status txt file
open(OUT, ">$results_dir/status.txt") or die "Failed to create status file: $results_dir/status.txt";
print OUT "$results_dir\n";
print OUT "$taxon_dir\n";
print OUT "$parent_err\n";
print OUT "$fragment_err\n";
print OUT "$parent_err_units\n";
print OUT "$frag_err_units\n";
print OUT "$parent_err2\n";
print OUT "$fragment_err2\n";
print OUT "$parent_err_units2\n";
print OUT "$frag_err_units2\n";
print OUT "$enzyme\n";
print OUT "$enzyme2\n";
print OUT "$cgi_results_dir\n";
print OUT "$show_score\n";
print OUT "$use_primers\n";
print OUT "$add_tails\n";
print OUT "$use_constants\n";
print OUT "$new_primers\n";
print OUT "$mgf_name\n";
print OUT "$mgf_name2\n";

close(OUT);
open(OUT, ">>$results_dir/status.txt") or die "Failed to edit status file: $results_dir/status.txt";


#get enzyme for search - indicates which pre-digested file to send to XTandem
my $num_xtandem_runs = 1;

my $file_suffix="";
if ($enzyme eq "trypsin")
{ $file_suffix = ""; }
elsif($enzyme eq "chymotrypsin")
{ $file_suffix = "_chymotrypsin"; }
elsif($enzyme eq "chymotrypsin2")
{ $file_suffix = "_chymotrypsin2"; }
else
{
    print OUT "ERROR: Unrecognized enzyme_1 type: '$enzyme'. Aborting run.";
    close_and_exit(0);
}

my $file_suffix2="";
if ($enzyme2 eq "trypsin")
{ $file_suffix2 = ""; $num_xtandem_runs=2; }
elsif($enzyme2 eq "chymotrypsin")
{ $file_suffix2 = "_chymotrypsin"; $num_xtandem_runs=2; }
elsif($enzyme2 eq "chymotrypsin2")
{ $file_suffix2 = "_chymotrypsin2"; $num_xtandem_runs=2; }
elsif($enzyme2 ne "-")
{ print OUT "ERROR: Unrecognized enzyme_2 type: '$enzyme2'. Run will continue with only enzyme_1."; }

# (1) run xtandem for the mgf/db, run also for the 2nd mgf/db if given
my $err = run_xtandem($enzyme, $file_suffix, $fragment_err, $parent_err, $frag_err_units, $parent_err_units, $mgf_name, "input1.xml", "output1.xml");
if ($err) { close_and_exit(0); }

if ($num_xtandem_runs == 2)
{
    $err = run_xtandem($enzyme2, $file_suffix2, $fragment_err2, $parent_err2, $frag_err_units2, $parent_err_units2, $mgf_name2, "input2.xml", "output2.xml");
    if ($err) { close_and_exit(0); }
}

# (2) parse xtandem to tab separated txt file, extract relevant columns
# the perl script called here will parse each xml file in the input directory
my $run_str = qq!"perl" "parse_xtandem_llama.pl" "$results_dir/tandem/results" "$expect_threshold"!;
print OUT "Calling: $run_str\n";

if (!$TEST)
{
    my $cmd_out = `$run_str 2>&1`;
    if ( $? == -1 )
    {
        print OUT "ERROR: command failed (parse_xtandem_llama.pl): $!\n";
        close_and_exit(0);
    }
    elsif($? >> 8 != 0) # exit value not 0, indicates error...
    {
        printf OUT "ERROR: command (parse_xtandem_llama.pl) exited with value %d\n", $? >> 8;
        print OUT "$cmd_out\n";
        close_and_exit(0);
    }
}

# (3) run map_peptides_to_proteins.pl (main LM algorithm)
my $pep_file = "$results_dir/tandem/results/output1.xml.peptide_list.$expect_threshold.txt";
my $index_file_str = "$taxon_dir/protein/protein_peptides" . "$file_suffix.fasta";
my $tandem_xml_file = "$cgi_results_dir/tandem/results/output1.xml";

my $pep_file2 = "";
my $index_file_str2 = "";
my $tandem_xml_file2 = "";
if ($num_xtandem_runs==2)
{
    $pep_file2="$results_dir/tandem/results/output2.xml.peptide_list.$expect_threshold.txt";
    $index_file_str2="$taxon_dir/protein/protein_peptides" . "$file_suffix2.fasta";
    $tandem_xml_file2="$cgi_results_dir/tandem/results/output2.xml";
}

$run_str = qq!"perl" "map_peptides_to_proteins.pl" "$pep_file" "$pep_file2" "$taxon_dir/protein/longest_nr.fasta" "$tandem_xml_file" "$tandem_xml_file2" "$index_file_str" "$index_file_str2" "$show_score" "$use_primers" "$add_tails" "$new_primers"!;
print OUT "Calling: $run_str\n";
if (!$TEST)
{
    my $cmd_out = `$run_str 2>&1`;
    if ( $? == -1 )
    {
        print OUT "ERROR: command failed (map_peptides_to_proteins.pl): $!\n";
        close_and_exit(0);
    }
    elsif($? >> 8 != 0) # exit value not 0, indicates error...
    {
        printf OUT "ERROR: command (map_peptides_to_proteins.pl) exited with value %d\n", $? >> 8;
        print OUT "$cmd_out\n";
        close_and_exit(0);
    }
}

close_and_exit(1);

###### SUBROUTINES ########
sub close_and_exit
{
    my $ret_code = shift;
    
    print OUT "DONE\n";
    close(OUT);
    
    exit($ret_code);
}

sub run_xtandem
{
    my $cur_enzyme = shift;
    my $cur_file_suffix = shift;
    my $cur_fragment_err = shift;
    my $cur_parent_err = shift;
    my $cur_frag_err_units = shift;
    my $cur_parent_err_units = shift;
    my $cur_mgf_name = shift;
    my $input_file_name = shift;
    my $output_file_name = shift;
    
    # (1) create xml input file for XTandem
    if(open(XML_OUT, ">$results_dir/tandem/$input_file_name"))
    {
        print XML_OUT <<XMLTEXT;
<?xml version="1.0"?>
    <bioml>
        <note type="input" label="spectrum, fragment monoisotopic mass error">$cur_fragment_err</note>
        <note type="input" label="spectrum, parent monoisotopic mass error plus">$cur_parent_err</note>
        <note type="input" label="spectrum, parent monoisotopic mass error minus">$cur_parent_err</note>
        <note type="input" label="spectrum, fragment monoisotopic mass error units">$cur_frag_err_units</note>
        <note type="input" label="spectrum, parent monoisotopic mass error units">$cur_parent_err_units</note>
        <note type="input" label="list path, default parameters">/gpfs/data/proteomics/projects/NCDIR/Llama/pipeline/results/default_input.xml</note>
XMLTEXT
    
        print XML_OUT qq(\n\t<note type="input" label="list path, taxonomy information">$taxon_dir/taxonomy);
        print XML_OUT qq($cur_file_suffix);
        print XML_OUT qq(.xml</note>);
      
        if($use_constants)
        {
            print XML_OUT qq(\n\t<note type="input" label="protein, taxon">llama_constant</note>);
        }
        else
        {
            print XML_OUT qq(\n\t<note type="input" label="protein, taxon">llama</note>);
        }
        
        print XML_OUT qq(\n\t<note type="input" label="spectrum, path">$results_dir/mgf/$cur_mgf_name</note>\n);
        
        print XML_OUT <<XMLTEXT;
        <note type="input" label="output, path">$results_dir/tandem/results/$output_file_name</note>
        <note type="input" label="output, results">valid</note>
    </bioml>
XMLTEXT
        close(XML_OUT);
    }
    else
    {
        print OUT "ERROR: Failed to create xtandem input file: $results_dir/tandem/$input_file_name";
        return 1;
    }
    
    # (2) run xtandem with uploaded mgf file
    my $run_str = qq!"tandem.exe" "$results_dir/tandem/$input_file_name"!;
    print OUT "Calling: $run_str\n";
    
    if (!$TEST)
    {
        my $cmd_out = `$run_str 2>&1`;
        if ( $? == -1 )
        {
            print OUT "ERROR: command failed (tandem.exe): $!\n";
            return 1;
        }
        elsif($? >> 8 != 0) # exit value not 0, indicates error...
        {
            printf OUT "ERROR: command (tandem.exe) exited with value %d\n", $? >> 8;
            print OUT "$cmd_out\n";
            return 1;
        } 
    }
    return 0;
}



