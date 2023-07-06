#!/usr/local/bin/perl

use File::Copy;
use warnings;
use strict;

#input params are mgf_dir, root_dir

#open parameters file
# index     mgf_name    db_id   parent_error    parent_units    frag_error  frag_units

#for reach row, read settings, name title to match  mgf name with errors in name
#create folder structure under root_dir
my $TEST=0;

my $params_file = "input_params.txt";
my $ms_list_file  = "ms_list.txt";

my $base_dir;
if ($ARGV[0] && $ARGV[0]=~/\w/) { $base_dir=$ARGV[0];}  
else { $base_dir = "/gpfs/data/proteomics/projects/NCDIR/Llama/pipeline"; }
#else { $base_dir = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/data_and_results/llama_pipeline_test"; }

my $pipeline_dir = "pipeline";
my $results_dir = "results";

open(IN, "$base_dir/$pipeline_dir/$params_file") || die("Failed to open file '$base_dir/$pipeline_dir/$params_file': $!");
my @lines = <IN>;
close(IN);
for(my $i = 0; $i <= $#lines; $i++)
{ 
    if ($lines[$i] =~ /(\d+)\t([^\t]*)\t([^\t]+)\t([^\t]*)\t(\d+)\t([\d\.]+)\t(\w+)\t([\d\.]+)\t(\w+)\t([\d\.]*)\t(\w*)\t([\d\.]*)\t(\w*)\t(trypsin|chymotrypsin|chymotrypsin2)\t(\w*)/)
    {
        my $i = $1;
        my $web_name = $2;
        my $mgf_name = $3;
        my $mgf_name2 = $4;
        my $db_id = $5;
        my $parent_err = $6;
        my $parent_err_units = $7;
        my $frag_err = $8;
        my $frag_err_units = $9;
        my $parent_err2 = $10;
        my $parent_err_units2 = $11;
        my $frag_err2 = $12;
        my $frag_err_units2 = $13;
        my $enzyme = $14;
        my $enzyme2 = $15;
        
        #read list of db searches, get next search id for the db id selected
        open(IN, "$base_dir/$results_dir/$ms_list_file") || die("Failed to open file: $!");
        my @lines = <IN>;
        close(IN);
        
        my $new_ms_id = 1;
        for(my $i = 0; $i <= $#lines; $i++)
        {
            if ($lines[$i] =~ /(\d+)\t(\d+)\t(.+)/)
            {
                if ($1 eq $db_id and $2 >= $new_ms_id) { $new_ms_id = $2+1; }
            }
        }
	
        #add 1000 since it's the pipeline - we don't want it to collide with any searches created manually
        if($new_ms_id < 1000) { $new_ms_id = $new_ms_id + 1000; }
            
        #create new dir structure for the search
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id") || die("Failed to make dir '$base_dir/$results_dir/$db_id/$new_ms_id': $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/mgf") || die("Failed to make dir '$base_dir/$results_dir/$db_id/$new_ms_id/mgf': $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/tandem") || die("Failed to make dir: '$base_dir/$results_dir/$db_id/$new_ms_id/tandem' $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/tandem/results") || die("Failed to make dir '$base_dir/$results_dir/$db_id/$new_ms_id/tandem/results': $!");
            
        #copy mgf file:
        copy("$base_dir/$pipeline_dir/$mgf_name","$base_dir/$results_dir/$db_id/$new_ms_id/mgf/$mgf_name") ||
            die("Failed to copy file from '$base_dir/$pipeline_dir/$mgf_name' to '$base_dir/$results_dir/$db_id/$new_ms_id/mgf/$mgf_name': $!");
        
        #copy mgf file 2, if it exists:
        $mgf_name2 =~ s/^\s+//;
        $mgf_name2 =~ s/\s+$//; 
        if($mgf_name2)
        {
            copy("$base_dir/$pipeline_dir/$mgf_name2","$base_dir/$results_dir/$db_id/$new_ms_id/mgf/$mgf_name2") ||
                die("Failed to copy file from '$base_dir/$pipeline_dir/$mgf_name2' to '$base_dir/$results_dir/$db_id/$new_ms_id/mgf/$mgf_name2': $!");
        }
        
        $web_name =~ s/^\s+//;
        $web_name =~ s/\s+$//;
        if (!$web_name)
        {
            $web_name = "$i.$mgf_name-$parent_err-$parent_err_units-$frag_err-$frag_err_units-$enzyme";
        }
        open(OUT, ">>$base_dir/$results_dir/$ms_list_file") || die("Failed to open file: $!");
        print OUT "$db_id\t$new_ms_id\t$web_name\n"; 
        close(OUT);
        
        #read db params to see if we should use primers and tail sequences
        my $USE_PRIMERS = 0; 
        my $ADD_TAIL_SEQUENCES = 0;
        my $NEW_PRIMERS = 0;
        if(open(IN, "<$base_dir/$results_dir/$db_id/db_params.txt"))
        { 
            while (<IN>)
            {
                if (/USE_PRIMERS=(\d)/) { $USE_PRIMERS = $1; }
                if (/ADD_TAIL_SEQUENCES=(\d)/) { $ADD_TAIL_SEQUENCES = $1; }
                if (/NEW_PRIMERS=(\d)/) { $NEW_PRIMERS = $1; }
            }
            close(IN);
        }
        else { die("Failed to open file: $!"); }
        
        #start qsub
        my $run_str = "";
        if ($mgf_name2)
        {
            $run_str=qq(sbatch -p cpu_short run_LM_script.sh "$base_dir/$results_dir/$db_id/$new_ms_id" "$base_dir/$results_dir/$db_id" "$parent_err" "$frag_err" "$parent_err_units" "$frag_err_units" "$parent_err2" "$frag_err2" "$parent_err_units2" "$frag_err_units2" "$enzyme" "$enzyme2" "/$results_dir/$db_id/$new_ms_id" "1" "$USE_PRIMERS" "$ADD_TAIL_SEQUENCES" "0" "$NEW_PRIMERS" "$mgf_name" "$mgf_name2");
        }
        else
        {
            $run_str=qq(sbatch -p cpu_short run_LM_script.sh "$base_dir/$results_dir/$db_id/$new_ms_id" "$base_dir/$results_dir/$db_id" "$parent_err" "$frag_err" "$parent_err_units" "$frag_err_units" "-" "-" "-" "-" "$enzyme" "-" "/$results_dir/$db_id/$new_ms_id" "1" "$USE_PRIMERS" "$ADD_TAIL_SEQUENCES" "0" "$NEW_PRIMERS" "$mgf_name" "-");
        }
        
        print "Running job: \n";
        print $run_str;
        print "\n";
        
        if (!$TEST)
        {
            system("$run_str");
        }
        
        
        #system("sbatch -p cpu_short run_LM_script.sh \"$base_dir/$results_dir/$db_id/$new_ms_id\" \"$base_dir/$results_dir/$db_id\" \"$parent_err\" \"$frag_err\" \"$parent_err_units\" \"$frag_err_units\"
        #\"$enzyme\" \"/$results_dir/$db_id/$new_ms_id\" \"1\" \"$USE_PRIMERS\" \"$ADD_TAIL_SEQUENCES\" \"0\" \"$NEW_PRIMERS\"");
                    
    }
}


