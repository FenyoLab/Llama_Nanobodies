#!/usr/bin/perl

# Web interface for viewing Llama Magic results

use strict;
use warnings;
use CGI ':standard';
use Proc::Background;

my $DBLIST_FILE = "db_list.txt";
my $MSLIST_FILE = "ms_list.txt";
my $SETTINGS_FILE = "../settings.txt"; 
my $img_source_plus = '/llama-magic-html/plus.gif';
my $img_source_minus = '/llama-magic-html/minus.gif';
my $img_source_star = '/llama-magic-html/greyplus.gif';
#my $img_source_new = '/llama-magic-html/add_item.png';
my $DEVELOPER_VERSION = 1;
#my $DEVELOPER_VERSION = 0;
my $DEVELOPER_LOGFILE = "cgi_dev_log.txt";
my $BASE_DIR = "/var/www/llama-magic"; #default, changed when settings.txt file is read
my $RESULTS_DIR = "results";

eval #for exception handling
{
    if($DEVELOPER_VERSION) { open(DEVEL_OUT, ">>$BASE_DIR/$DEVELOPER_LOGFILE"); }
    if($DEVELOPER_VERSION) { print DEVEL_OUT "Opened developer log file...\n"; }
    
    my $err = "";
    if($err = read_settings())
    {
		display_error_page("Cannot load settings file: $err");
		if($DEVELOPER_VERSION) { print DEVEL_OUT "Cannot load settings file: $err\n"; }
    }
    else
    {
		if($DEVELOPER_VERSION) { print DEVEL_OUT "Loaded settings file: BASE_DIR = $BASE_DIR\n"; }
		if(!param())
		{#no posted data
	    	display_home_page();
		}
		else
		{
			display_home_page();
		}
    }
};
if ($@)
{
	if($DEVELOPER_VERSION) { print DEVEL_OUT "Exception thrown: $@\n"; }
	display_error_page("$@"); 
}

if($DEVELOPER_VERSION) { close(DEVEL_OUT); }

###########################

sub display_home_page
{
    top_header();
    
    print "<tr class='main_list'>";

    #get list of dbs existing in repository any ms files + 'nanobody candidate results' under it
    open(IN, "$BASE_DIR/$RESULTS_DIR/$DBLIST_FILE");
    my @db_lines = <IN>;
    close(IN);
    
    open(IN, "$BASE_DIR/$RESULTS_DIR/$MSLIST_FILE");
    my @ms_lines = <IN>;
    close(IN);
    
    print "<td>", qq!<br/><b><u>Databases:</u></b> <br /><br />!;
    for(my $i = 0; $i <= $#db_lines; $i++)
    {
		chomp($db_lines[$i]);
		if($db_lines[$i] =~ /^(\d+)\t([^\t]+)$/)
		{
	    	my $id = $1; my $name = $2;

	    	my @ms_list;
	    	#check status, if db translated and digested, create link, else if errors, create msg, else if not done, create disabled 'link' and msg
	    	#my $status = check_status_file("$id/protein");
			my $status = "DONE";
	    	if ($status eq 'DONE')
	    	{
				#get ms searches for this db
				for(my $j = 0; $j <= $#ms_lines; $j++)
				{
					if ($ms_lines[$j] =~ /^$id\t(\d+)\t([^\t]+)$/)
					{
						push(@ms_list, [$1, $2]);
					}
				}
		
				if ($#ms_list >= 0)
				{
					print qq!<img src="$img_source_plus" onclick="ec('text_$id', 'img_$id')" id="img_$id" style="cursor:hand;" alt="+" />!;
				}
				else
				{
					print qq!<img src="$img_source_star" alt="*" />!;
				}

				print qq!$name<br />!;

	    	}
	    	elsif($status eq '') { print qq!<img src="$img_source_star" alt="*" />!, u($name), " &nbsp;(in progress...)", br(); }
	    	else { print qq!<img src="$img_source_star" alt="*" />!,  u($name), " &nbsp;<a href='#' onclick='alert(\"" . $status . "\")'>*Error!</a>", br(), br(); }
	    
	    	#print the ms sub projects
	    	print qq!<br/><div id="text_$id" style="display:none">!;
	    	for(my $j = 0; $j <= $#ms_list; $j++)
	    	{
				my $ms_id = $ms_list[$j][0]; my $ms_name = $ms_list[$j][1];
		
				print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
		
				#check status, if tandem run/candidate list created, create link, else if errors, create msg, else if not done, create disabled 'link' and msg
				#my $status = check_status_file("$id/$ms_id");
				my @files = <$BASE_DIR/$RESULTS_DIR/$id/$ms_id/tandem/results/*.1.cdr_coverage.html>;
				#my $status = "DONE";
				#if ($status eq 'DONE')
				if(scalar @files >= 1)
				{
					#my @files = <$BASE_DIR/$RESULTS_DIR/$id/$ms_id/tandem/results/*.1.cdr_coverage.html>;
					#strip everything but file name....
					$files[0] =~ /[\/\\]([^\/\\]+)$/;
					print a({href=>"../llama-magic/$id/$ms_id/tandem/results/$1", target=>'_blank'}, "$ms_name"), br(), br();
				}
				elsif($status eq '') { print u($ms_name), " &nbsp;(in progress...)", br(), br(); }
				else { print u($ms_name), " &nbsp;<a href='#' onclick='alert(\"" . $status . "\")'>*Error!</a>", br(), br(); }
		
			}
			if($#ms_list == -1) { print br(); }
			print '</div>';
		}
    }
    
    print "</td>";
    print "</tr>";
    display_footer();
}

sub get_db_name
{
    my $id = shift;
    
    open(IN, "$BASE_DIR/$RESULTS_DIR/$DBLIST_FILE");
    my @lines = <IN>;
    close(IN);
    
    for(my $i = 0; $i <= $#lines; $i++)
    {
	chomp($lines[$i]);
	if ($lines[$i] =~ /^$id\t(.+)$/)
	{
	    return $1;
	}
    }
    
    return '';
}

sub display_footer
{
	print '</table>',
	      end_multipart_form(),
	      end_html();
}

sub display_error_page
{
    my $msg = shift;
    top_header(1);
    print "<tr class='main_help'><td>";
    print p($msg);
    print "</td></tr>";
    
    display_footer();
}

sub display_message_page
{
    my $msg = shift;
    top_header(1);
    print "<tr class='main_help'><td>";
    print p($msg),
    a({href=>"../llama-magic-cgi/llama_magic.pl?submit=Home"}, "Back To Home Page");
    print "</td></tr>";
    display_footer();
}

#########################

sub read_settings
{
	open(IN, "$SETTINGS_FILE") || return $!;
	my $found = 0;
	while(<IN>)
	{
		chomp();
		if(/^INSTALL_DIR=(.*)$/) { $BASE_DIR = $1; $found++; }
	}
	close(IN);
	if($found == 1) { return ""; }
	else { return "Information missing from settings file: $SETTINGS_FILE.\n"; }
}

sub top_header
{
    print header(),
	  start_html(-title => 'Llama Magic Viewer',
		     -style => [{ -src=>'../llama-magic-html/smoothness/jquery-ui-1.10.3.custom.min.css' }, #download and fix this
				{ -src=>'../llama-magic-html/main.css' }], 
		     -script => [{ -type=>'javascript', -src=>'../llama-magic-html/jquery-1.9.1.min.js' },
				{ -type=>'javascript', -src=>'../llama-magic-html/jquery-ui-1.10.3.custom.min.js' }, #download and fix this
				{ -type=>'javascript', -src=>'../llama-magic-html/jquery.MultiFile.pack.js' },
				{ -type=>'javascript', -src=>'../llama-magic-html/setup.js' }]);
	  

	  print "<table class='maintable'><tr class='banner' ><td colspan='1'>",
	  qq!<h1><img id="llama" src="/llama-magic-html/llama.jpg">Llama Magic Viewer</h1> </td></tr>!;
}

sub check_status_file
{
    my $dir = shift;

    #try to open file:
    if (-e "$BASE_DIR/$RESULTS_DIR/$dir/status.txt")
    {
	open(IN, "$BASE_DIR/$RESULTS_DIR/$dir/status.txt");
	my @lines = <IN>;
	close(IN);
	my $err_str = '';
	my $done = 0;
	for(my $i = 0; $i <= $#lines; $i++)
	{
	    chomp($lines[$i]);
	    if ($lines[$i] =~ /^ERROR:/)
	    {
		if ($lines[$i] =~ /exited with value/)
		{
		    chomp($lines[$i+1]);
		    if($lines[$i+1] ne ''){ $err_str .= $lines[$i] . ': ' . $lines[$i+1]; }
		    else { $err_str .= $lines[$i] }
		}
		else { $err_str .= $lines[$i]; }
		$err_str .= '\n';
	    }
	    elsif ($lines[$i] eq 'DONE')
	    {
		$done = 1;
		last;
	    }
	}
	if ($done)
	{
	    if ($err_str) { return $err_str; }
	    else { return 'DONE'; }
	}
	else { return ""; }
    }
    else { return ""; }
}
