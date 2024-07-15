#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use File::Spec;

my(@Options, $verbose, $runid, $lane, $rootdir, $force, $to_email, $from_email);
setOptions();

#bail("Invalid --lane [1-8]") unless defined $lane and $lane >= 1 and $lane <= 8;
#bail("Invalid --runid [YYMMDD]") unless defined $runid and $runid =~ m/^\d{6}$/;

bail("Invalid customer --email") unless defined $to_email and $to_email =~ m/\w\@\w/;
bail("Can't write to --rootdir $rootdir") unless -w $rootdir;

my @files = sort @ARGV;
bail("Please specify files to put on web") unless @files;
for my $f (@files) {
  bail("Can't read '$f'") unless -r $f;
  print STDERR "Including: $f\n";
}
my $md5 = md5_hex(@files, $$, localtime);
#my $md5 = md5_hex(@files, $$.localtime);
#my $md5 = md5_hex(@files);
print STDERR "MD5: $md5\n";
$md5 = substr $md5, 0, 16;
my $dir = File::Spec->catdir($rootdir,$md5);
print STDERR "Making dir: $dir\n";
bail("$dir already exists") if -r $dir;
mkdir $dir or bail("unable to create dir: $dir");

# read the template
print STDERR "Reading template...\n";
my $page = do {local $/;<DATA>}; # slurp it all in

my $insert = '';
$insert .=  "<h2>Customer</h2><p>$to_email\n";
$insert .=  "<h2>Files</h2>\n";
#$insert .=  "<p>You have ".(0+@files)." sequence files available to download:\n";
$insert .=  "<ul>\n";
for my $f (@files) {
  my $size = -s File::Spec->catfile($dir, $f);
  bail("Size of $f was undefined. Not sure why. Exiting.") unless defined $size;
  $size = int($size/1E6);
  $insert .=  "<li><a href='$f'><tt>$f</tt></a> ($size Mb)\n";
}
$insert .=  "</ul>\n";

$insert .=  "<P>Email queries to <a href='mailto:$from_email'>$from_email</a>\n";

$insert .=  "<h2>Upload date</h2><p>".qx(date);

my $html = File::Spec->catfile($dir, 'index.html');
print STDERR "Writing: $html\n";
open HTML, '>', $html;
$page =~ s/INSERT_HERE/$insert/;
print HTML $page;
close HTML;

my($hostname) = qx(hostname); chomp $hostname;
#my($uname) = qx(whoami); chomp $uname;
print STDERR "URL: http://$hostname/results/$md5/\n";

sub bail {
  print STDERR "ERROR: @_\n";
  exit -1;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
#    {OPT=>"r|runid=s",  VAR=>\$runid, DEFAULT=>'090729', DESC=>"Run date in YYMMDD format"},
#    {OPT=>"l|lane=i",  VAR=>\$lane, DEFAULT=>0, DESC=>"Lane"},
#    {OPT=>"r|rootdir=s",  VAR=>\$rootdir, DEFAULT=>$ENV{HOME}.'/public_html/', DESC=>"Root folder"},
    {OPT=>"r|rootdir=s",  VAR=>\$rootdir, DEFAULT=>$ENV{HOME}.'/public_html/results/', DESC=>"Root folder"},
#    {OPT=>"f|force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force overwrite"},
#    {OPT=>"from=s",  VAR=>\$from_email, DEFAULT=>'enquiries@dna.med.monash.edu.au', DESC=>"From email"},
    {OPT=>"from=s",  VAR=>\$from_email, DEFAULT=>'scott.coutts@monash.edu', DESC=>"From email"},
    {OPT=>"e|email=s",  VAR=>\$to_email, DEFAULT=>'', DESC=>"Email address of customer"},
  );

#  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 <file1.fq.gz file2.fq.gz ...>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<!-- DW6 -->
<head>
<!-- Copyright 2005 Macromedia, Inc. All rights reserved. -->
<title>Micromon - High Throughput Sequencing Service</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link rel="icon" type="image/jpg" href="/MicromonLogo-Icon.jpg"/>
<link rel="stylesheet" href="/mm_training.css" type="text/css" />
<style type="text/css">
<!--
.style2 {
	color: #003399
}
-->
</style>
</head>
<body bgcolor="#64748B">

<!-- Start of StatCounter Code -->
<script type="text/javascript">
var sc_project=4872716; 
var sc_invisible=1; 
var sc_partition=57; 
var sc_click_stat=1; 
var sc_security="54f4760f"; 
</script>

<script type="text/javascript"
src="http://www.statcounter.com/counter/counter.js"></script><noscript><div
class="statcounter"><a title="profile counter myspace"
href="http://www.statcounter.com/myspace/"
target="_blank"><img class="statcounter"
src="http://c.statcounter.com/4872716/0/54f4760f/1/"
alt="profile counter myspace" ></a></div></noscript>
<!-- End of StatCounter Code -->

<p>&nbsp;</p>
<table width="1020" border="0" align="center" bordercolor="#E0F4FE" bgcolor="#E0F4FE">
  <tr>
    <td><div align="center">
      <table width="100%" border="0" cellpadding="0" cellspacing="0" bordercolor="#000099">
        <tr bgcolor="#26354A">
          <td height="70" colspan="2" align="left" valign="middle" nowrap="nowrap" bgcolor="#FFFFFF"><h1><img src="../mm_spacer.gif" alt="G" width="15" height="1" border="0" /><img src="http://dna.med.monash.edu.au/MicromonLogo-Small-White.jpg" alt="MICROMON" width="211" height="80" /></h1></td>
          <td height="70" align="left" valign="middle" nowrap="nowrap" bgcolor="#FFFFFF">&nbsp;</td>
          <td height="70" align="right" valign="middle" nowrap="nowrap" bgcolor="#FFFFFF" class="pageName"><strong>HIGH THROUGHPUT SEQUENCING SERVICE</strong></td>
          <td bgcolor="#FFFFFF">&nbsp;</td>
        </tr>

        <tr bgcolor="#D3DCE6">
          <td colspan="5"><img src="../mm_spacer.gif" alt="E" width="1" height="1" border="0" /></td>
        </tr>
        <tr bgcolor="#75BDEB">
          <td width="229" nowrap="nowrap">&nbsp;</td>
          <td height="24" colspan="3"><table border="0" cellpadding="0" cellspacing="0" id="navigation">
            <tr>
              <td height="34" align="center" nowrap="nowrap" class="navText"><a href="/index.html">HOME</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/sequencing_services.html">SEQUENCING<br />
                SERVICES</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/informatics_services.html">INFORMATICS<br />
                SERVICES</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/sample.html">SAMPLE <br />
                SUBMISSION</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/pricing.html" class="sidebarText style2">PRICING</a></td>
              <td class="navText" align="center" nowrap="nowrap"><a href="/news.html">NEWS</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/resources.html">RESOURCES</a></td>
              <td align="center" nowrap="nowrap" bgcolor="#FF9900" class="navText"><a href="/results.html">RESULTS</a></td>
              <td class="navText" align="center" nowrap="nowrap"><a href="/faq.html">F.A.Q.</a></td>
              <td align="center" nowrap="nowrap" class="navText"><a href="/contact.html">CONTACT<br />
                US</a></td>
            </tr>
          </table></td>
          <td width="3">&nbsp;</td>
        </tr>
        <tr bordercolor="#0000FF" bgcolor="#D3DCE6">
          <td colspan="5"><img src="../mm_spacer.gif" alt="D" width="1" height="1" border="0" /></td>
        </tr>
        <tr bgcolor="#FF6600">
          <td colspan="5" bgcolor="#FFFFFF"><img src="../mm_spacer.gif" alt="C" width="1" height="4" border="0" /></td>
        </tr>
        <tr bgcolor="#D3DCE6">
          <td colspan="2" valign="top" bgcolor="#B1D9F4"><img src="http://dna.med.monash.edu.au/sequence.jpg" alt="" width="230" height="286" /><br />
              <table border="0" cellspacing="0" cellpadding="0" width="230">
                <tr>
                  <td width="230" class="sidebarText" id="padding"><p><br />                  
                  </p>
                  </td>
                </tr>
            </table></td>
          <td width="50" valign="top" bgcolor="#EAEFFA"><img src="../mm_spacer.gif" alt="B" width="50" height="1" border="0" /></td>
          <td width="778" valign="top" bgcolor="#EAEFFA"><br />
              <br />
              <table border="0" cellspacing="0" cellpadding="0" width="90%">
                <tr>
                  <td colspan="2" class="pageName"><p align="left">Results                  </p>                  </td>
                </tr>
                <tr>
                  <td width="7%" class="bodyText"><div align="justify"></div></td>
                  <td width="93%" class="bodyText"><p align="justify">INSERT_HERE</p>
                  <p align="justify">&nbsp;</p></td>
                </tr>
             </table>
                     <br />
            &nbsp;<br />          </td>
          <td width="3" bgcolor="#EAEFFA">&nbsp;</td>
        </tr>
        <tr bgcolor="#D3DCE6">
          <td colspan="5" bgcolor="#FFFFFF"><table width="100%" border="0">
              <tr>
                <td width="37%" height="24">&nbsp;</td>
                <td width="30%">&nbsp;</td>
                <td width="33%">&nbsp;</td>
              </tr>
            </table></td>
        </tr>
        <tr bgcolor="#D3DCE6">
          <td colspan="5"><img src="../mm_spacer.gif" alt="A" width="1" height="1" border="0" /></td>
        </tr>
      </table>
    </div></td>
  </tr>
</table>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p></p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
</body>
</html>
