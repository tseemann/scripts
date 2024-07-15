#!/usr/bin/env perl
use strict;
use Data::Dumper;

my(@Options, $verbose, $output, $force, $dirs, $header, $footer);
setOptions();

my @files = @ARGV ? @ARGV : (<*>);
@files = grep { $_ !~ m/^($output|$header|$footer)$/ } @files;
die "no files found in current dir or in argument list" unless @files;
die "$output already exists" if (!$force and -e $output);

print STDERR "Writing: $output ...\n";
open my $OUT, '>', $output;
select $OUT;

if (-r $header) {
  print "<!--#include virtual='$header'-->\n";
}
else {
  print "<h1>\n", qx(basename `pwd`), "</h1>\n";
}

print "<ul>\n";
print map { "<li><a href='$_'>$_</a>\n" } @files;
print "</ul>\n";

if (-r $footer) {
  print "<!--#include virtual='$footer'-->\n";
}
else {
  print "<p>Published by ", $ENV{USER}, " at ", qx(date);
}

my $host = qx(hostname); chomp $host;
my $uid = qx(whoami); chomp $uid;
my $dir = qx(pwd); chomp $dir; $dir =~ s/^.*public_html\/?//;
my $url = "http://$host/~$uid/$dir/$output";
print STDERR "URL: $url\n";

#print "<p><small><a href='$url'><tt>$url</tt></a></small>\n";
close $OUT;
print STDERR "Done.\n";


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"output=s",  VAR=>\$output, DEFAULT=>'index.shtml', DESC=>"Output file"},
#    {OPT=>"d|dirs!",  VAR=>\$dirs, DEFAULT=>0, DESC=>"Include folder links too?"},
    {OPT=>"force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force overwrite"},
    {OPT=>"header=s",  VAR=>\$header, DEFAULT=>'header.inc', DESC=>"SSI-include this at top"},
    {OPT=>"footer=s",  VAR=>\$footer, DEFAULT=>'footer.inc', DESC=>"SSI-include this at bottom"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
