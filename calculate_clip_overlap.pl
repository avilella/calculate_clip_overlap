#!/usr/bin/perl
my $version = '1.0.1';
use strict;
use warnings;
use Getopt::Long;

my $self = bless {};

my $inputfile;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
my $iter = 100000;
my $genome_size = 3000000000;
my $range = "75:100:125:150:250";
my $tag = 'clip_overlap';
my $obs = '2x75';
my $print_version;

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'iter:s'              => \$iter,
           'genome_size:s'       => \$genome_size,
	   'range:s'             => \$range,
	   'tag:s'               => \$tag,
	   'obs|observed:s'      => \$obs,
           'debug'    => \$debug,
           'verbose'  => \$verbose,
           'version'  => \$print_version,
	   'simulate' => \$simulate,
          );

if ($print_version) {
  print "$version\n";
  exit(0);
}

$cmd = "samtools view $inputfile |";
my $stats;
my @ranges = split(":",$range);
open FILE, $cmd or die $!;
my $count = 1;
my $avg_tlen = 0; my $num_tlen = 0;
while (<FILE>) {
  my $line = $_;
  my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual) = split(" ",$line);
  $avg_tlen += abs($tlen); $num_tlen++;
  foreach my $cycles (@ranges) {
    my $gap = $cycles - abs($tlen/2);
    my $mapped = 0;
    my @c = $cigar =~ /(\d+)M/g;
    map { $mapped += $_ } @c;
    $stats->{mapped_nuc}{$cycles} += $mapped;
    if ($gap > 0) {
      # We are losing bases
      $stats->{lost_nuc}{$cycles} += $gap;
    }
    $stats->{total}{$cycles} += $cycles;
  }
  $self->summary if ($verbose && $count >= $iter && 0 == $count % $iter); $count++;
}

my $fh;
my $outfile = "$inputfile.$tag.txt";
open ($fh, ">" . "$outfile") or die $!;
$self->summary(1,$fh);
close $fh;

print "# Output file:\n";
print "$outfile\n";

1;

########################################

sub summary {
  my $self = shift;
  my $stdout  = shift;
  my $fh  = shift;

  my $outh = *STDERR; $outh = $fh if ($stdout);

  my $this_avg_tlen = sprintf("%.02f",$avg_tlen/$num_tlen);
  print $outh "#inputfile=$inputfile\n";
  print $outh "#avg_tlen=$this_avg_tlen\n";
  print $outh "SEQ\tread_PE\tscov\ttotal_nuc\tmapped_nuc\tlost_nuc\teffic_pc\n";
  foreach my $cycles (@ranges) {
    my $str = $count . "\t" .
      sprintf("%.05f",($stats->{total}{$cycles})/$genome_size) . "\t" .
      $stats->{total}{$cycles} . "\t" .
      $stats->{mapped_nuc}{$cycles} . "\t" .
      $stats->{lost_nuc}{$cycles} . "\t" .
      sprintf("%.02f",100-((100*($stats->{lost_nuc}{$cycles}))/$stats->{total}{$cycles})) . "\n";
    print $outh "obs" . "\t" . $str if ($obs eq "2x$cycles");
    print $outh "2x$cycles" . "\t" . $str;
  }

  return;
}
