#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use vars qw();
getopts('');

my $SWITCH = 0.0001;
my $ERROR = 0.001;

my $usage = "
usage: $0 [optons] <datafile.gz>
options:
  -g <float> probabilty of genotype switch  [$SWITCH]
  -n <float> probabilty of nucleotide error [$ERROR]
";

die $usage unless @ARGV == 1;
my ($file) = @ARGV;

open(my $fh, "gunzip -c $file |") or die;

my $hs = <$fh>;
my ($strains) = $hs =~ /strains: (\d+)/;
my @strain;
for (my $i = 0; $i < $strains; $i++) {
	my $s = <$fh>;
	chomp $s;
	push @strain, $s;
}

my $hc = <$fh>;
my ($chroms) = $hc =~ /chroms: (\d+)/;
my @chrom;
for (my $i = 0; $i < $chroms; $i++) {
	my $s = <$fh>;
	chomp $s;
	push @chrom, $s;
}

my %genotype;
my $hg = <$fh>;
while (<$fh>) {
	chomp;
	my ($chr, $var) = split;
	push @{$genotype{$chr}}, $var;
}

my $strain = int rand(@strain);
my @history = "$strain";
foreach my $chr (@chrom) {
	print "$chr: ";
	for (my $i = 0; $i < @{$genotype{$chr}}; $i++) {
		if (rand(1) < $SWITCH) {
			my $newstrain = int rand(@strain);
			if ($newstrain != $strain) {
				$strain = $newstrain;
				push @history, "$strain:$chr:$i"
			}
		}
		if (rand(1) < $ERROR) {
			print "n";
		} else {
			print substr($genotype{$chr}[$i], $strain, 1);
		}
	}
	print "\n";
}
print "# @history\n";





