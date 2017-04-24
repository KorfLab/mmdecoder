#!/usr/bin/perl
use strict; use warnings;

# strains
my $header = <>;
$header =~ s/\s+//g;
my ($id, $chr, $b37, $b38, @strain) = split(/,/, $header);
#print "strains: ";
#print scalar(@strain), "\n";
#foreach my $strain (@strain) {print "$strain\n"}

# read
my %g;
my $count = 0;
while (<>) {
	s/\s+//g; # possible windows style line endings
	my ($id, $chr, $b37, $b38, @var) = split(",", $_);
	push @{$g{$chr}}, join("", @var);
	$count++;
}

# chromosomes
no warnings;
my @chr = sort {$a <=> $b or $a cmp $b} keys %g;
use warnings;
#print "chroms: ";
#print scalar keys %g, "\n";
#foreach my $chr (@chr) {print "$chr\n"}

# sequences
my @seq;
foreach my $chr (@chr) {
	foreach my $gt (@{$g{$chr}}) {
		for (my $i = 0; $i < length($gt); $i++) {
			$seq[$i]{$chr} .= substr($gt, $i, 1);	
		}
	}
}

# output
for (my $i = 0; $i < @seq; $i++) {
	for (my $j = 0; $j < @chr; $j++) {
		print "$strain[$i] $chr[$j] $seq[$i]{$chr[$j]}\n";
	}
}
