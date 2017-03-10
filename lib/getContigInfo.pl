#! /usr/bin/perl -w
####################################################################
#
#	Get contigs info
#
#	Author: Ching-Hsin Liu, <crocodilepp@gmail.com>, 2010.03.10
#
#	Input:
#		1. contigs.fa (fasta file)
#		2. minimum contig length (integer)
#
####################################################################
use strict;

if($#ARGV + 1 < 2) {
	print  <<EOF;
Usage:
	$0 <contigs.fa> <minimum contig length>

EOF
	exit;
}

my $input = $ARGV[0];
my $minLength = $ARGV[1];

open IN, "$input" or die;

my %contig_length = ();
my $contig_id = "";

while(<IN>) {
	if(/^>(.+)/) {
		$contig_id = $1;
	} else {
		chomp $_;
		$contig_length{$contig_id} += length $_;
	}
}

close IN;

my @length_order = sort by_number values(%contig_length);
my $contig_counts = 0;
my $total_read_length = 0;

my $hundred = 0;
my $large = 0;
my $thousand = 0;
my $fivehund = 0;

foreach my $l(@length_order) {
	$l >= $minLength and $total_read_length += $l and $contig_counts++;
        $l >= $minLength and print $l."\n";
	$l >= 100 and $hundred++;
	$l >= 200 and $large++;
	$l >= 500 and $fivehund++;
	$l >= 1000 and $thousand++;
}

my $N50 = 0;
my $N50_ctg = 0;

for(my $i = 0;$i <= $#length_order;$i++) {
	$length_order[$i] < $minLength and next;
	$N50 += $length_order[$i];
	if($N50 > ($total_read_length / 2)) {
		$N50_ctg = $i - 1;
		last;
	}
}

print "\nminimum length: $length_order[0]\n";
print "maximum length: $length_order[$#length_order]\n";
print "2nd long contig: $length_order[$#length_order-1]\n";
print "3rd long contig: $length_order[$#length_order-2]\n";
print "total length: $total_read_length\n";
print "avg. length: ".round ($total_read_length/$contig_counts,1)."\n";
print "N50: $length_order[$N50_ctg]\n\n";

print "contig counts:\n";
print " > $minLength bp: $contig_counts\n";
print "( >  100 bp: $hundred)\n";
print "( >  200 bp: $large)\n";
print "( >  500 bp: $fivehund)\n";
print "( > 1000 bp: $thousand)\n\n";

#exec "~keepp/perl_scripts/base_content.pl $ARGV[0]";

sub by_number {
	$a<=>$b;
}
sub round {
	my $val = shift;
	my $col = shift;
	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;
	return int($val * $r + $a) / $r;
}

=cut
use Data::Dumper;
print Dumper(\@length_order);
print Dumper(\%contig_length);
