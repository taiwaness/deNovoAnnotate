#!/usr/bin/perl -w
use strict;
if (!@ARGV){ print "Usage: $0 -h\n"; exit;}
use Getopt::Std;

sub usage {
	print <<EOF;
Usage: get_besthit.pl [OPTION]
  -i  input file (Blast+ outfrm=6 or Blast -m=8 , default:blastout)
  -o  output file (default:outfile)
  -e  E-value threshold (default:0.00001) (example: 1e-05 , 0.00001)
  -d  identity(0~1 ex. 0.85)
  -h  print this help, then exit

Report bugs to <tinin\@mars.csie.ntu.edu.tw>.
EOF
	exit;
}
my %opt = (); &getopts( "i:o:he:d:", \%opt );

# default values of parameters
my $in = defined $opt{i} ? $opt{i} : "blastout"; # blastout file
my $out = defined $opt{o} ? $opt{o} : "outfile"; # output file
my $eth = defined $opt{e} ? $opt{e} : "0.00001"; # E-value threshold
($eth !~/^\d+\.\d+$/ and $eth !~/^\d+(\.\d+)?e\-?\d+$/ and $eth!~/^\d+$/) and die print "E-value format error\n";

# special modesa
my $identity = defined $opt{d} ? $opt{d} : "0"; #identity threshold
defined $opt{h} and &usage();

my %id; 
open FH,"$in" or die print "$in not exist\n";
for (<FH>)
{
    my $string = $_;
    chomp;
    /#/ and next;
    my @t=split /\t/,$_;
    
    
    ## Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    $#t!=14 and next;
    my ($qid,$sid,$idn,$aln,$ev,$qcov,$hcov,$taxo) = ($t[0],$t[1],$t[2],$t[3],$t[10],$t[12],$t[13],$t[14]);
    $ev > $eth and next;
    $idn < $identity*100 and next;
    ($idn/100 * $aln) < 25 and next;
    defined $id{$qid}->{EV} and $id{$qid}->{EV} < $ev and next;
    (defined $id{$qid}->{EV} and $id{$qid}->{EV} == $ev ) and (defined $id{$qid}->{IDN} and $id{$qid}->{IDN} > $idn) and next;
    (defined $id{$qid}->{EV} and $id{$qid}->{EV} == $ev ) and (defined $id{$qid}->{IDN} and $id{$qid}->{IDN} == $idn) and (defined    $id{$qid}->{ALN} and $id{$qid}->{ALN} > $aln) and next;
    (defined $id{$qid}->{EV} and $id{$qid}->{EV} == $ev ) and (defined $id{$qid}->{IDN} and $id{$qid}->{IDN} == $idn) and (defined $id{$qid}->{ALN} and $id{$qid}->{ALN} == $aln) and (defined $id{$qid}->{QCOV} and $id{$qid}->{QCOV} > $qcov) and next;
    (defined $id{$qid}->{EV} and $id{$qid}->{EV} == $ev ) and (defined $id{$qid}->{IDN} and $id{$qid}->{IDN} == $idn) and (defined $id{$qid}->{ALN} and $id{$qid}->{ALN} == $aln) and (defined $id{$qid}->{QCOV} and $id{$qid}->{QCOV} == $qcov) and (defined $id{$qid}->{HCOV} and $id{$qid}->{HCOV} > $hcov) and next;

    $id{$qid}={SID=>$sid,EV=>$ev,IDN=>$idn,ALN=>$aln,QCOV=>$qcov,HCOV=>$hcov,BEST=>$string,TAXO=>$taxo};
    
}
close FH;
my ($i,$j)=(0,0); # number of isotigs and genes in blastout
my %iso;
my %gene;
open OT,">$out";
open OUT,">$out\_besthit";
foreach (sort keys %id)
{
    print OT "$_\t$id{$_}->{SID}\t$id{$_}->{EV}\t$id{$_}->{IDN}\n";
    print OUT "$id{$_}->{BEST}";
    not defined $iso{$_} and $i++;
    not defined $gene{$id{$_}->{SID}} and $j++;
    $iso{$_}="o";
    $gene{$id{$_}->{SID}}=$id{$_}->{TAXO};
}
close OT;
close OUT;
print "# of gene-hit isotigs($eth) : $i\n";
print "# of reference gene-hit($eth) : $j\n";

my %taxonomy;
foreach my $k (keys %gene){
    not defined $taxonomy{$gene{$k}} and $taxonomy{$gene{$k}}=0;
    $taxonomy{$gene{$k}}++;
}

open FO, ">$out\_taxonomy";
foreach my $k(sort keys %taxonomy){
    my $v = $taxonomy{$k};
    print FO "$k\t$v\n";
}
close FO;
