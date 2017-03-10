#!/usr/bin/perl -w

use strict;

die "

    Example: $0 blastx_trnity_arthopod.txt Trinity.fasta 1e-10 0.5 30 Filtered_Trinity
    Command: $0 [BLAST result] [FASTA sequences of de novo assembled transcripts] [e-value cutoff] [identity cutoff] [cutoff of isoform number for a gene] [output prefix]

" if !@ARGV;

$0 =~ /(.+?)\/run.pl/;
my $lib_path = "$1/lib";

my ($blast, $fasta, $ev, $idn, $isoNum, $out_prefix)=@ARGV;

print "\n";
my $cmd1 = "perl $lib_path/get_besthit_idn_25aa.pl -i $blast -o $blast.best_$idn -e $ev -d $idn &> $blast.best_$idn.log";
print "[CMD] Selete best hit reference ID for each transcript: $cmd1 ...\n";
`$cmd1`;

my $cmd2 = "perl $lib_path/extract_seq.pl $blast.best_$idn\_besthit $fasta $out_prefix $isoNum";
print "[CMD] Extract the transcripts with blast hit homologs from the reference sequences: $cmd2 ...\n";
`$cmd2`;

print "\n";
