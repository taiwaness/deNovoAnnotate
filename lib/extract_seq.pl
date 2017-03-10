#!/usr/bin/perl -w

use strict;

die "

    Example: $0 blastx_AA_trnity_new.txt.best_0.5_idn_25aa_besthit sequences/AA.fa sequences/AA 60
    Command: $0 blast_best_hit fasta output number_of_isotigs

" if !@ARGV;

my ($blast, $fasta, $out, $cutoff) = @ARGV;

my %fa = ();
my %all_isoform_num=();
my ($gid, $tid) = ((),());
open FI, "$fasta";
while (<FI>){
    chomp $_;
    if ($_ =~ /^>(.+?)_i\d+\s+/){
        $gid = $1;
        if (defined $fa{$gid}){
            $fa{$gid} .= "\n$_";
            $all_isoform_num{$gid} += 1;
        }else{
            $fa{$gid} = $_;
            $all_isoform_num{$gid} = 1;
        }
    }else{
        if (defined $fa{$gid}){
            $fa{$gid} .= "\n$_";
        }else{
            print "[ERROR1]\tNo ID for the sequence, $_\n";
            exit;
        }
    }
}
close FI;

my %hspcount = ();
open FG, "> $out\_isogroup.fa";
open FT, "> $out\_isotig.fa";
open FLT, "> $out\_longest_isotig.fa";
open FI, "$blast";
while (<FI>){
    chomp $_;
    my @t = split("\t", $_);
    $t[0] =~ /(.+?)_i\d+$/;
    my $gid = $1;
    my $tid = $t[0];
    if ($all_isoform_num{$gid} <= $cutoff){
        my %tmpfa = ();
        my ($maxlen, $maxid) = (0,());
        my $id = ();
        my @tfa = split(/\n/, $fa{$gid});
        foreach my $i (0 .. $#tfa){
            if ($tfa[$i] =~ /^>(\S+)/){
                $id = $1;
                $tfa[$i] =~ /len=(\d+)/;
                if ($1 > $maxlen){ $maxlen = $1; $maxid = $id; }
                if (defined $tmpfa{$id}){}else{ $tmpfa{$id}=$tfa[$i]; }
            }else{
                if (defined $tmpfa{$id}){
                    $tmpfa{$id} .= "\n$tfa[$i]";
                }else{
                    print "[ERROR2]\tNo ID for the sequence, $_\n";
                    exit;
                }
            }
        }

        if (defined $tmpfa{$tid}){
            print FT "$tmpfa{$tid}\n";
        }

        if (defined $hspcount{$gid}){ 
            $hspcount{$gid} += 1;
        }else{
            $hspcount{$gid} = 1;
        }

        if ($hspcount{$gid} > 1){ next; }
        if (defined $tmpfa{$maxid}){
            print FLT "$tmpfa{$maxid}\n";
        }
        print FG "$fa{$gid}\n";
    }
}
close FI;
close FG;
close FT;
close FLT;
