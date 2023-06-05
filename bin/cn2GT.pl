#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $usage = "\tperl cn2GT.pl -cn cn.txt -vcf vcf.file\n";
my $cn;
my $vcf;

die("not enough argyments. $usage\n") unless ( @ARGV );

GetOptions( 'cn|cn:s' => \$cn,'vcf|vcf:s' => \$vcf );

my %h;
open(CN,"$cn");
while(<CN>){
        chomp;
        my @t=split;
        $h{$t[0]}=$t[1];
}

open(VCF,"$vcf");
while(<VCF>){
        chomp;
        if(/^#/){
                print "$_\n";
                next;
        }
        my @t=split;
        $t[8]="$t[8]:CN";
        $t[7]=~/END=(\d+);/;
        my $e=$1;
        my $m="$t[0]:$t[1]-$e";
        next if(!exists($h{$m}));
        $t[9]="$t[9]:$h{$m}";
        print join("\t",@t),"\n";
}
