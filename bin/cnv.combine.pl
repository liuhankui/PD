#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $usage = "\tperl cnv.combine.pl -ped ped.file\n";
my $ped;

die("not enough argyments. $usage\n") unless ( @ARGV );

GetOptions( 'ped|ped:s' => \$ped );

my %h;
my $i=0;
open(ID,"$ped");
while(<ID>){
        chomp;
        my @t=split;
        $h{$t[1]}=$i;
        $i++;
}

while(<>){
        chomp;
        my @t=split;
        next if($t[3]<2);
        my @g=(split/,/,$t[-2]);
        my %c;
        foreach my $k(@g){
                my @info=(split/_/,$k);
                $c{$info[0]}+=$info[1];
        }
        print "$t[0]_$t[1]_$t[2]_$t[-1]";
        foreach my $k(sort {$h{$a} <=> $h{$b}} keys %h){
                my $w1=$c{$k} || 0;
                print " $w1"
        }
        print "\n";
}
