#!/usr/bin/perl
use strict;

while(<>){
        chomp;
        next if(/\#/);
        my @t=(split/\t/,$_,10);
        next if($t[4] eq '.');
        my @g=(split/\t/,$t[-1]);

        $t[7]=~/END=(\d+);REF=(\d+);RL=(\d+);RU=([A|T|C|G|N|R]+);/;
        my $e=$1;
        my $r=$2;
        my $u=$4;

        $t[4]=~s/STR|\>|\<//g;
        my @a=(split/,/,$t[4]);
        unshift @a,$r;
        print "$t[0]_$t[1]_${e}_${u}";

        for(my $i=0;$i<=$#g;$i++){
                if($g[$i]=~/\.\/\./){
                        print " NA";
                }else{
                        my @str=(split/\//,(split/\:/,$g[$i])[0]);
                        my $add=$a[$str[0]]+$a[$str[1]];
                        print " $add";
                }
        }
        print "\n";
}
