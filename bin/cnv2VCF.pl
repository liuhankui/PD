#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $usage = "\tperl cnv2vcf.pl -prefix prefix -fai GRCh38.fa.fai\n";

my $prefix;
my $fai;
my %h;
die("not enough argyments. $usage\n") unless ( @ARGV );

GetOptions( 'p|prefix:s' => \$prefix,'f|fai:s' => \$fai );
open(FI,"$fai");
while(<FI>){
        chomp;
        my @t=split;
        $h{$t[0]}=$t[1];
}

print '##fileformat=VCFv4.1',"\n";
print '##fileDate='.`date '+%Y%m%d'`;
print '##source=CNVpytor',"\n";
print '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',"\n";
print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',"\n";
print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',"\n";
print '##INFO=<ID=RD,Number=1,Type=Float,Description="Normalized RD">',"\n";
print '##INFO=<ID=P1,Number=1,Type=Float,Description="e-value (p-value multiplied by genome size divided by bin size) calculated using t-test statistics between RD statistics in the region and global">',"\n";
print '##INFO=<ID=P2,Number=1,Type=Float,Description="e-value (p-value multiplied by genome size divided by bin size) from the probability of RD values within the region to be in the tails of a gaussian distribution of binned RD">',"\n";
print '##INFO=<ID=P3,Number=1,Type=Float,Description="same as e-val1 but for the middle of CNV">',"\n";
print '##INFO=<ID=P4,Number=1,Type=Float,Description="same as e-val2 but for the middle of CNV">',"\n";
print '##INFO=<ID=Q0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality in call region">',"\n";
print '##INFO=<ID=N,Number=1,Type=Integer,Description="Fraction of reference genome gaps in call region">',"\n";
print '##INFO=<ID=D,Number=1,Type=Integer,Description="distance from closest large (>100bp) gap in reference genome">',"\n";
print '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">',"\n";
print '##ALT=<ID=DEL,Description="Deletion">',"\n";
print '##ALT=<ID=DUP,Description="Duplication">',"\n";
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',"\n";
print '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',"\n";
print '##FORMAT=<ID=RD,Number=1,Type=Float,Description="Normalized RD">',"\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$prefix\n";
my $count = 0;
while (my $line = <>) {
    my ($type,$coor,$len,$rd,$p1,$p2,$p3,$p4,$q0,$pe,$dis) = split(/\s+/,$line);
    my ($chrom,$start_end) = split(/:/,$coor);
    my ($start,$end) = split(/-/,$start_end);
    $end=$h{$chrom} if ($end>$h{$chrom});
    my $isDel = ($type eq "deletion");
    my $isDup = ($type eq "duplication");
    if ($isDup) {
    } elsif ($isDel) {
    } else {
        next;
    }
    print $chrom,"\t",$start,"\t",$coor,"\t.\t";
    if    ($isDel) { print "<DEL>"; }
    elsif ($isDup) { print "<DUP>"; }
    print "\t.\tPASS\t";
    my $INFO = "END=".$end;
    if    ($isDel) {
        $INFO .= ";SVTYPE=DEL";
        $INFO .= ";SVLEN=-".int($len);
    } elsif ($isDup) {
        $INFO .= ";SVTYPE=DUP";
        $INFO .= ";SVLEN=".int($len);
    }
    if (defined($rd) && ($rd ne "")) { $INFO .= ";RD=".$rd; }
    if (defined($p1) && ($p1 ne "")) { $INFO .= ";P1=".$p1; }
    if (defined($p2) && ($p2 ne "")) { $INFO .= ";P2=".$p2; }
    if (defined($p3) && ($p3 ne "")) { $INFO .= ";P3=".$p3; }
    if (defined($p4) && ($p4 ne "")) { $INFO .= ";P4=".$p4; }
    if (defined($q0) && ($q0 ne "")) { $INFO .= ";Q0=".$q0; }
    if (defined($pe) && ($pe ne "")) { $INFO .= ";N=".$pe; }
    if (defined($dis) && ($dis ne "")) { $INFO .= ";D=".$dis; }
    print $INFO;

    my $GT="GT";

    if(defined($rd) && ($rd ne "")) {
        $GT.=":CN";
        if(defined($rd) && ($rd ne "")) {
            $GT.=":RD";
        }
        $GT.="\t";

        if ($isDel && $rd < 0.25) {
            $GT .= "1/1:0";
        } elsif ($isDel && $rd >= 0.25) {
            $GT .= "0/1:1";
        } elsif ($isDup && $rd <= 1.75) {
            $GT .= "0/1:2";
        } elsif ($isDup && $rd > 1.75 && $rd <= 2.25) {
            $GT .= "1/1:2";
        } elsif ($isDup && $rd > 2.25) {
            my $cn = sprintf("%.0f",$rd);
            $GT.="./1:$cn"; # w/o other data, we can't really say if this is
                            # a hom dup, or het dup with higher copy number.
        } else {
            $GT = "GT\t./.";
        }

        if(defined($rd) && ($rd ne "")) {
            $GT.=":$rd";
        }
    } else {
        $GT.="\t./.";
    }
    print "\t$GT\n";
}

exit;
