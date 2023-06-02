# PD WGS (structure variants, copy number variants,short tandem repeats) GWAS

## fq2cram

```
fastp -q 20 -u 10 -n 5 --in1 fq1.gz --in2 fq2.gz --out1 clean.fq1.gz --out2 clean.fq2.gz
bwa mem -t 4 -R '@RG\tID:id\tPL:illumina\tPU:id\tLB:sample\tSM:id\tCN:BGI' GRCh38.fa clean.fq1.gz clean.fq2.gz|samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam|samtools view -Sb -|samtools sort - -O CRAM -o sort.cram --reference GRCh38.fa
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam|sambamba view -S -f bam /dev/stdin|samtools sort - -O CRAM -o discordants.cram --reference GRCh38.fa
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam|sambamba view -S -f bam -l 0 /dev/stdin|samtools sort - -O CRAM -o splitters.cram --reference GRCh38.fa
```

## cram2depth

```
mosdepth qc sort.cram -f GRCh38.fa --fast-mode --no-per-base --by 1000000 --thresholds 1,2,4,10
```

## cran2variants

```
# SV
lumpyexpress -P -T ./tmp -R GRCh38.fa -B sort.cram -S splitters.bam -D discordants.bam -x exclude.bed -o sv.vcf && svtyper -B sort.cram -T GRCh38.fa -i sv.vcf |gzip -f > sv.gt.vcf.gz
# CNV
cnvpytor -root root.pytor -rd sort.cram -T GRCh38.fa && cnvpytor -root root.pytor -his 100 && cnvpytor -root root.pytor -partition 100 && cnvpytor -root root.pytor -call 100 | perl cnv2vcf.pl -prefix id -fai GRCh38.fa.fai|bgzip -f > cnv.vcf.gz
# STR
ExpansionHunter --reads sort.cram --reference GRCh38.fa --variant-catalog variant_catalog.json --output-prefix ./kSTR && vawk --header '{$8="SVTYPE=DUP;"$8;print}'|bgzip -f > kSTR.vcf.gz && tabix -f -s 1 -b 2 -e 2 kSTR.vcf.gz
ExpansionHunterDenovo profile --reads sort.cram --reference GRCh38.fa --output-prefix dSTR --min-anchor-mapq 50 --max-irr-mapq 40
# SNV/INDEL
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar HaplotypeCaller --QUIET true -R GRCh38.fa -I sort.cram -O gatk.gvcf.gz -ERC GVCF -A ClippingRankSumTest -A LikelihoodRankSumTest -A MappingQualityZero && java -Xmx3g -jar gatk-package-4.1.2.0-local.jar GenotypeGVCFs --QUIET true -R GRCh38.fa --variant gatk.gvcf.gz -O gatk.vcf.gz
```

## combined variants

### SV

```
svtools lsort -f sv.list -t ./tmp -b 100 > sorted.sv.vcf
svtools lmerge -i sorted.sv.vcf -f 20 -t ./tmp --sum > merged.sv.vcf

cat sv.list | while read sv
dir=`dirname $sv`
vawk '{if(I$END-$2>100)print $1":"$2"-"I$END}'|cnvpytor -root $dir/root.pytor -genotype 100 > $sv.genotype
vawk '{if(I$END-$2>100)print}' merged.sv.vcf|svtyper -l $sv.lib.json -B $dir/sort.cram -T GRCh38.fa -i /dev/stdin|perl cn2GT.pl > $sv.gt.vcf
python del_pe_resolution.py $sv.lib.json | awk 'NR==2{print}'
done > sv.lib.size

cat sv.list|awk '{print $0".gt.vcf"}'|svtools vcfpaste -f /dev/stdin -q|svtools afreq|svtools vcftobedpe|svtools bedpesort|svtools prune -s -d 100 -e "AF"|svtools bedpetovcf|svtools classify -g sex.txt -a repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.GRCh38.sorted.bed.gz -m large_sample|python geno_refine_12.py -i - -g sex.txt -d refine.dfile.txt|python filter_del.py -i - -t sv.lib.size -s 0.1|resvtyper.py -|vawk --header '{if(I$MSQ!="." && !((I$SVTYPE=="BND" && I$MSQ<250) || (I$SVTYPE=="INV" && I$MSQ<150))){print}}'|gzip -f > sv.gt.filter.vcf.gz
zcat sv.gt.filter.vcf.gz|vawk '{if($8!~/SECONDARY/)print $1"_"$2"_"I$END"_"I$SVTYPE,S$*$GT}'|gzip -f > sv.genotype.gz
```

sv.list file format：
```
path1/sv.vcf
path2/sv.vcf
......
path9/sv.vcf
```

### CNV

```
cat cnv.list|while read cnv
do
id=`dirname $cnv|awk -F '/' '{print $NF}'`
zcat $cnv|vawk -v id=$id '{if(I$P1<0.0001 && I$Q0<0.5 && I$N<0.2 && I$D>10000)print id,$1,$2,I$END,I$SVTYPE}'
done > cnv.txt

awk '$5=="DUP"{print $2"\t"$3"\t"$4"\t"$1"_"$4-$3}' cnv.txt|sort -k1,1 -k2,2n -k3,3n > DUP.sort.txt
awk '$5=="DEL"{print $2"\t"$3"\t"$4"\t"$1"_"$4-$3}' cnv.txt|sort -k1,1 -k2,2n -k3,3n > DEL.sort.txt
bedtools merge -i DUP.sort.txt -c 1,4 -o count,collapse > cnv.overlap.txt
bedtools merge -i DEL.sort.txt -c 1,4 -o count,collapse >> cnv.overlap.txt
perl cnv.combine.pl cnv.list cnv.overlap.txt|gzip -f > cnv.size.gz
```

cnv.list file format
```
path1/cnv.vcf.gz
path2/cnv.vcf.gz
......
path9/cnv.vcf.gz
```

### STR

```
bcftools merge -l kSTR.list -m all -O z -o kSTR.vcf.gz
zcat kSTR.vcf.gz|perl str.genotype.pl|gzip -f > kSTR.genotype.gz
ExpansionHunterDenovo merge --reference GRCh38.fa --manifest dSTR.list --output-prefix dSTR
```

kSTR.list file format
```
path1/kSTR.vcf.gz
path2/kSTR.vcf.gz
......
path9/kSTR.vcf.gz
```

dSTR.list file format
```
id1 case path1/dSTR_profile.json
id2 control path2/dSTR_profile.json
......
id9 case path9/dSTR_profile.json
```

## GWAS

### SV & CNV

```
fdf<-read.table('plink.fam')
y<-fdf[,6]
sex<-fdf[,5]
age<-fdf[,7]
pc1<-fdf[,8]
pc2<-fdf[,9]
pc3<-fdf[,10]

#sv
gdf<-read.table(gzfile('sv.genotype.gz'))
gdf[gdf=='./.']<-NA
gdf[,-1]<-as.numeric(factor(gdf[,-1],levels=c('0/0','0/1','1/1'),order=T))
pdf<-as.data.frame(t(apply(gdf[,-1],1,function(x){coef<-summary(lm(y ~ sex + age + pc1 + pc2 + pc3 + x))$coefficients;return(coef[7,-3])})))
pdf<-cbind(gdf[,1],pdf)
write.table(pdf,file='sv.gwas.txt',row.names=F,col.names=F,quote=F)

#cnv
gdf<-read.table(gzfile('cnv.size.gz'))
pdf<-as.data.frame(t(apply(gdf[,-1],1,function(x){coef<-summary(lm(y ~ sex + age + pc1 + pc2 + pc3 + x))$coefficients;return(coef[7,-3])})))
pdf<-cbind(gdf[,1],pdf)
write.table(pdf,file='cnv.gwas.txt',row.names=F,col.names=F,quote=F)
```

### STR

kSTR
```
gdf<-read.table(gzfile('kSTR.genotype.gz'))
gdf[,2]<-apply(gdf[,-1],1,function(x){wilcox.test(x~y,alternative='less')$p.value})
write.table(gdf[,c(1,2)],file='kSTR.gwas.txt',row.names=F,col.names=F,quote=F)
```

dSTR
```
casecontrol.py locus --manifest dSTR.list --multisample-profile dSTR.multisample_profile.json --output dSTR.GWAS.tsv
```

