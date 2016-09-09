#! /usr/bin/perl
#Author:jianhuiguo@126.com,this Program is used for check raw reads align rate with a refer fa file;need bowtie2 suggest
use utf8;
use strict;
use warnings;

die "perl $0 <raw_fq1> <raw_fq2> <genome.fa> <prefix> <reads_num>\n use top (suggest:100000) reads to check align rate\n perl $0 raw_xx_1.fq.gz raw_xx_2.fq.gz genome.fa xxx 10000" if @ARGV ne 5;

my($name1) = "/$ARGV[0]" =~ /.*\/(.*)\.gz/;
my($name2) = "/$ARGV[1]" =~ /.*\/(.*)\.gz/;
my($name3) = "/$ARGV[0]" =~ /.*\/(.*?)_/;
print $name1."\t".$name2."\n";
my $num = $ARGV[4]*4;
unless(-s $name1){
`zcat $ARGV[0]|head -n $num > ./$name1`;
`gzip ./$name1`;
}
unless(-s $name2){
`zcat $ARGV[1]|head -n $num > ./$name2`;
`gzip ./$name2`;
}
if(-s "$ARGV[3].1.bt2"){
	print "bowtie2 index exist";
}else{
	`bowtie2-build $ARGV[2] ./$ARGV[3]`;
}
`bowtie2 --local -S /dev/null -p 4 --mm  -x ./$ARGV[3] -1 ./$name1.gz -2 ./$name2.gz --un-conc-gz ./$name3.%.fq.gz 2> ./$name3.align.log`
