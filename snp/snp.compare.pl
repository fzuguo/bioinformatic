#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use lib "/home/guojianhui/lib";
use Fzuguo;
use File::Basename;
my $pwd = dirname(__FILE__);
#use Getopt::Long;
#use Bio::SeqIO;
#use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;
#Author: jhguo@genedenovo.com
die "perl $0 <snp.annot.xls> <compare.txt>\n" unless(@ARGV eq 2);
open FILE,$ARGV[0] or die "failed to open file1";
my $head = <FILE>;
chomp($head);
my @heads = split("\t",$head);
my $num=0;
my %id2name;
my %ts = (
    'AG' => 1,
    'GA' => 1,
    'TC' => 1,
    'CT' => 1,
);
my %tv = (
    'AC' => 1,
    'CA' => 1,
    'AT' => 1,
    'TA' => 1,
    'GC' => 1,
    'CG' => 1,
    'GT' => 1,
    'TG' => 1,
);
my %compairs;
my %fileh;
open CMP,$ARGV[1] or die $!;
while(<CMP>){
	chomp;
	my @as = split("\t",$_);
	$compairs{$as[0]} = [$as[1],$as[2]];
	if(exists $fileh{$as[0]}){
	}else{
		open my $tmp,">$as[0].xls" or die $!;
		$fileh{$as[0]} = $tmp;
	}
}


my %id2col;
for(my $i=8;$i<=scalar(@heads)-8;$i++){
	if($heads[$i]=~ /(.*?)_type/){
		$num++;
		$id2name{$i} = $1;
		$id2col{$1} = $i;
		print $i."\t".$1."\n";
	}
}
foreach my $compair(keys %compairs){
	my $tmp = $fileh{$compair};
#	open $tmp,">$compair.xls" or die $!;
	print $tmp join("\t",@heads[0..7]);
	my @g1 = split(",",$compairs{$compair}->[0]);
	my @g2 = split(",",$compairs{$compair}->[1]);
	foreach my $s((@g1,@g2)){
		print $tmp  "\t".$heads[$id2col{$s}];
	}
	foreach my $s((@g1,@g2)){
		print $tmp "\t".$heads[$id2col{$s}+$num];
	}
	print $tmp "\t".join("\t",@heads[-4..-1])."\ttype\n"; 
}
while(<FILE>){
	chomp;
	my @as = split("\t",$_);
	foreach my $compair(keys %compairs){
		my $tmp = $fileh{$compair};
		my @g1 = split(",",$compairs{$compair}->[0]);
		my @g2 = split(",",$compairs{$compair}->[1]);
		my $type_g1 = "0";
		my $type_g2 = "0"; 
		foreach my $s(@g1){
			if($as[$id2col{$s}]=~ /[1-9]/){
				$type_g1 = "1";
				last;
			}
		}
		foreach my $s(@g2){
			if($as[$id2col{$s}]=~ /[1-9]/){
				$type_g2 = "1";
				last;
			}
		}
		if($type_g1 ne $type_g2){
			print $tmp join("\t",@as[0..7]);
			foreach my $s((@g1,@g2)){
				print $tmp  "\t".$as[$id2col{$s}];
			}
			foreach my $s((@g1,@g2)){
				print $tmp  "\t".$as[$id2col{$s}+$num];
			}
			print $tmp "\t".join("\t",@as[-4..-1]);
			my $type ="";
			my @alts = split(",",$as[4]);
			foreach my $alt(@alts){
				if($tv{$as[3].$alt}){
					$type = $type.","."transversion";
				}elsif($ts{$as[3].$alt}){
					$type = $type.","."transition";
				}
			}
			$type = "indel" if ($type eq "");
			print $tmp "\t$type\n";
		}
	}
}
