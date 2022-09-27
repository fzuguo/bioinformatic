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
die "perl $0 <IN.bam> <out.bam>\n" unless(@ARGV eq 2);
open FILE,"samtools view -h $ARGV[0]|" or die "failed to open file1";
open OFILE,"|samtools view -Sb - > $ARGV[1]" or die "failed to open file1";
my $flag="";
my @record=();
while(<FILE>){
	if(/^@/){
		print OFILE $_;
		next;
	}
	my $line1 = $_;
	my @as = split("\t",$_);
	if($flag eq ""){
		$flag = $as[0];
	}
	if($flag ne $as[0]){
		my $map = 0;
		foreach my $r(@record){
			my @as2 = split("\t",$r);
			if($as2[2] ne "*"){
				$map = 1;
			}
		}
		if($map == 1){
			foreach my $r(@record){
				print OFILE $r;
			}
		}
		@record = ();
		push(@record,$_);
		$flag = $as[0];
	}else{
		push(@record,$_);
	}
}
my $map = 0;
foreach my $r(@record){
	my @as2 = split("\t",$r);
	if($as2[2] ne "*"){
		$map = 1;
	}
}
if($map == 1){
	foreach my $r(@record){
		print OFILE $r;
	}
}
