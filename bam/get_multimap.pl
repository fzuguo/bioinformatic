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
die "perl $0 <input.order.bam> <output.multialing.bam>\n" unless(@ARGV eq 2);
open FILE,"samtools view -h $ARGV[0]|" or die "failed to open file1";
open OFILE,"|samtools view -Sb - > $ARGV[1]" or die "failed to wirte bam";
my $flag="";;
my $reads1=0;
my $reads2=0;
my @record=();
while(<FILE>){
	if(/^@/){
		print OFILE $_;
	}else{
#		chomp;
		my @as = split("\t",$_);
		my $map = sprintf ("%012b",$as[1]);
		if($flag eq ""){
			$flag = $as[0];
			push(@record,$_);
		}else{
			if($flag eq $as[0]){
				if(substr($map,4,1) eq "1"){
					$reads1++;
				}else{
					$reads2++;
				}
				push(@record,$_)
			}else{
				if($reads1>1 || $reads2 > 1){
					print OFILE join("",@record);
				}
#				print OFILE "\n";
				$flag = $as[0];
				$reads1=0;
				$reads2=0;
				@record = ();
				push(@record,$_)
			}
		}
	}
}
if($reads1>1 || $reads2 > 1){
	print OFILE join("\n",@record);
}
