#! /usr/bin/perl

use utf8;
use strict;
use warnings;
#use Getopt::Long;
#use Bio::SeqIO;
#use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;
#Author: jianhuiguo@126.com,To know this Program,you shoud know sam file format[https://samtools.github.io/hts-specs/SAMv1.pdf]
die "perl $0 <align.sam> <out_prefix>\nused to get failed align paired reads from samfile" unless(@ARGV eq 2);
open SAM,$ARGV[0] or die "failed to open SAM file";
open READS1,"| gzip > $ARGV[1]_1.fq.gz" or die "failed to write output file1";
open READS2,"| gzip > $ARGV[1]_2.fq.gz" or die "failed to write output file1";
my %flag;
my %reads;
my $lastread;
while(<SAM>){
	next if (/^[@#]/);
	chomp;
	my @array = split (/\t/,$_);
	if($array[2] eq "*"){
		#check if paired reads is all failed aligned to ref fa
		$flag{$array[0]}++;
		if($flag{$array[0]} eq "1"){
			$reads{$array[0]."/1"} = "@".$array[0]."/1\n".$array[9]."\n+\n$array[10]\n";
		}elsif($flag{$array[0]} eq "2"){
			$reads{$array[0]."/2"} = "@".$array[0]."/2\n".$array[9]."\n+\n$array[10]\n";
			print READS1 $reads{$array[0]."/1"};
			print READS2 $reads{$array[0]."/2"};
			#if flag eq 2,means paired reads is all failed aligned to ref fa,then write out;and release related hash
			delete $reads{$array[0]."/1"};
			delete $reads{$array[0]."/2"};
		}
		
	}
	#delete related hash 
	if($lastread ne "" && $lastread ne $array[0]){
		delete $reads{$lastread."/1"};
		delete $reads{$lastread."/2"};
	}
	#rememer last read;
	$lastread = $array[0];
}

