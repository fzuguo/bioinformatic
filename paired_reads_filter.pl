#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use Fzuguo;
#use Getopt::Long;
#use Bio::SeqIO;
#use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;

die "perl $0 <reads1.fq.gz> <reads2/fq/gz>\n" unless(@ARGV eq 2);
open FILE,"gzip -dc $ARGV[0]|" or die "failed to open file1";
open FILE2,"gzip -dc $ARGV[1]|" or die "failed to open file2";
open OUT,"|gzip > Reads1.fq.gz" or die "failed to write reads1";
open OUT2,"|gzip > Reads2.fq.gz" or die "failed to write reads2";
#on:origin_num; dn:dealed_num; dupn,duplicated_num; r1:reads1; r2:reads2
my $r1on=0;
my $r1dn=0;
my $r2on=0;
my $r2dn=0;
my $r1dupn=0;
my $r2dupn=0;
my %hash_reads1;
while(<FILE>){
	chomp;
	my $line1 = $_;
	my $line2 = <FILE>;
	my $line3 = <FILE>;
	my $line4 = <FILE>;
	$line1 =~ s/\/1$//;
#	$hash_reads1{$line1} = $line1."\/1\n".$line2.$line3.$line4;
	$hash_reads1{$line1} = 1;
	$r1on++;
}
my %exist_reads;
while(<FILE2>){
	chomp;
	my $line1 = $_;
	my $line2 = <FILE2>;
	my $line3 = <FILE2>;
	my $line4 = <FILE2>;
	$line1 =~ s/\/2$//;
	if(defined($hash_reads1{$line1})){
		if($hash_reads1{$line1} == 1){
			print OUT2 $line1."\/2\n".$line2.$line3.$line4;
			$exist_reads{$line1} = 1;
			$hash_reads1{$line1} = 1;
			$r2dn++
		}else{
			$r2dupn++
		}
	}
	$r2on++;
}
open FILE3,"gzip -dc $ARGV[0]|" or die "failed to open file1";
while(<FILE3>){
	chomp;
	my $line1 = $_;
	my $line2 = <FILE3>;
	my $line3 = <FILE3>;
	my $line4 = <FILE3>;
	$line1 =~ s/\/1$//;
	if(defined($exist_reads{$line1})){
		if($exist_reads{$line1} == 1){
			print OUT $line1."\/1\n".$line2.$line3.$line4;
			$exist_reads{$line1}=0;
			$r1dn++;
		}else{
			$r1dupn++;
		}
	}
	
}

open REPORT,">./report.txt";
print REPORT "origin_reads1_num: $r1on\tdealed_reads1_num:$r1dn\tfilter_paired_reads1_dup_num:$r1dupn\n";
print REPORT "origin_reads2_num: $r2on\tdealed_reads2_num:$r2dn\tfilter_paired_reads2_dup_num:$r2dupn\n";

