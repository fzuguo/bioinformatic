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
#Author: jmuguo@foxmail.com
die "perl $0 <NCBI_GFF> <output_GTF>\n" unless(@ARGV eq 2);
#this program use id as ref index
open GFF,$ARGV[0] or die "failed to open gtf file";
my (%gene_type,%gene_id,%gene_name);
my (%mrna_type,%mrna_geneid,%mrna_name);
#build gene info include:type name geneid
while(<GFF>){
	next if /^#/;
	chomp;
	my $line = $_;
	$line =~ s/;/ ;/g;
	$line =~ s/,/ ,/g;
	my @array = split(/\t/,$line);
	if($array[2] =~ /^gene$/i){
		my ($id) = $line =~ /ID=(\S+)\s/i;
		if($line =~ /gene_biotype=(\S+)\s/){
			$gene_type{$id} = $1;
		}else{
			$gene_type{$id} = "unKnown";
		}
		if($line =~ /GeneID:(\S+)\s/i){
			$gene_id{$id} = "GeneID_$1";
			$line =~ /Name=(\S+)\s/i;
			$gene_name{$id} = $1;
		}elsif($line =~ /Name=(\S+)\s/i){
			$gene_id{$id} = $1;
			$gene_name{$id} = $gene_id{$id};
		}
	}
}
seek GFF,0,0;
#build mrna info,include mrna_name mrna_type mrna_id
while(<GFF>){
	next if /^#/;
	chomp;
	my $line = $_;
	$line =~ s/;/ ;/g;
	$line =~ s/,/ ,/g;
	my @array = split(/\t/,$line);
	my ($parent) = $line =~ /Parent=(\S+)\s/i;
	if($array[2] =~ /^mrna$/i || $array[2] =~ /^rRNA$/i || $array[2] =~ /^tRNA$/i){
		my ($id) = $line =~ /ID=(\S+)\s/i;
		my ($parent) = $line =~ /Parent=(\S+)\s/i; 
		#mrna biotype same as it's parent gene biotype
		$mrna_type{$id} = $gene_type{$parent};
		#get mrna's parent's geneid
		$mrna_geneid{$id} = $gene_id{$parent};
		#mrna_name
		if($line =~ /Genbank:(\S+)\s/i){
			$mrna_name{$id} = $1;
		}elsif($line =~ /Name=(\S+)\s/i){
			#if didn't exist Genbank name,use Name replaced
			$mrna_name{$id} = $1;
		}else{
			$mrna_name{$id} = $gene_id{$parent};
		}
	}
}
seek GFF,0,0;
#sometimes  genes or mrna lost Genename info then get the message from 
while(<GFF>){
	next if /^#/;
	chomp;
	my $line = $_;
	$line =~ s/;/ ;/g;
	$line =~ s/,/ ,/g;
	my @array = split(/\t/,$line);
	if($array[2] =~ /^cds$/i || $array[2] =~ /^exon$/i){
		my ($parent) = $line =~ /Parent=(\S+)\s/i;
		if($parent =~ /gene/i){
			if(!defined $mrna_name{$parent}){
				if($line =~ /Genbank:(\S+)\s/i){
					$mrna_name{$parent} = $1;
				}elsif($line =~ /Name:(\S+)\s/i){
					$mrna_name{$parent} = $1;
				}else{
					$mrna_name{$parent} = $gene_name{$parent};
				}
			}
			if(!defined $mrna_geneid{$parent}){
				$mrna_geneid{$parent} = $gene_id{$parent};
			}
			if(!defined $mrna_type{$parent}){
				$mrna_type{$parent} = $gene_type{$parent};
			}
		}elsif($parent =~ /rna/i){
			if(!defined $mrna_name{$parent}){
				if($line =~ /Genbank:(\S+)\s/i){
					$mrna_name{$parent} = $1;
				}elsif($line =~ /Name:(\S+)\s/i){
					$mrna_name{$parent} = $1;
				}
			}
			if(!defined $mrna_geneid{$parent}){
				$mrna_geneid{$parent} = $gene_id{$parent};
			}
			if(!defined $mrna_type{$parent}){
				$mrna_type{$parent} = $gene_type{$parent};
			}
		}
	}
}
#report
open OUTGTF,">./$ARGV[1]" or die "failed to write output file";
seek GFF,0,0;
while(<GFF>){
	next if /^#/;
	chomp;
#print "deal:"$_.;
	my $line = $_;
	$line =~ s/;/ ;/g;
	$line =~ s/,/ ,/g;
	my @array = split(/\t/,$line);
	my $tail = pop @array;
	my $text = join("\t",@array);
	if($array[2] =~ /^cds$/i || $array[2] =~ /^exon$/i){
		my ($parent) = $tail =~ /Parent=(\S+)\s/i;
			print OUTGTF $text."\ttranscript_id \"$mrna_name{$parent}\"; gene_id \"$mrna_geneid{$parent}\"; gene_name \"$mrna_geneid{$parent}\"; transcript_biotype \"$mrna_type{$parent}\";\n"
	}
}

