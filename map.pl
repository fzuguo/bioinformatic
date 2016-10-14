#! /usr/bin/perl

use utf8;
use strict;
use warnings;

#
# this perl script is used to summary ssr locus type:CDS/UTR/intron/intergenic 
# need file: gtf,nr_aligh_table,ssr misa file 
# Contact:jmuguo@foxmail.com
#


die "perl $0 <gtf> <align_table> <ssr_misa> <report> <new_align_file>\n" unless(@ARGV eq 5);
my %hash_misa;

#record misa ssr infomation: unigene:num_type_startlocu
open MISA,$ARGV[2] or die "failed to open misa file";
while(<MISA>){
	chomp;
	my @array = split (/\t/,$_);
	if($hash_misa{$array[0]}){
		$hash_misa{$array[0]} = $hash_misa{$array[0]}."/".$array[1]."_".$array[2]."_".$array[-2];
	}else{
		$hash_misa{$array[0]} = $array[1]."_".$array[2]."_".$array[-2];
	}
}

#generate new hash: (num_types_newloc_chr)/($1)/...
my %hash_gene_start;
open ALIGN,$ARGV[1] or die "failed to open align table";
while (<ALIGN>){
	chomp;
	my @array = split (/\t/,$_);
	if($array[-4] < $array[-3]){
		if($hash_misa{$array[0]} =~ /\//){
			my @locus = split(/\//,$hash_misa{$array[0]});
			foreach my $values(@locus){
				my ($loc) = $values =~ /.*_(.*)/;
				my $newloc = $array[-4] - ($array[-6] - $loc);
				my @deal = split(/_/,$hash_misa{$array[0]});
				if($hash_gene_start{$array[0]}){
					$hash_gene_start{$array[0]} = $hash_gene_start{$array[0]}."/".$deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
				}else{
					$hash_gene_start{$array[0]} = $deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
				}
			}
		}else{
			my ($loc) = $hash_misa{$array[0]} =~ /.*_(.*)/;
			my $newloc = $array[-4] - ($array[-6] - $loc);
			my @deal = split(/_/,$hash_misa{$array[0]});
			$hash_gene_start{$array[0]} = $deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
		}
	}else{
		if($hash_misa{$array[0]} =~ /\//){
			my @locus = split(/\//,$hash_misa{$array[0]});
			foreach my $values(@locus){
				my ($loc) = $values =~ /.*_(.*)/;
				my $newloc = $array[-4] + ($array[-6] - $loc);
				my @deal = split(/_/,$hash_misa{$array[0]});
				if($hash_gene_start{$array[0]}){
					$hash_gene_start{$array[0]} = $hash_gene_start{$array[0]}."/".$deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
				}else{
					$hash_gene_start{$array[0]} = $deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
				}
			}
		}else{
			my ($loc) = $hash_misa{$array[0]} =~ /.*_(.*)/;
			my $newloc = $array[-4] + ($array[-6] - $loc);
			my @deal = split(/_/,$hash_misa{$array[0]});
			$hash_gene_start{$array[0]} = $deal[0]."_".$deal[1]."_".$newloc."_".$array[1];
		}
	}
}
#defind gene zone and type
open GTF,$ARGV[0] or die "failed to open gtf file";
my %hash_gene;
my %hash_trans;
my %hash_cds;
my %hash_3utr;
my %hash_5utr;
while(<GTF>){
	chomp;
	next if ($_ =~ /^#.*/); 
	my @array = split (/\t/,$_);
	if($array[2] eq "gene"){
		if($hash_gene{$array[0]}){
			$hash_gene{$array[0]} = $hash_gene{$array[0]}."/".$array[3]."_".$array[4];
		}else{
			$hash_gene{$array[0]} = $array[3]."_".$array[4];
		}
	}
	if($array[2] eq "transcript"){
		if($hash_trans{$array[0]}){
			$hash_trans{$array[0]} = $hash_trans{$array[0]}."/".$array[3]."_".$array[4];
		}else{
			$hash_trans{$array[0]} = $array[3]."_".$array[4];
		}
	}
	if($array[2] eq "CDS"){
		if($hash_cds{$array[0]}){
			$hash_cds{$array[0]} = $hash_cds{$array[0]}."/".$array[3]."_".$array[4];
		}else{
			$hash_cds{$array[0]} = $array[3]."_".$array[4];
		}
	}
	if($array[2] eq "five_prime_utr"){
		if($hash_5utr{$array[0]}){
			$hash_5utr{$array[0]} = $hash_5utr{$array[0]}."/".$array[3]."_".$array[4];
		}else{
			$hash_5utr{$array[0]} = $array[3]."_".$array[4];
		}
	}
	if($array[2] eq "three_prime_utr"){
		if($hash_3utr{$array[0]}){
			$hash_3utr{$array[0]} = $hash_3utr{$array[0]}."/".$array[3]."_".$array[4];
		}else{
			$hash_3utr{$array[0]} = $array[3]."_".$array[4];
		}
	}
}
my $count_3utr = 0;
my $count_5utr = 0;
my $count_intergenic = 0;
my $count_cds = 0;
my $count_intron = 0;
my $total_align = 0;
open NEW_ALIGN,">./$ARGV[4]" or die "faild to write new align file";
foreach my $keys(keys %hash_gene_start){
	my @array_gene = split(/\//,$hash_gene_start{$keys});
	foreach my $values(@array_gene){
		my ($locus,$chr) = $values =~ /.*?_.*?_(.*?)_(.*)/;
		if(&judge($hash_gene{$chr},$locus) eq "in"){

			if(&judge($hash_3utr{$chr},$locus) eq "in"){
				$count_3utr ++;
				print NEW_ALIGN $keys."\t".$values."\t3utr\n";
			}elsif(&judge($hash_5utr{$chr},$locus) eq "in"){
				$count_5utr ++;
				print NEW_ALIGN $keys."\t".$values."\t5utr\n";
			}elsif(&judge($hash_cds{$chr},$locus) eq "in"){
				print NEW_ALIGN $keys."\t".$values."\tcds\n";
				$count_cds ++;
			}else{
				print NEW_ALIGN $keys."\t".$values."\tintron\n";
				$count_intron ++;
			}
		}else{
			print NEW_ALIGN $keys."\t".$values."\tintergenic\n";
			$count_intergenic ++ ;
		}
		$total_align ++;
	}
}
#report
my $total = `cat $ARGV[2]|sed 1d|cut -f 1|sort|uniq|wc -l`;
my $unmap = $total - $total_align;
open REPORT,">./$ARGV[3]" or die "failed to write out report file";
print REPORT "3UTR_count\t$count_3utr\n5UTR_count\t$count_5utr\nintergenic_count\t$count_intergenic\ncds_count\t$count_cds\nintron_count\t$count_intron\nunmap\t$unmap\n";
sub judge{
	my $zone = shift;
	my $loc = shift;
	my @zones = split(/\//,$zone);
	foreach my $value(@zones){
		my ($start,$end) = $value =~ /(.*?)_(.*)/;
		if($loc >= $start && $loc <= $end){
			return "in";
		}
	}
	return "out";
}
