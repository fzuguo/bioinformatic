#!/usr/bin/perl
use warnings;
#use strict;
#this script used for summary gene methylation rate,including upstream_2K,genebody,downstream_2K

#used to generate ME ratio,UPstream_2k,GeneBody,Downstream_2k

die "perl $0 <DIR:count_dir> <FILE:chr_start_end_gene>" if @ARGV ne 2;

#open OUT,">./result.xls";
open LOC,$ARGV[1] or die "failed to read FILE";

#get chr list (for analyze) 
my %chr;
while(<LOC>){
	chomp;
	my @array = split(/\t/,$_);
	$chr{$array[0]} = $array[0];
	#start_end_(+/-)
	$chr{$array[0]}{$array[4]} = $array[1]."_".$array[2]."_".$array[3];
}
my ($basename)  = $ARGV[0] =~ /.*\/(.*)/;
open OUT,">./$1.result.xls" or die "failed to write output file";
#initialize:analysis chr one by one
foreach my $value(keys %chr){
	#get chr every loc message on chr:$value
	print "$ARGV[0]/$value.cout";
	open CHR,"$ARGV[0]/$value.cout" or die "failed to open chr $value file";
	my %loc_exist;
	my %loc_checkme;
	my %loc_direct;
	while(<CHR>){
		chomp;
		my @array = split (/\t/,$_);
		#deep filter (suggest + oppose)>=4
		if($array[6] + $array[7] >= 4){
			$loc_exist{$array[1]} = 1;
			$loc_direct{$array[1]} = $array[2];
			if($array[6] ne "0"){
				$loc_checkme{$array[1]} = 1;
			}else{
				$loc_checkme{$array[1]} = 0;
			}
		}
	}
	close CHR;
	#calculate every gene on chr:$value----chr
	foreach my $gene(keys %{$chr{$value}}){
		my @array = split(/_/,$chr{$value}{$gene});
#genbody
		my $sum = 0;
		my $sum_valid = 0;
		for(my $i=$array[0];$i <= $array[1];$i++){
			if($loc_exist{$i}){
				#calculate same direction
				if($loc_direct{$i} eq $array[2]){
					$sum++;
					$sum_valid = $sum_valid + $loc_checkme{$i};
				}
			}
		}
		my $genebody = "";
		if($sum>0){
			$genebody = substr(100*$sum_valid/$sum,0,5)."%($sum_valid/$sum)";
		}else{
			$genebody = "NA($sum_valid/$sum)";
		}
#upstream2K
		my $sum_up = 0;
		my $sum_valid_up = 0;
		if($array[2] eq "+"){
			for(my $i=$array[0]-2000;$i < $array[0];$i++){
				if($loc_exist{$i}){
					if($loc_direct{$i} eq $array[2]){
						$sum_up++;
						$sum_valid_up = $sum_valid_up + $loc_checkme{$i};
					}
				}
			}
		}else{
			for(my $i=$array[1]+1;$i < $array[1]+2000;$i++){
				if($loc_exist{$i}){
					if($loc_direct{$i} eq $array[2]){
						$sum_up++;
						$sum_valid_up = $sum_valid_up + $loc_checkme{$i};
					}
				}
			}
		}
		my $upstream;
		if($sum_up>0){
			$upstream = substr(100*$sum_valid_up/$sum_up,0,5)."%($sum_valid_up/$sum_up)";
		}else{
			$upstream = "NA($sum_valid_up/$sum_up)";
		}
#downstream2K
		my $sum_down = 0;
		my $sum_valid_down = 0;
		if($array[2] eq "+"){
			for(my $i=$array[1]+1;$i < $array[1]+2000;$i++){
				if($loc_exist{$i}){
					if($loc_direct{$i} eq $array[2]){
						$sum_down++;
						$sum_valid_down = $sum_valid_down + $loc_checkme{$i};
					}
				}
			}
		}else{
			for(my $i=$array[0]-2000;$i < $array[0];$i++){
				if($loc_exist{$i}){
					if($loc_direct{$i} eq $array[2]){
						$sum_down++;
						$sum_valid_down = $sum_valid_down + $loc_checkme{$i};
					}
				}
			}
		}
		my $downstream;
		if($sum_down>0){
			$downstream = substr(100*$sum_valid_down/$sum_down,0,5)."%($sum_valid_down/$sum_down)";
		}else{
			$downstream = "NA($sum_valid_down/$sum_down)";
		}
		print OUT $gene."\t".$upstream."\t".$genebody."\t".$downstream."\n";
	}
}

