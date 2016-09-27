#!/usr/bin/perl
use warnings;
use strict;

die "perl $0 <gene_loc_info> <GTF> <new_GTF> <prefix>" if @ARGV ne 4;
######
# used to get domain locus on chromosome. 
# <gene_loc_info> transcript_name \t start_acid_locus \t end_acid_locus
# <GTF> original gtf 
# <new_GTF> outputfile contain Domain locus info
# <prefix> user defined  domain name 
######


#get CDS loc info & direction info from gtf file
my %CDS_loc;
my %CDS_direc;
my %CDS_chr;
my %CDS_source;
my %CDS_geneid;
my %CDS_start_code;
open GTF,$ARGV[1] or die "failed to open GTF file";
while(<GTF>){
	chomp;
	my @array = split (/\t/,$_);
	if($array[2] =~ /CDS/i){
		my($trans) = $array[8] =~ /transcript_id "(.*?)"/;
		if($CDS_loc{$trans}){
			$CDS_loc{$trans} = $CDS_loc{$trans}.$array[3]."_".$array[4]."/";
			if($CDS_direc{$trans} eq "+"){
				$CDS_start_code{$trans} = ($array[3] + $array[7]) if (($array[3] + $array[7]) < $CDS_start_code{$trans});
			}else{
				$CDS_start_code{$trans} = ($array[4] - $array[7]) if(($array[4] - $array[7]) > $CDS_start_code{$trans});
			}
		}else{
			$CDS_loc{$trans} = $array[3]."_".$array[4]."/";
			$CDS_direc{$trans} = $array[6];
			$CDS_chr{$trans} = $array[0];
			$CDS_source{$trans} = $array[1];
			$array[8] =~ /gene_id "(.*?)"/;
			$CDS_geneid{$trans} = $1;
			if($CDS_direc{$trans} eq "+"){
				$CDS_start_code{$trans} = $array[3] + $array[7];
			}else{
				$CDS_start_code{$trans} = $array[4] - $array[7];
			}
		}
	}
}
#finsh get CDS loc info

#sort CDS loc hash and deal start code
foreach my $keys(keys %CDS_loc){
	my @array = split(/\//, $CDS_loc{$keys});
	my @new_array = sort{(split /_/,$a)[0] cmp (split /_/,$b)[0]} @array;
	if($CDS_direc{$keys} eq "+"){
		my @s_e = split(/_/,$new_array[0]);
		if($s_e[0] > $CDS_start_code{$keys}){
			$new_array[0] = $CDS_start_code{$keys}."_".$s_e[1];
			print "$keys start code change to $CDS_start_code{$keys}\n";
		}
	}else{
		my @s_e = split(/_/,$new_array[-1]);
		if($s_e[1] < $CDS_start_code{$keys}){
			$new_array[-1] = $s_e[0]."_".$CDS_start_code{$keys};
			print "$keys start code change to $CDS_start_code{$keys}\n";
		}
	}
	my $string = join("/",@new_array);
	$CDS_loc{$keys} = $string;
}
#finish sort

#finsh

#start deal protein align file
open  ALIGN,$ARGV[0] or die "failed to open protein_align_info file";
my %domain;
while(<ALIGN>){
	chomp;
	my @array = split(/\t/, $_);
	my $locus = &get_Domain_info($array[0],$array[1],$array[2]);
	$domain{$array[0]} = $locus;
}
#finish deal align file

#begin write out new gtf
open NEWGTF,">./$ARGV[2]" or die "failed to write out new gtf file";
seek GTF,0,0;
my $flag  = "";
my $last_line = "";
while(<GTF>){
	if($_ =~ /transcript_id "(.*?)"/){
		if($flag ne $1 && $flag ne ""){
			if($domain{$1}){
				my @array = split(/\//,$domain{$1});
				foreach my $values(@array){
					my @s_e = split(/_/,$values);
					print NEWGTF $CDS_chr{$1}."\t".$CDS_source{$1}."\t$ARGV[3]\t".$s_e[0]."\t".$s_e[1]."\t.\t".$CDS_direc{$1}."\t.\ttranscript_id \"$1\"; gene_id \"$CDS_geneid{$1}\";\n";
				}
			delete $domain{$1};
			}
		}
	$flag = $1;
	$last_line = $_;
	}
	print NEWGTF $_;
}
if($last_line =~ /transcript_id "(.*?)"/){
	if($domain{$1}){
		my @array = split(/\//,$domain{$1});
		foreach my $values(@array){
			my @s_e = split(/_/,$values);
			print NEWGTF $CDS_chr{$1}."\t".$CDS_source{$1}."\t$ARGV[3]\t".$s_e[0]."\t".$s_e[1]."\t.\t".$CDS_direc{$1}."\t.\ttranscript_id \"$1\"; gene_id \"$CDS_geneid{$1}\";\n";
		}
	}
}

#sub programs

#used for get Domain loc info
sub get_Domain_info{
	my $trans_name = shift;
	my $trans_start = shift;
	my $trans_end = shift;
	#if transcription is not exists return;
	if(!$CDS_loc{$trans_name} or !$CDS_direc{$trans_name}){
		return "unknown transcript or direction erro!";
	}
	
	my $locus_start = "0";
	my $locus_end = "0";
	my $zone = "";
	if($CDS_direc{$trans_name} eq "+"){
		$locus_start = &plus($CDS_loc{$trans_name},$trans_start);
		$locus_end = &plus($CDS_loc{$trans_name},$trans_end) + 2;
		$zone = &newzone($CDS_loc{$trans_name},$locus_start,$locus_end);
	}elsif($CDS_direc{$trans_name} eq "-"){
		$locus_start = &minus($CDS_loc{$trans_name},$trans_start);
		$locus_end = &minus($CDS_loc{$trans_name},$trans_end) - 2;
		$zone = &newzone($CDS_loc{$trans_name},$locus_end,$locus_start);
	}
	return $zone;
}


#used to confirm loc on plus chain
sub plus{
	my $zone = shift;
	my $p_loc = shift;
	my $length = ($p_loc-1)*3;
	my @array = split(/\//,$zone);
	foreach my $values(@array){
		my @start_end = split(/\_/,$values);
		my $last_length = $length;
		$length = $length - ($start_end[1] - $start_end[0] + 1);
		if($length < 0){
			return ($last_length+$start_end[0])
		}
	}
}

#used to confirm loc on minus chain
sub minus{
	my $zone = shift;
	my $p_loc = shift;
	my $length = ($p_loc-1)*3;
	my @array = split(/\//,$zone);
	my $i = @array;
	for($i = @array-1; $i >= 0; $i--){
		my @start_end = split(/\_/,$array[$i]);
		my $last_length = $length;
		$length = $length - ($start_end[1] - $start_end[0] + 1);
		if($length < 0){
			return ($start_end[1] - $last_length);
		}
	}
}

#generate newzone
sub newzone{
	my $zone = shift;
	my $start = shift;
	my $end = shift;
	my @array = split(/\//,$zone);
	my @newarray;
	for(my $i=0; $i<@array; $i++){
		my @s_e = split (/_/,$array[$i]);
		if($s_e[1] < $start){$array[$i] = ""};
		if($s_e[0] > $end){$array[$i] = ""};
		if($start > ($s_e[0]-1) && $start <= ($s_e[1]+1)){
			$array[$i] = $start."_".$s_e[1];
			@s_e = split (/_/,$array[$i]);
		}
		if($end > ($s_e[0]-1) && $end <= ($s_e[1]+1)){
			$array[$i] = $s_e[0]."_".$end;
		}
		if($array[$i] ne ""){
		 push @newarray,$array[$i];
		}
	}
	my $string = join("/",@newarray);
	return $string;
}
