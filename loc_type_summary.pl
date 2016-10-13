#!/usr/bin/perl
use warnings;
use strict;
die "perl $0 <GTF_gene_info> <accept_bam> <RUN_log> > result.txt" if @ARGV ne 3;
my %gene_loc_hash;
open GTF,$ARGV[0] or die "failed to open GTF file";
while(<GTF>){
	chomp;
	my @array = split(/\t/,$_);
	my ($biotype) = $_ =~ /gene_biotype "(.*?)"/;
	my $loc = $array[3]."_".$array[4];
	#defined loc type
	$gene_loc_hash{"chr_".$array[0]}{$loc} = $biotype;
}
#defined type hash
my %type_hash;
open BAM,"samtools view $ARGV[1] |" or die "failed to open BAM file";
open LOG,">./$ARGV[2]" or die "failed to open log file";
while(<BAM>){
	chomp;
	my @array2 = split(/\t/,$_);
	my $mark = "unknown";
	foreach my $keys(keys $gene_loc_hash{"chr_".$array2[2]}){
		#get start & end loc;
		my ($loc_start,$loc_end) = $keys =~ /(.*?)_(.*)/;
		if($array2[3] >= $loc_start && $array2[3] <= $loc_end){
			$type_hash{$gene_loc_hash{"chr_".$array2[2]}{$keys}}++;
			$mark = $gene_loc_hash{"chr_".$array2[2]}{$keys};
			last;
		}
	}
	print LOG $array2[0]."\t$mark\n";
}
#result output
foreach my $keys(keys %type_hash){
	print $keys."\t".$type_hash{$keys}."\n";
}
