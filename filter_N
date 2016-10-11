#! /usr/bin/perl

use utf8;
use strict;
use warnings;
#use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;

die "perl $0 <fa>\nused to filter fa file N rate>80%" unless(@ARGV eq 1);

my $in = Bio::SeqIO->new(-file=>"< $ARGV[0]", -format=>"fasta");
open  FLIST,">./filter.list" or die "failed to write out filter list";
while (my $seq = $in->next_seq()){
	my $sequence = $seq->seq;
	my $seq_id = $seq->id;
	my $seq_length = $seq->length;
	my @array = split(//,$sequence);
	my $n = 0;
	foreach my $values(@array){
		if($values eq "N" || $values eq "n"){
			$n++;
		}
	}
	if($n/$seq_length < 0.8){
		print ">".$seq_id."\n".$sequence."\n";
	}else{
		print FLIST $seq_id."\n";
	}
}


