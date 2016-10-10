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

#die "perl $0 <arg1> <arg2> <arg3>\n" unless(@ARGV eq 3);
die "perl $0 <fa>\nused to count fa length" unless(@ARGV eq 1);
my $in = Bio::SeqIO->new(-file=>"< $ARGV[0]", -format=>"fasta");
while (my $seq = $in->next_seq()){
	my $seq_id = $seq->id;
	my $seq_length = $seq->length;
	print $seq_id."\t".$seq_length."\n";
}
