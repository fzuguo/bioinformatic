#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
die "perl $0 <FAfile> <FAoutfile>" if @ARGV ne 2;
my %existhash;
open FAOUT,">./$ARGV[1]" or die "failed to write output file";
my $in = Bio::SeqIO->new(-file=>"< $ARGV[0]", -format=>"fasta");
while (my $seq = $in->next_seq()){
	my $seq_id = $seq->id;
	my $seq_fa = $seq->seq;
	if($existhash{$seq_id}){
	}else{
		print FAOUT ">".$seq_id."\n".$seq_fa."\n";
		$existhash{$seq_id} = 1;
	}
}
