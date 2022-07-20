#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
#use lib "$FindBin::Bin";
use FindBin;
my $EXE = $FindBin::RealScript;
my $hmmout = $ARGV[0];

my %seqHash = ();
my $in = Bio::SearchIO->new(-format => 'hmmer',
                           -file   =>$hmmout);
my $crispr = 0;

while( my $result = $in->next_result ) {
    
    # this is a Bio::Search::Result::HMMERResult object
    my $qname = $result->query_name();
    my $qacc = $result->query_accession();
    
    #print "$qname for HMM ", $result->hmm_name(), "\n";
    
    while( my $hit = $result->next_hit ) {
	my $hitname = $hit->name();
        #print "$hitname\n";
	$crispr++;
        while( my $hsp = $hit->next_hsp ) {
            #print "length is ", $hsp->length();
	    #print ", qstart is ", $hsp->query->start;
	    #print ", qend is ", $hsp->query->end;
	    #print ", hstart is ", $hsp->hit->start;
	    #print ", hend is ", $hsp->hit->end;
	    #print ", score is ", $hsp->score;
	    #print ", evalue is ", $hsp->evalue, "\n";
	    my $strand = ".";
	   
	    push @{$seqHash{$hitname}{FEATURE}}, Bio::SeqFeature::Generic->new(
		-primary    => "casgenes",
		-seq_id     => $hitname,
		-source     => $EXE,
		-start      => $hsp->hit->start,
		-end        => $hsp->hit->end,
		-strand     => $strand,
		-score      => $hsp->evalue,
		-frame      => ".",
		-tag        => {
		    "name" => $qname,
		    "accession" => $qacc,
		    "ID" => "CASGENE\_$crispr"
		}
		);
	    
	}
    }
}
my $gffver = 3;
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
print "##gff-version $gffver\n";
#for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
for my $sid (sort keys %seqHash) { 
    for my $f ( sort { $a->start <=> $b->start } @{ $seqHash{$sid}{FEATURE} }) {
	print $f->gff_string($gff_factory),"\n";
    }
}
