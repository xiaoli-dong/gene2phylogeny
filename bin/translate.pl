#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::PrimarySeqI;
use Getopt::Long;

my ($fasta, $id);
$id = 1;

&GetOptions(
    "i=s" =>\$fasta,
    "c=i" =>\$id
    );

($fasta && $id) ||
    die "usage: $0 OPTIONS where options are:\n".
    "-i <fasta dna file>\n".
    "-c <codontable_id: default is 1 standard code>\n".
    "see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c for genetic code\n";

# This is the main sequence object in Bioperl

# gets a sequence from a file
my $seqio  = Bio::SeqIO->new( '-format' => 'fasta' , -file =>$fasta);

#frame 4, 5 6, is the frame 1, 2, 3, of the reverse complement dna 
while (my $seq = $seqio->next_seq) {
    
    foreach my $frame (0..2) {
	my $trans = $seq->translate(-frame=>$frame, -codontable_id =>$id);
	my $seqstr = $trans->seq;
	$seqstr =~ s/\*/X/g;
	my $id = $trans->id;
	print ">$id\_", $frame+1, "\n$seqstr\n";
    }
    
    my $revseq = $seq->revcom;
    
    foreach my $frame (0..2) {
	my $trans = $revseq->translate(-frame=>$frame, -codontable_id =>$id);
	my $seqstr = $trans->seq;
	$seqstr =~ s/\*/X/g;
	my $id = $trans->id;
	print ">$id\_", $frame+4, "\n$seqstr\n";
    }
}

