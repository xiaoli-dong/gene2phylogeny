#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

open(TAB, $ARGV[0]);
open(FASTA, $ARGV[1]);
my %id2Features = ();
#1  Azotobacter vinelandii  hoxG  AAA82506 P21949   HoxG
while(<TAB>){
    chomp;
    my @line = split(/\s+/, $_);
    $id2Features{$line[4]} = "$line[2]_$line[1]_$line[0]";
    print STDERR "$line[4]\t$line[2]_$line[1]_$line[0]\n";
}


$/ = "\n>";
while(<FASTA>){
    chomp;
    
    if(my ($seqname,$other, $seq) =  /^>?(\w+)(;?.*?)\n(.*)/s){
	print STDERR "match $seqname\n";
        if(exists $id2Features{$seqname}){
            print  ">$seqname\_$id2Features{$seqname}\n$seq\n";
        }

    }

}

$/ = "\n";

close(FASTA);
