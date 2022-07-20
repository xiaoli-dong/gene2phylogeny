#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($refGenes, $mGeneFrags, $aligned, $bootstrap, $threads);

$aligned = 0;
$bootstrap = 100;
$threads = 2;

&GetOptions(
    "r=s" =>\$refGenes,
    "m=s" =>\$mGeneFrags,
    "a=i" =>\$aligned,
    "b=i" =>\$bootstrap,
    "t=i" =>\$threads
    );

($refGenes && $mGeneFrags && $aligned) ||
    die "usage: $0 OPTIONS where options are:\n".
    "-r <reference AA gene seqs in fasta file>\n".
    "-m <AA gene fragments from contigs in fasta file>\n".
    "-a <aligned reference gene seqs 0|1, default 0: not aligned>\n".
    "-b <the number of the bootstrap replicates>\n".
    "-t <the number of the threads>\n";

##reference alignment
my $ref = "";
if($aligned == 1){
    $ref = $refGenes;
}
if($aligned == 0){
    $ref = makeMSA($refGenes);
}

#ref + mGeneFrags MSA
my $whole = addSeq2AlignedSeqs($ref, $mGeneFrags);

#multi-step approach utilizing refeence alignments and trees in order to 
#minimize errors and biases introduced by the fragmentary and non-overlapping nature of the metagenome

#build bootstrapped maximum-likelihood tree from reference alignment

use Cwd 'abs_path';
my $abs_path = abs_path($ref);
use File::Basename;
my($file, $dir, $ext) = fileparse($abs_path);
$dir =~ s/\/$//;
print STDERR "$dir\n";
my $cmd = "";
if(! -e "$dir/RAxML_bestTree.reftree"){
    $cmd = "raxmlHPC-PTHREADS  -f a -s $ref -w $dir -n reftree -m PROTGAMMAWAGF -x 0123 -# 100 -p 012345  -T $threads";
    print STDERR "cmd = $cmd\n";
    system ($cmd);
}

$cmd = "raxmlHPC-PTHREADS  -f a -s $whole -n whole -m PROTGAMMAWAGF -x 0123 -# $bootstrap -p 012345  -T $threads -r $dir/RAxML_bestTree.reftree";
print STDERR "cmd = $cmd\n";
system ($cmd);




#get multiple sequence alignment for the reference gene sequences
sub makeMSA{
    my $fasta = @_;
    my $cmd = "mafft --genafpair --maxiterate 1000 $fasta > ref.aligned.fasta";
    
    system ($cmd);
    return "ref.aligned.$fasta";
}

#add the new sequences to the existing alignment and the alignment length is kept unchanged 
sub addSeq2AlignedSeqs{
    my ($refFasta, $fragFasta) = @_;
    my $cmd = "mafft --addfragments $fragFasta --keeplength --reorder --thread 8 $refFasta > all.aligned.fasta";
    print STDERR "$cmd\n";
    system ($cmd);
    return "all.aligned.fasta";
}
