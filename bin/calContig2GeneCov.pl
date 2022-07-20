#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
my ($ref, $fastq1, $fastq2, $rname, $qname, $gff);
my  $threads = 8;
my $JAVA_OPT = "-Xmx20g ";
my $outdir = ".";
my $picard = "/export/data/programs/Picard";
my $mtool = "bowtie2";
my $bo = "";

&GetOptions(
    "r=s" =>\$ref,
    "rn=s" =>\$rname,
    "1=s" =>\$fastq1,
    "2=s" =>\$fastq2,
    "qn=s" =>\$qname,
    "m=s" =>\$mtool,
    "t=i" =>\$threads,
    "o=s" =>\$outdir,
    "g=s" =>\$gff,
    "bo=s" =>\$bo
    );

($ref && $fastq1 && $fastq2 && $rname && $qname && $gff) ||
    die "usage: $0 OPTIONS where options are:\n ".
    "-r  <fasta reference file>\n".
    "-rn <reference name>\n".
    "-1 <fastq read1>\n".
    "-2 <fastq read2>\n".
    "-g <gff file>\n".
    "-qn <reads file name>\n".
    "-m <mapping tool bowtie2|bbamp, default is bowtie2>\n".
    "-bo <bowtie2 option >\n".
    "-t <number of threads to use, default is 8>\n".
    "-o <output directory, default is ./>\n";

# -r PF00374.NiFe.hmm.hit.dna.fasta -rn NiFe -1 out1.fastq -2 out2.fastq -qn Encana
msg("calContig2GeneCov start");
mkdir($outdir) unless(-d $outdir);
my $perbasemfile = doMapping($ref, $threads, $fastq1, $fastq2, $outdir, $rname, $qname, $mtool);
my $genes = gff2geneFeature($gff);
calCov($genes, $perbasemfile);
msg("calContig2GeneCov end");

#clear();
sub doMapping{
    
    #$ref is the reference genes

    my ($ref, $threads, $fastq1, $fastq2, $outdir, $rname, $qname, $mtool, $bo) = @_;
    my $cmd = "";

    if($mtool eq "bowtie2"){
	# Index reference, Burrows-Wheeler Transform
	if(! -e "$outdir/$ref.1.bt2"){
	    $cmd .= "bowtie2-build $ref $outdir/$ref;";
	    #runcmd($cmd);
	}
    }
    # Align Paired end and bam it
    if(! -e  "$outdir/$rname\_$qname.sam"){
	if($mtool eq "bbmap"){
	    $cmd .= "bbmap.sh $JAVA_OPT threads=$threads ambiguous=random ref=$ref in=$fastq1 in2=$fastq2 out=$outdir/$rname\_$qname.sam covstats=$outdir/bbmap_covstats.txt scafstats=$outdir/bbmap_scafstats.txt;";
	}
	else{
	    $cmd .= "bowtie2 -p $threads $bo -x $outdir/$ref -1 $fastq1 -2 $fastq2 -S $outdir/$rname\_$qname.sam;";
	}
    }
    
    # Index reference for samtools
    if(! -e "$ref.fai"){
	$cmd .= "samtools faidx $outdir/$ref;";
    }
    if(! -e "$outdir/$rname\_$qname.bam"){
	#covert sam to bam and filter out the unmapped reads to make the downstream process faster
	$cmd .= "samtools view -F4 -bt $outdir/$ref.fai $outdir/$rname\_$qname.sam > $outdir/$rname\_$qname.bam;";
    }
    if(! -e "$outdir/$rname\_$qname-s.bam"){
	$cmd .= "samtools sort $outdir/$rname\_$qname.bam $outdir/$rname\_$qname-s;";
    }
    if(! -e "$outdir/$rname\_$qname-s.bam.bai"){
	$cmd .= "samtools index $outdir/$rname\_$qname-s.bam;";
    }
    if(! -e "$outdir/$rname\_$qname-smd.bam"){
	# Mark duplicates and sort
	$cmd .= "java $JAVA_OPT -jar $picard/picard.jar MarkDuplicates  ".
	    "INPUT=$outdir/$rname\_$qname-s.bam ".
	    "OUTPUT=$outdir/$rname\_$qname-smd.bam ".
	    "METRICS_FILE=$outdir/$rname\_$qname-smd.metrics ".
	    "AS=TRUE ".
	    "VALIDATION_STRINGENCY=LENIENT ".
	    "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ".
	    "REMOVE_DUPLICATES=TRUE;";
    }

    if(! -e "$outdir/$rname\_$qname-smds.bam.bai"){
	$cmd .= "samtools sort $outdir/$rname\_$qname-smd.bam $outdir/$rname\_$qname-smds;";
	$cmd .= "samtools index $outdir/$rname\_$qname-smds.bam;";
    }
    if(! -e "$outdir/$rname\_$qname-smds.flagstat"){   
	$cmd .= "samtools flagstat $outdir/$rname\_$qname-smds.bam > $outdir/$rname\_$qname-smds.flagstat;";
    }
# Determine Genome Coverage and mean coverage per contig
    if(! -e "$outdir/$rname\_$qname-smds.perbase.coverage"){
	$cmd .= "genomeCoverageBed -ibam $outdir/$rname\_$qname-smds.bam -dz > $outdir/$rname\_$qname-smds.perbase.coverage;";
    }
    # generate table with length and coverage stats per contig (From http://github.com/BinPro/CONCOCT)
    #$cmd .= "python gen_contig_cov_per_bam_table.py --isbedfiles $ref $outdir/$rname_$rname-smds.coverage > $outdir/$rname\_$qname-smds.coverage.percontig";
    msg("cmds:", $cmd);
    runcmd($cmd) if $cmd ne "";
    return "$outdir/$rname\_$qname-smds.perbase.coverage";
    
}

sub calCov{
    
    my ($genes, $perbasemfile) = @_;
    
#The “per-base” output format is as follows:

    #chromosome
    #chromosome position
    #depth (number) of features overlapping this chromosome position.
    #zero-based coordinates
    my %refBaseCovs = ();
    open(COV,$perbasemfile) or die "Could not open $perbasemfile, $!\n";
    while(<COV>){
	next if /^genome/;
	my @line = split(/\t/, $_);
	$refBaseCovs{$line[0]}->{$line[1]} = $line[2];
    }
    close(COV);
    #calculate contig coverage
    open(CONTIGCOV, ">$outdir/$rname\_$qname-smds.contigcov") or die "Could not open $outdir/$rname\_$qname-smds.contigcov to write, $!\n";
    print CONTIGCOV "#contigid\ttotalBases\tgenelength\tcov\n";
    foreach my $contig (keys %refBaseCovs){
	my $tbp = 0;
	my $len = 0;
	foreach my $base (keys %{$refBaseCovs{$contig}}){
	    $len++;
	    $tbp += $refBaseCovs{$contig}->{$base};
	}
	my $cov = sprintf("%.4f", $tbp/$len);
		
	print CONTIGCOV "$contig\t$tbp\t$len\t$cov\n";
    }
    close(CONTIGCOV);

    #calculate gene coverage
    open(GENECOV, ">$outdir/$rname\_$qname-smds.genecov") or die "Could not open $outdir/$rname\_$qname-smds.genecov to write, $!\n";
    print GENECOV "#geneid\tparent\ttotalBases\tgenelength\tcov\n";
    foreach my $gene (sort keys %{$genes}){
	my $start = $genes->{$gene}->{start};
	my $end = $genes->{$gene}->{end};
	my $contig = $genes->{$gene}->{parent};
	my $len = $end - $start + 1;
	my $tbp = 0;
	for my $i ($start...$end){
	    if(exists $refBaseCovs{$contig}->{$i}){
		$tbp += $refBaseCovs{$contig}->{$i};
	    }
	}
	my $cov = sprintf("%.4f", $tbp/$len);
	
	print GENECOV "$gene\t$contig\t$tbp\t$len\t$cov\n";
    }
    close(GENECOV);
}

sub clear{
# Remove temp files
    my @tempfiles = ("$outdir/$rname\_$qname.sam", "$outdir/$rname\_$qname.bam", "$outdir/$rname\_$qname-smd.bam", "$outdir/$rname\_$qname-s.bam", "$outdir/$rname\_$qname-s.bam.bai");

    foreach (@tempfiles){
	delfile($_);
    }

}

sub gff2geneFeature{
    my ($gff) = @_;
    my %genes = ();
    open(GFF, $gff) or die "Could not open $gff to read, $!\n"; 
    
    while(<GFF>){
	next if /^#/;
	my $taxon = "";
	my $geneid = "";
	if(/ID=(\S+?)[;\n]/){
	    $geneid = $1;
	    #print "$geneid*****************$_\n";
	}
	
	if(/taxon=(.*)/){
	    $taxon = $1;
	}
	chomp;
	my @line =split(/\t/, $_);
	my $tags = $line[8];
	if($geneid ne ""){
	    $genes{$geneid}->{parent} = $line[0];
	    $genes{$geneid}->{start} = $line[3];
	    $genes{$geneid}->{end} = $line[4];
	    $genes{$geneid}->{taxon} = $taxon;
	}
    }
    close(GFF);
    return \%genes;
    
}
#----------------------------------------------------------------------

sub runcmd {
    msg("Running:", @_);
    system(@_)==0 or err("Could not run command:", @_);
}

sub msg {
    my $t = localtime;
    my $line = "[".$t->hms."] @_\n";
    #print STDERR $line if openhandle(\*LOG);
    print STDERR $line;
}

#----------------------------------------------------------------------

sub delfile {
    for my $file (@_) {
	#if ($debug) {
	 #   msg("In --debug mode, saving temporary file:", $file);
	#}
	#else {
	msg("Deleting unwanted file:", $file);
	unlink $file or warn "Could not unlink $file: $!\n";;
	#}
    }
}

