#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
my ($ref, $fastq1, $fastq2, $rname, $qname);
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
    "bo=s" =>\$bo
    );

($ref && $fastq1 && $fastq2 && $rname && $qname) ||
    die "usage: $0 OPTIONS where options are:\n ".
    "-r  <fasta reference file>\n".
    "-rn <reference name>\n".
    "-1 <fastq read1>\n".
    "-2 <fastq read2>\n".
    "-qn <reads file name>\n".
    "-m <mapping tool bowtie2|bbamp, default is bowtie2>\n".
    "-bo <bowtie2 option >\n".
    "-t <number of threads to use, default is 8>\n".
    "-o <output directory, default is ./>\n";

# -r PF00374.NiFe.hmm.hit.dna.fasta -rn NiFe -1 out1.fastq -2 out2.fastq -qn Encana
mkdir($outdir) unless(-d $outdir);
geneAmount($ref, $threads, $fastq1, $fastq2, $outdir, $rname, $qname, $mtool);

sub geneAmount{

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
    if(! -e "$outdir/$rname\_$qname-smds.coverage"){
	$cmd .= "genomeCoverageBed -ibam $outdir/$rname\_$qname-smds.bam > $outdir/$rname\_$qname-smds.coverage;";
    }
    # generate table with length and coverage stats per contig (From http://github.com/BinPro/CONCOCT)
    #$cmd .= "python gen_contig_cov_per_bam_table.py --isbedfiles $ref $outdir/$rname_$rname-smds.coverage > $outdir/$rname\_$qname-smds.coverage.percontig";
    msg("cmds:", $cmd);
    runcmd($cmd);

#1. chromosome (or entire genome)
#2. depth of coverage from features in input file
#3. number of bases on chromosome (or genome) with depth equal to column 2.
#4. size of chromosome (or entire genome) in base pairs
#5. fraction of bases on chromosome (or entire genome) with depth equal to column 2.
#PairedContig_1054_1_156_524     0       12      1107    0.0108401
#PairedContig_1054_1_156_524     1       25      1107    0.0225836

    #open(FASTA, $ref) or die "Could not open $outdir/$ref to read, $!\n";
    my %genes = ();
    #while(<FASTA>){
#	if(/^>(\S+?)\s+?dnalength=(\d+)/){
#	    $genes{$1}->{bases} = 0;
#	    $genes{$1}->{length} = $2;
#	}
 #   }
  #  close(FASTA);

    open(COV, "$outdir/$rname\_$qname-smds.coverage") or die "Could not open $outdir/$rname\_$qname-smds.coverage to read, $!\n";
    while(<COV>){
	next if /^genome/;
	my @line = split(/\t/, $_);
	#print STDERR "$line[0]\n";
	if(not exists $genes{$line[0]}->{length}){
	    $genes{$line[0]}->{length} = $line[3];
	}
	$genes{$line[0]}->{bases} += $line[1] * $line[2];


    }
    close(COV);

    open(GENECOV, ">$outdir/$rname\_$qname-smds.genecov") or die "Could not open $outdir/$rname\_$qname-smds.genecov to write, $!\n";
    print GENECOV "#geneid\ttotalBases\tgenelength\tcov\n";
    foreach my $gene (keys %genes){
	print GENECOV "$gene\t", $genes{$gene}->{bases}, "\t", $genes{$gene}->{length}, "\t", sprintf ("%.4f",$genes{$gene}->{bases}/$genes{$gene}->{length}), "\n";
    }
    close(GENECOV);

# Remove temp files
    my @tempfiles = ("$outdir/$rname\_$qname.sam", "$outdir/$rname\_$qname.bam", "$outdir/$rname\_$qname-smd.bam", "$outdir/$rname\_$qname-s.bam", "$outdir/$rname\_$qname-s.bam.bai");

    foreach (@tempfiles){
	#delfile($_);
    }

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

