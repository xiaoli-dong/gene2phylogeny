#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use threads;
use threads::shared;
use Time::Piece;
use Scalar::Util qw(openhandle);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SearchIO;
use List::Util qw(min max sum);
#use Data::Printer;
#...........................................................................
# Global variabies
my $starttime = localtime;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.1";
my $EXE = $FindBin::RealScript;
my $bin = $FindBin::Bin;


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options

my(@Options, $quiet, $debug,$force,$qtype, $outdir, $len, $cpus, $evalue);

setOptions();

my $logfile = "$outdir/$EXE.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

my $BIDEC = '(\d+\.\d+)';  # pattern of NN.NN for versions that can be compared


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare fasta, remove short contig, replace ambiguous with Ns
my $query = shift @ARGV or err("Please supply a query file: dna/aa fasta or hmm file");
my $contig = shift @ARGV or err("Please supply a contig fasta file on the command line.");

my %contigs = ();
open(CONTIG, "$contig") or die "Could not open $contig to read. $!\n";

$/ = "\n>";
while(<CONTIG>){
    chomp;
    next if !/\S/;
    if(my ($seqname,$seq) =  /^>?(\S+).*?\n(.*)/s){
	$seq =~ tr/ \r\n\t//d;
	my $seqlen = length($seq);
	$contigs{$seqname} = $seqlen;
    }
}
close(CONTIG);
$/ = "\n";

my $db = "";
my %genesDNA = ();

if($qtype eq "hmm"){
    $db = "$outdir/$contig.6frame.translated.faa";
    my $cmd = "";
    if(-e $db){
	$cmd = "hmmsearch --domtblout $outdir/hmm.domtblout -A $outdir/hmm.align -E $evalue --notextw --cpu $cpus -o /dev/null $query $db";
    }else{
	 #translate($contig, $db);
	system("$^X $bin/translate.pl -i $contig > $db");
	$cmd = "hmmsearch --domtblout $outdir/hmm.domtblout -A $outdir/hmm.align -E $evalue --notextw --cpu $cpus -o /dev/null $query $db";
    }
    if(! -e "$outdir/hmm.domtblout"){
	runcmd($cmd);
    }
    open(HMMOUT, "$outdir/hmm.domtblout");
    my %aacontig = ();

#PairedContig_480034_6 -          43388 NiFeSe_Hases         PF00374.15   506  2.2e-225  756.7   0.0   1   1  3.8e-229  3.9e-225  755.8   0.0     1   506 34184 34702 34184 34702 0.99 -
    my $cascount = 0;
while(<HMMOUT>){
	chomp;
	next if /^#/;
	my @line = split(/\s+/, $_);
	my ($contigid, $frame) = $line[0] =~ /(\S+?)\_(\d+)$/;
	my $envStart = $line[19];
	my $envEnd = $line[20];
	my $contiglen = $contigs{$contigid};
	my $qname = $line[3];
	my $acc = $line[4];
	my $qlen = $line[5];
	my $fevalue = $line[6];
	my $cevalue = $line[11];
	next if $cevalue > $evalue;
	my $ievalue = $line[12];
	$aacontig{$line[0]}->{LENGTH} = $line[2];
	my $dnaStart = 0;
	my $dnaEnd = 0;
	my $dnaLen = 0;
	my $domainid = $line[9];
        #get DNA coordinates:
    my $strand = "+";
	#forward
	if($frame >= 1 && $frame <=3){
	    $dnaStart = $envStart * 3 - 2 + ($frame - 1);
	    $dnaEnd = $envEnd * 3  + ($frame - 1);
	    $dnaLen = $dnaStart - $dnaEnd + 1;

	}

	#reverse complement strand
	elsif($frame >= 4 && $frame <=6){
	    $dnaEnd = $contiglen - ($envStart * 3 - 2 + ($frame - 4)) + 1;
	    $dnaStart = $contiglen - ($envEnd * 3  + ($frame - 4)) + 1;
	    $dnaLen = $dnaStart - $dnaEnd + 1;
       $strand = "-";
	}
	$cascount++;
	my $gff_out_string = "$contigid\t$EXE\tcasgenes\t$dnaStart\t$dnaEnd\t$evalue\t$strand\t\.\tID=cas_$cascount;qname=$qname;accession=$acc;cevalue=$cevalue;ievalue=$ievalue";
	push(@{$aacontig{$line[0]}->{qname}->{$acc}->{DOMAIN}}, $gff_out_string);
	push(@{$aacontig{$line[0]}->{qname}->{$acc}->{DOMAIN_DNA}}, $gff_out_string);
	$aacontig{$line[0]}->{qname}->{$acc}->{LENGTH} = $qlen;
    }
     open(GFF, ">$outdir/hmm.dna_coordinate.gff.txt");

    foreach my $contigid (keys %aacontig){
	foreach my $myqname (keys %{$aacontig{$contigid}->{qname}}){
	    foreach my $domain (sort @{$aacontig{$contigid}->{qname}->{$myqname}->{DOMAIN_DNA}}){
		print GFF $domain, "\n";
	    }
	}
    }
    close(GFF);

    open(DB, "$db");
    open(AA, ">$outdir/hmm.hit.aa.fasta");
    $/ = "\n>";
    while(<DB>){
	chomp;

	chomp;

	if(my ($seqname,$seq) =  /^>?(\S+).*?\n(.*)/s){
	    $seq =~ tr/ \r\n\t//d;

	    if(exists $aacontig{$seqname}){
		foreach my $myqname (keys %{$aacontig{$seqname}->{qname}}){
		    foreach my $domain (@{$aacontig{$seqname}->{qname}->{$myqname}->{DOMAIN}}){
			my ($start, $end) = $domain =~ /(\d+?)\:(\d+)/;
			my $length = $end - $start + 1;

			if($length >= $len){
			    my $subseq = substr($seq,$start-1,$length);
			print  AA ">$seqname\_$start\_$end gene=$myqname aalength=$length\n$subseq\n";
		    }
		}
		}
	    }
	}
    }
    $/ = "\n";
    close(DB);
    close(AA);


    open(CONTIG, "$contig");
    open(DNA, ">$outdir/hmm.hit.dna.fasta");
    $/ = "\n>";
    while(<CONTIG>){
	chomp;

	if(my ($seqname,$seq) =  /^>?(\S+).*?\n(.*)/s){
	    $seq =~ tr/ \r\n\t//d;
	    my $contiglen = length($seq);

	    foreach my $f (1..6){
		my $dbEntryName = "$seqname\_$f";
		if(exists $aacontig{$dbEntryName}){
		    my $index = 0;
           foreach my $myqname (keys %{$aacontig{$dbEntryName}}){

		    foreach my $domain (@{$aacontig{$dbEntryName}->{qname}->{$myqname}->{DOMAIN_DNA}}){
			my $aa_domain = ${$aacontig{$dbEntryName}->{qname}->{$myqname}->{DOMAIN}}[$index++];
			print STDERR "$dbEntryName\tdna_domain_coord=$domain\taa_domain_coord=$aa_domain\n";
                        #AA coordinates
			my ($start, $end) = $domain =~ /(\d+?)\:(\d+)/;
			my ($aastart, $aaend) = $aa_domain =~ /(\d+?)\:(\d+)/;
			my $dnaLength = $end - $start + 1;

			if($dnaLength >= $len*3){
			    my $subseq = substr($seq,$start-1,$dnaLength);
			    print  DNA ">$seqname\_$f\_$aastart\_$aaend gene=$myqname dnalength=$dnaLength\n$subseq\n";
			}
		    }
		}}
	    }
	}
    }
    $/ = "\n";
    close(CONTIG);
    close(DNA);
}



##gff-version 3

#NODE_1256            -          euk_18SrRNA          -              894    1848    1110     259    1126     256   16974    -     1.2e-96  320.9  12.7  -


#my $gff_out_string = "$seqid\t$EXE\t$gene\t$begin\t$end\t$evalue\t$strand\t\.\tID=$ID;Name=$gene\n";



my $endtime = localtime;
my $walltime = $endtime - $starttime;
#msg("Walltime used:", $walltime->pretty);  # Heikki says this method only came with 1.20
my $pretty = sprintf "%.2f minutes", $walltime->minutes;
msg("Walltime used: $pretty");
msg($walltime % 2 ? "Share and enjoy!" : "Thank you, come again.");
#EXIT


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
    use Getopt::Long;

    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},
	{OPT=>"quiet!",  VAR=>\$quiet, DEFAULT=>0, DESC=>"No screen output"},
	{OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug mode: keep all temporary files"},
	'Inputs:',
	{OPT=>"qtype=s",  VAR=>\$qtype, DEFAULT=>'dna', DESC=>"query type: dna, aa or hmm"},
	'Outputs:',
	{OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'', DESC=>"Output folder [auto]"},
   {OPT=>"len=s",  VAR=>\$len, DEFAULT=>100, DESC=>"gene length cutoff"},
	'Computation:',
	{OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>8, DESC=>"Number of CPUs to use"},
	{OPT=>"evalue=f",  VAR=>\$evalue, DEFAULT=>1E-2, DESC=>"Similarity e-value cut-off"},
	);

    (!@ARGV) && (usage());

    &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @Options) || usage();

    # Now setup default values.
    foreach (@Options) {
	if (ref $_ && defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	    ${$_->{VAR}} = $_->{DEFAULT};
	}
    }
}

#----------------------------------------------------------------------
sub usage {
    print STDERR
	"Name:\n  ", ucfirst($EXE), " $VERSION by $AUTHOR\n",
	"Synopsis:\n  Detecting genes from genomic/metagenomic assembles\n",
	"Usage:\n  $EXE [options] <dna/aa in fasta or hmm profile> <contigs.fasta>\n";
    foreach (@Options) {
	if (ref) {
	    my $def = defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	    $def = ($def ? ' (default OFF)' : '(default ON)') if $_->{OPT} =~ m/!$/;
	    my $opt = $_->{OPT};
	    $opt =~ s/!$//;
	    $opt =~ s/=f$/ [n.n]/;
	    printf STDERR "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
	}
	else {
	    print STDERR "$_\n";
	}
    }
    exit(1);
}

#----------------------------------------------------------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

#----------------------------------------------------------------------

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}

#----------------------------------------------------------------------

sub delfile {
  for my $file (@_) {
    if ($debug) {
      msg("In --debug mode, saving temporary file:", $file);
    }
    else {
      msg("Deleting unwanted file:", $file);
      unlink $file or warn "Could not unlink $file: $!\n";;
    }
  }
}

#----------------------------------------------------------------------

sub version {
  print STDERR "$EXE $VERSION\n";
  exit;
}

#----------------------------------------------------------------------

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}



#----------------------------------------------------------------------
sub find_exe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe;
  }
  return;
}
