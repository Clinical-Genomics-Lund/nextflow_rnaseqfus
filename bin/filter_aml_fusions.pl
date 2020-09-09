#!/usr/bin/env perl

use strict;
use File::Find::Rule;
use Data::Dumper;

my $dir = shift;

my @files = File::Find::Rule->file->name("final-list_candidate-fusion-genes.txt")->in($dir);

#print Dumper @files;

my @results;
my $header;

foreach(@files){

    my @path = split("/",$_);

    my $sample = $path[-3];
#    print Dumper @path;#die;
    open(my $data, '<', $_) or die "buhu";
    $header=<$data>;
    
    while (my $line = <$data>) {
	chomp($line);
	if ($line=~ /(ABL1|ABL2|ABL|AFF1|AML1|BCR|CBFB|CSF1R|DEK|ELL|ETV6|KMT2A|MLL|MECOM|MLL|MLLT1|MLLT4|MYH11|NUP214|PBX1|PML|RARA|RUNX1|RUNX1|RUNX1T1|RUNX1|TEL|TCF3|PAX|PDGFRB|HLF|FIP1L1|PDGFRA|FGFR1|JAK2|PCM1)/){
	    push(@results,$sample."\t".$line."\n");
	}
    }
}

chomp($header);
print "Sample\t$header\n";
foreach(@results){print;}
