#!/usr/bin/env perl
use strict;
use CMD::tsv qw( read_tsv );
use Data::Dumper;
use Getopt::Long;
use JSON;

## Main ##
my %opt;    
GetOptions( \%opt, 
			'base=s', 
			'fusioncatcher=s', 
			'starfusion=s', 
			'fuseq=s', 
			'jaffa=s',
            'arriba=s', 
			'priority=s' );

die "--priority must be given" unless $opt{priority};

my @caller_priority = split /,/, $opt{priority};

# Read blacklist
my %blacklist = read_blacklist("/data/bnf/ref/fusion_blacklist");

# Read the fusion caller output files
my %aggregated;
read_fusioncatcher( $opt{fusioncatcher}, \%aggregated ) if $opt{fusioncatcher};
read_starfusion( $opt{starfusion}, \%aggregated ) if $opt{starfusion};
read_arriba( $opt{arriba}, \%aggregated ) if $opt{arriba};
read_fuseq( $opt{fuseq}, \%aggregated ) if $opt{fuseq};
read_jaffa( $opt{jaffa}, \%aggregated ) if $opt{jaffa};

# Select one representative "best" call for the fusion.
foreach my $genes (keys %aggregated ) {
    my @fus = @{ $aggregated{$genes} };

    # Find "best" call for each caller (most supporting reads)
    my %best;
    foreach my $f ( @fus ) {
		if( !$best{$f->{caller}} or ($f->{spanreads}+$f->{spanpairs}) > ($best{$f->{caller}}->{spanreads}+$best{$f->{caller}}->{spanpairs}) ) {
	    $best{$f->{caller}} = $f;
		}
    }

    # Select the caller with highest priority
    foreach my $caller ( @caller_priority ) {
		if( $best{$caller} ) {
	    	$best{$caller}->{selected} = 1;
	    last;
		}
    }
}

my @data;
foreach my $genes ( keys %aggregated ) {

    my $f;
    $f->{genes} = $genes;
    my @g = split /\^/, $genes;
    $f->{gene1} = $g[0];
    $f->{gene2} = $g[1];
    $f->{calls} = \@{$aggregated{$genes}};
    $f->{blacklisted}=1 if $blacklist{$genes};
    push( @data, $f );
}

# Print json to STDOUT
my $json = JSON->new;
print $json->pretty->encode(\@data);

## Sub routines used in the Main script

sub noversion {
    my $gene = shift;
    $gene =~ s/\.\d+$//;
    return $gene;
}


sub read_blacklist {
    my $fn = shift;

    open( BL, $fn );
    my %bl;
    while( <BL> ) {
	chomp;
	$bl{$_} = 1;
    }
    return %bl;
}


sub read_fusioncatcher {
    my $fn = shift;
    my $agg = shift;

    my @fcatcher = read_tsv($fn);
    
    foreach my $fus ( @fcatcher ) {
	# Get gene symbol pair
	my $gene1 = $fus->{'Gene_1_symbol(5end_fusion_partner)'};
	my $gene2 = $fus->{'Gene_2_symbol(3end_fusion_partner)'};
	my $genes = noversion($gene1)."^".noversion($gene2);
	
	my %fusion_info;
	
	# Get breakpoints 
	$fusion_info{breakpoint1} = $fus->{'Fusion_point_for_gene_1(5end_fusion_partner)'};
	$fusion_info{breakpoint2} = $fus->{'Fusion_point_for_gene_2(3end_fusion_partner)'};
	
	# Get spanning reads and pairs
	$fusion_info{spanreads} = $fus->{Spanning_unique_reads};
	$fusion_info{spanpairs} = $fus->{Spanning_pairs};
	$fusion_info{longestanchor} = $fus->{Longest_anchor_found};
	$fusion_info{commonreads} = $fus->{Counts_of_common_mapping_reads};
	
	$fusion_info{desc} = $fus->{Fusion_description};
	$fusion_info{effect} = $fus->{Predicted_effect};
	
	$fusion_info{caller} = 'fusioncatcher';
	
	unless( $fusion_info{desc} =~ /banned/ ) {
	    push @{$agg->{$genes}}, \%fusion_info;
	}
    }
}


sub read_arriba {
	my $fn = shift;
    my $agg = shift;

    my @fcatcher = read_tsv($fn);


    foreach my $fus ( @fcatcher ) {
	# Get gene symbol pair
	my $gene1 = $fus->{'#gene1'};
	my $gene2 = $fus->{'gene2'};
	my $genes = noversion($gene1)."^".noversion($gene2);
	my %fusion_info;
	
	# Get breakpoints 
	$fusion_info{breakpoint1} = $fus->{'breakpoint1'};
	$fusion_info{breakpoint2} = $fus->{'breakpoint2'}; 
	my $read1 =  $fus->{'split_reads1'};
	my $read2 = $fus->{'split_reads2'};

	# Get spanning reads and pairs
	$fusion_info{spanreads} = $read1 + $read2;
	$fusion_info{spanpairs} = 0;
    $fusion_info{longestanchor} = ( $read1 + $read2 > 25 ? ">25" : "<25" );
	
	$fusion_info{desc} = $fus->{'confidence'};
	$fusion_info{effect} = $fus->{'reading_frame'};
	
	$fusion_info{caller} = 'arriba';
	
	unless( $fusion_info{desc} =~ /banned/ ) {
	    push @{$agg->{$genes}}, \%fusion_info;
	}
    }	
}

sub read_starfusion {
    my $fn = shift;
    my $agg = shift;
    
    my @starf    = read_tsv($fn);    
    foreach my $fus ( @starf ) {

		# Get gene symbol pair
		my $gene1 = (split /\^/, $fus->{LeftGene})[0];
		my $gene2 = (split /\^/, $fus->{RightGene})[0];
		my $genes = noversion($gene1)."^".noversion($gene2);
		my %fusion_info;

		# Get breakpoints 
		$fus->{LeftBreakpoint} =~ s/chr//;
		$fus->{RightBreakpoint} =~ s/chr//;
		$fusion_info{breakpoint1} = $fus->{LeftBreakpoint};
		$fusion_info{breakpoint2} = $fus->{RightBreakpoint};
		
		# Get spanning reads and pairs
		$fusion_info{spanreads} = $fus->{JunctionReadCount};
		$fusion_info{spanpairs} = $fus->{SpanningFragCount};
		$fusion_info{FFPM} = $fus->{FFPM};
		$fusion_info{longestanchor} = ( $fus->{LargeAnchorSupport} eq "YES_LDAS" ? ">25" : "<25" );
	
		$fusion_info{caller} = 'starfusion';
		push @{$agg->{$genes}}, \%fusion_info;
    }
}


sub read_fuseq {
    my $fn = shift;
    my $agg = shift;
    
    my @fuseq = read_tsv($fn, "_ARRAY", "lenient");    
    foreach my $fus ( @fuseq ) {

	# Get gene symbol pair
	my $gene1 = $fus->{symbol5};
	my $gene2 = $fus->{symbol3};
	my $genes = noversion($gene1)."^".noversion($gene2);
	
    
	my %fusion_info;
	
	# Get breakpoints 
	$fusion_info{breakpoint1} = $fus->{chrom5}.":".$fus->{brpos5}.":".$fus->{strand5};
	$fusion_info{breakpoint2} = $fus->{chrom3}.":".$fus->{brpos3}.":".$fus->{strand3};
	
	# Get spanning reads and pairs
	$fusion_info{spanreads} = $fus->{'MR.passed'};
	$fusion_info{spanpairs} = $fus->{'SR.passed'};
#	$fusion_info{supportreads} = $fus->{supportReads};

	$fusion_info{desc} = ($fus->{info} or '');	
	
	$fusion_info{caller} = 'fuseq';
	push @{$agg->{$genes}}, \%fusion_info;
    }
}


    
sub read_jaffa {
    my $fn = shift;
    my $agg = shift;

    system("sed 's/,/\t/g' $fn > $fn.tsv.tmp");
    system("sed -i 's/\"//g' $fn.tsv.tmp");
    
    my @jaffa = read_tsv($fn.".tsv.tmp", "_ARRAY", "lenient");    
    foreach my $fus ( @jaffa ) {

	
	# Get gene symbol pair
	my @genes = split /:/, $fus->{'fusion genes'};
	my $gene1 = $genes[0];
	my $gene2 = $genes[1];
	my $genes = noversion($gene1)."^".noversion($gene2);

	next if $fus->{classification} ne "HighConfidence" and $gene1 ne "MECOM" and $gene2 ne "MECOM";

    
	my %fusion_info;
	
	# Get breakpoints 
	$fusion_info{breakpoint1} = $fus->{chrom1}.":".$fus->{base1}.":".$fus->{strand1};
	$fusion_info{breakpoint2} = $fus->{chrom2}.":".$fus->{base2}.":".$fus->{strand2};
	
	# Get spanning reads and pairs
	$fusion_info{spanreads} = $fus->{'spanning reads'};
	$fusion_info{spanpairs} = $fus->{'spanning pairs'};

	$fusion_info{effect} = "in-frame" if $fus->{inframe} eq "TRUE";
	$fusion_info{desc} = "mitelman" if $fus->{known} eq "Yes";	
	
	$fusion_info{caller} = 'jaffa';
	push @{$agg->{$genes}}, \%fusion_info;
    }
}