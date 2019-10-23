#!/usr/bin/env perl
use strict;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use CMD::tsv qw( read_tsv );

use Getopt::Long;
use JSON;

my %opt;
GetOptions( \%opt, 'fusions=s', 'id=s', 'clarity-sample-id=s', 'clarity-pool-id=s', 'group=s', 'qc=s', 'classification=s', 'expr=s' );


my( $fus_file, $name ) = ( $opt{fusions}, $opt{id} );
my @groups = split /,/, $opt{group};


# Read QC data
my @QC;
my $vcf;
my $samples;
if( $opt{qc} ) {
    my @qc_files = split /,/, $opt{qc};
    foreach( @qc_files ) {
	if( -s $_ ) {
	    push @QC, read_json($_)
	}
	else {
	    print STDERR "WARNING: QC-json does not exist: $_\n.";
	}
    }
}

my $classification_data;
if( $opt{classification} ) {
    $classification_data = read_json($opt{classification});
}

my $expression_data;
if( $opt{expr} ) {
    $expression_data = read_json($opt{expr});
}



#################
# INSERT SAMPLE #
#################

# Connect to mongodb
my $client = MongoDB->connect();


# Prepare data to insert into sample collection
my $SAMP_COLL = $client->ns("coyote.samples");

# Find the sample, if it already exists
my $sample_results = $SAMP_COLL->find( {'name'=>$name } );

my $SAMPLE_ID;
my $cnt = 0;
while (my $s = $sample_results->next) {
    $SAMPLE_ID = $s->{'_id'}->{value};
    $cnt++;
}
warn "Sample with name $name already exists.\n" if $cnt > 0;


# Insert sample document if it doesn't already exist
if( $cnt == 0 ) {
    my %sample_data = ( 'name'=>$name, 'groups'=>\@groups, 'time_added'=>DateTime->now, 'fusion_files'=>[$fus_file] );
    if ( scalar @QC > 0 ) {
	$sample_data{QC} = \@QC;
    }
    
    # Add clarity information if specified
    if( $opt{'clarity-sample-id'} ) {
	$sample_data{'clarity-sample-id'} =  $opt{'clarity-sample-id'};
    }
    if( $opt{'clarity-pool-id'} ) {
	$sample_data{'clarity-pool-id'} =  $opt{'clarity-pool-id'};
    }
    
    # Insert into collection
    my $result2 = $SAMP_COLL->insert_one(\%sample_data); 
    $SAMPLE_ID = $result2->inserted_id->value;
}


# Add expression classification data to samples document
if( $classification_data ) {
    my $class_result = $SAMP_COLL->update_one( {'_id' => MongoDB::OID->new(value => $SAMPLE_ID)}, {'$set' => {'classification' => $classification_data}} );
}

# Add expression classification data to samples document
if( $expression_data ) {
    my $expr_result = $SAMP_COLL->update_one( {'_id' => MongoDB::OID->new(value => $SAMPLE_ID)}, {'$set' => {'expr' => $expression_data}} );
}


# Add fusions to mongodb collection 
if( $fus_file and -s $fus_file ) {
    my $fusions = read_json($fus_file);
    foreach my $f ( @$fusions ) {
	$f->{SAMPLE_ID} = $SAMPLE_ID;
    }

    my $fusions_coll = $client->ns("coyote.fusions");
    my $fus = $fusions_coll->with_codec( prefer_numeric => 1 );
    my $result = $fus->insert_many($fusions);
}


sub fix {
    my $str = shift;
    #$str =~ s/-/_/g;
    return $str;
}




sub read_json {
    my $fn = shift;

    print STDERR "Reading json $fn\n";

    open( JSON, $fn );
    my @json = <JSON>;
    my $decoded = decode_json( join("", @json ) );
    close JSON;

    return $decoded;
}
