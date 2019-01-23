=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

    RDFDump - Hive Process to start with a species name and produce triples

=head1 DESCRIPTION

    This process is now based on the ensembl-io unified framework.
    It gets the appropriate data from BulkFetcher, and serialises to RDF.

=cut

package Bio::EnsEMBL::Production::Pipeline::RDF::RDFDump;

use strict;
use warnings;
use FindBin qw($Bin);

use parent ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

use IO::File;
use File::Spec::Functions qw/catdir/;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;

use Bio::EnsEMBL::Production::DBSQL::BulkFetcher;

use Bio::EnsEMBL::IO::Object::RDF;
use Bio::EnsEMBL::IO::Translator::Slice;
use Bio::EnsEMBL::IO::Translator::BulkFetcherFeature;
use Bio::EnsEMBL::IO::Writer::RDF;
use Bio::EnsEMBL::IO::Writer::RDF::XRefs;

sub fetch_input {
    my $self = shift;
    $self->param_required('species');   # just make sure it has been passed
    $self->param_required('config_file');
    $self->param_required('release');
    my $eg = $self->param('eg');
    $self->param('eg', $eg);

    if($eg){
      my $base_path  = $self->build_base_directory();
      $self->param('base_path', $base_path);
      my $release = $self->param('eg_version');
      $self->param('release', $release);
    }

}


sub run {
  my $self = shift;

  my $species = $self->param('species');
  my $release = $self->param('release');

  ### Fetch the data, i.e. genes-transcripts-translations and sequence regions
  #
  # configure bulk extractor to go all the way down to protein features.
  # can also be told to stop at transcript level as well as others.
  # my $bulk = Bio::EnsEMBL::Production::DBSQL::BulkFetcher->new(-level => 'protein_feature');
  my $dba = $self->get_DBAdaptor;
  my $formatter = Bio::EnsEMBL::Production::Search::RDFFormatter->new(
      -config_file  => $self->param_required('config_file'),
      -ontology_dba => $self->get_DBAdaptor('ontology')
  );
  my $sub_dir = $self->get_data_path('json');
  my $genome_in_file = File::Spec->catfile($sub_dir, $species . '_genome.json');
  my $genes_in_file = File::Spec->catfile($sub_dir, $species . '_genes.json');

  my $hive_dbc = $self->dbc;
  $hive_dbc->disconnect_if_idle() if defined $hive_dbc;

  my $slices = $self->get_Slices(undef, ($species eq 'homo_sapiens')?1:0);
  #
  ############

  ### Dump RDF
  #
  my $path = $self->data_path();
  $path = $self->get_dir($release) unless defined $path && $path ne '';

  # my $core_rdf_file = catdir($path, $species . ".ttl");
  my $core_rdf_file = File::Spec->catfile($path, $species . ".ttl");
  my $core_dba = $self->get_DBAdaptor();
  $formatter->reformat_core($genome_in_file, $genes_in_file, $core_dba, $core_rdf_file);

  # xrefs are optional
  # my $xrefs_rdf_file = catdir();
  my $xrefs_rdf_file = File::Spec->catfile($path, $species . "_xrefs.ttl");
  $formatter->reformat_xrefs($genome_in_file, $genes_in_file, $core_dba, $xrefs_rdf_file) if $self->param('xref');
  # $self->dump_xrefs_rdf($xrefs_rdf_file, $genes) if $self->param('xref');

  #
  ###########

  ### Add a graph file for Virtuoso loading.
  my $graph_path = $path;
  $self->param('dir', $graph_path);
  $graph_path = $self->get_dir($release) unless $graph_path;

  # graph files need to be named exactly the same as the underlying
  # data for Virtuoso to correctly namespace the data
  $self->create_virtuoso_file(sprintf("%s/%s.ttl.gz.graph", $graph_path, $self->production_name));
  $self->create_virtuoso_file(sprintf("%s/%s_xrefs.ttl.gz.graph", $graph_path, $self->production_name));
  #
  ###########

  ### Compress the files
  system("gzip -n $core_rdf_file");
  system("gzip -n $xrefs_rdf_file");

  ### Create list of files to validate
  my @files_to_validate = ($core_rdf_file . '.gz');
  push @files_to_validate, $xrefs_rdf_file . '.gz' if $self->param('xref');
  $self->param('validate_me', \@files_to_validate);

  $dba->dbc()->disconnect_if_idle();
}

sub write_output {  # store and dataflow
    my $self = shift;
    my $files = $self->param('validate_me');
    my $dir = $self->param('dir');
    while (my $file = shift @$files) {
        $self->dataflow_output_id({filename => $file},2);
    }
    $self->dataflow_output_id({dir => $dir},3);
}

sub data_path {
  my ($self) = @_;

  $self->throw("No 'species' parameter specified")
    unless $self->param('species');

  return $self->get_dir('rdf', $self->param('species'));
}


sub create_virtuoso_file {
  my $self = shift;
  my $path = shift; # a .graph file, named after the rdf file.

  my $version_graph_uri = sprintf "http://rdf.ebi.ac.uk/dataset/%s/%d", $self->param('eg')?'ensemblgenomes':'ensembl', $self->param('release');
  my $graph_uri = sprintf "%s/%s", $version_graph_uri, $self->production_name;

  my $graph_fh = IO::File->new($path, 'w') or die "Cannot open virtuoso file $path: $!\n";
  print $graph_fh $graph_uri . "\n";
  $graph_fh->close();

}

1;
