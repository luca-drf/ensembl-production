=pod

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::Production::Pipeline::EBeye::DumpTypeFactory

=head1 DESCRIPTION

Small extension of the job factory to do default database type submission.

Allowed parameters are:

=over 8

=item types - The database types to use; defaults to core and vega

=back

=cut

package Bio::EnsEMBL::Production::Pipeline::EBeye::DumpTypeFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Hive::RunnableDB::JobFactory/;
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;

sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults()},
    column_names => ['type','species'],
    default_types => [qw/core vega/],
  };
}

sub fetch_input {
  my ($self) = @_;
  my $user_types = $self->param('types');
  my $types = (defined $user_types && @{$user_types}) ? $user_types : $self->param('default_types');
  $types = wrap_array($types);
  my @inputlist = map { [ $_, $self->param('species') ] } @{$types};
  $self->param('inputlist', \@inputlist);
  return;
}

1;
