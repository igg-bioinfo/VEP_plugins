=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

gnomADpLI - Add gnomAD pLI and other scores to the VEP output 

=head1 SYNOPSIS

  mv gnomADpLI.pm ~/.vep/Plugins
  mv gnomADpLI_values.txt ~/.vep/Plugins
  ./vep -i variants.vcf --plugin gnomADpLI

=head1 DESCRIPTION

  TO BE UPDATED

  A VEP plugin that adds the probabililty of a gene being 
  loss-of-function intolerant (pLI) to the VEP output.
  
  Lek et al. (2016) estimated pLI using the expectation-maximization 
  (EM) algorithm and data from 60,706 individuals from 
  ExAC (http://exac.broadinstitute.org/about). The closer pLI is to 1, 
  the more likely the gene is loss-of-function (LoF) intolerant. 
  
  Note: the pLI was calculated using a representative transcript and
  is reported by gene in the plugin.

  The data for the plugin is provided by Kaitlin Samocha and Daniel MacArthur. 
  See https://www.ncbi.nlm.nih.gov/pubmed/27535533 for a description 
  of the dataset and analysis.

  The ExACpLI_values.txt file is found alongside the plugin in the 
  VEP_plugins GitHub repository. The file contains the fields gene and pLI 
  extracted from the file at 
    
    ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/
      fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt

  To use another values file, add it as a parameter i.e.

     ./vep -i variants.vcf --plugin gnomADpLI,values_file.txt


=cut

package gnomADpLI;

use strict;
use warnings;

use DBI;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  my $file = $self->params->[0];

  if(!$file) {
    my $plugin_dir = $INC{'gnomADpLI.pm'};
    $plugin_dir =~ s/gnomADpLI\.pm//i;
    $file = $plugin_dir.'/gnomADpLI_values.txt';
  }
  
  die("ERROR: gnomADpLI values file $file not found\n") unless $file && -e $file;
  
  open my $fh, "<",  $file;
  my %pli_scores;
  my %syn_z_scores;
  my %mis_z_scores;
  my %lof_z_scores;
  my %oe_mis_scores;
  my %oe_syn_scores;
  my %oe_lof_scores;
  my %oe_syn_upper_scores;
  my %oe_mis_upper_scores;
  my %oe_lof_upper_scores;
 
  while(<$fh>) {
    chomp;
    my ($gene, $oe_mis, $oe_syn, $pli, $oe_lof, $oe_syn_upper, $oe_mis_upper, $oe_lof_upper, $syn_z, $mis_z, $lof_z) = split;
    next if $pli eq 'pLI';
    ###$scores{lc($gene)} = sprintf("%.2f", $score);
    $pli_scores{lc($gene)} = $pli eq 'NA' ? 'NA' : sprintf("%.2f", $pli);
    $syn_z_scores{lc($gene)} = $syn_z;
    $mis_z_scores{lc($gene)} = $mis_z;
    $lof_z_scores{lc($gene)} = $lof_z;
    $oe_mis_scores{lc($gene)} = $oe_mis;
    $oe_syn_scores{lc($gene)} = $oe_syn;
    $oe_lof_scores{lc($gene)} = $oe_lof;
    $oe_syn_upper_scores{lc($gene)} = $oe_syn_upper;
    $oe_mis_upper_scores{lc($gene)} = $oe_mis_upper;
    $oe_lof_upper_scores{lc($gene)} = $oe_lof_upper;


  }
  
  close $fh;
  
  die("ERROR: No scores read from $file\n") unless scalar keys %pli_scores;
  
  $self->{pli} = \%pli_scores;
  $self->{syn_z} = \%syn_z_scores;
  $self->{mis_z} = \%mis_z_scores;
  $self->{lof_z} = \%lof_z_scores;
  $self->{oe_mis} = \%oe_mis_scores;
  $self->{oe_syn} = \%oe_syn_scores;
  $self->{oe_lof} = \%oe_lof_scores;
  $self->{oe_syn_upper} = \%oe_syn_upper_scores;
  $self->{oe_mis_upper} = \%oe_mis_upper_scores;
  $self->{oe_lof_upper} = \%oe_lof_upper_scores;
  
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    gnomADpLI => "gnomAD pLI value for gene",
    gnomADsyn_z => "gnomAD syn_z value for gene",
    gnomADmis_z => "gnomAD mis_z value for gene",
    gnomADlof_z => "gnomAD lof_z value for gene",
    gnomADoe_mis => "gnomAD oe_mis value for gene",
    gnomADoe_syn => "gnomAD oe_syn value for gene",
    gnomADoe_lof => "gnomAD oe_lof value for gene",
    gnomADoe_syn_upper => "gnomAD oe_syn_upper value for gene",
    gnomADoe_mis_upper => "gnomAD oe_mis_upper value for gene",
    gnomADoe_lof_upper => "gnomAD oe_lof_upper value for gene"
  };
}

sub run {
  my $self = shift;
  my $tva = shift;
  
  my $symbol = $tva->transcript->{_gene_symbol} || $tva->transcript->{_gene_hgnc};
  return {} unless $symbol;
  
  my $return = {};
  $return = {
      gnomADpLI => $self->{pli}->{lc($symbol)} ? $self->{pli}->{lc($symbol)} : '',
      gnomADsyn_z => $self->{syn_z}->{lc($symbol)} ? $self->{syn_z}->{lc($symbol)} : '',
      gnomADmis_z => $self->{mis_z}->{lc($symbol)} ? $self->{mis_z}->{lc($symbol)} : '',
      gnomADlof_z => $self->{lof_z}->{lc($symbol)} ? $self->{lof_z}->{lc($symbol)} : '',
      gnomADoe_mis => $self->{oe_mis}->{lc($symbol)} ? $self->{oe_mis}->{lc($symbol)} : '',
      gnomADoe_syn => $self->{oe_syn}->{lc($symbol)} ? $self->{oe_syn}->{lc($symbol)} : '',
      gnomADoe_lof => $self->{oe_lof}->{lc($symbol)} ? $self->{oe_lof}->{lc($symbol)} : '',
      gnomADoe_syn_upper => $self->{oe_syn_upper}->{lc($symbol)} ? $self->{oe_syn_upper}->{lc($symbol)} : '',
      gnomADoe_mis_upper => $self->{oe_mis_upper}->{lc($symbol)} ? $self->{oe_mis_upper}->{lc($symbol)} : '',
      gnomADoe_lof_upper => $self->{oe_lof_upper}->{lc($symbol)} ? $self->{oe_lof_upper}->{lc($symbol)} : ''
  };
  #return $self->{scores}->{lc($symbol)} ? { ExACpLI_mis_z => $self->{scores}->{lc($symbol)}} : {};
  return $return;
}

1;

