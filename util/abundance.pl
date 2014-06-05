#!/usr/bin/perl

use Modern::Perl qw(2010);
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use File::Basename qw(fileparse);
use Term::ANSIColor qw(:constants);
use Text::CSV;
use Bio::SeqIO;
use List::Util qw(sum);

my $base_dir = "./";
my $method = "residual";

say BOLD, BLUE, "Info: ", RESET, "Computing residual sum of squares for abundance estimation.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"method=s" => \$method,
);

my @experiments;
my $find_experiments = sub {
	push @experiments, $File::Find::name if (-d $_ and $_ =~ /-contigs$/);
};
find ($find_experiments, $base_dir);

foreach my $experiment (@experiments) {
	my ($filename, $directory, undef) = fileparse($experiment, '.fna');
	# Load the filtered contigs into memory:
	my @contig_files;
	my $find_contigs = sub {
		push @contig_files, $File::Find::name if ($_ =~ /^\d+-filtered.fna$/);
	};
	find ($find_contigs, $experiment);

	my %contigs;
	my $abundance_sum = 0;
	foreach my $contig_file (@contig_files) {
		my $seq_in = Bio::SeqIO->new(-file => $contig_file);
		while (my $seq = $seq_in->next_seq()) {
			my ($abundance) = $seq->description =~ /markov-chain: ([\d\.]+)/;
			my $id = $seq->primary_id();
			$contigs{$id} = $abundance;
			$abundance_sum += $abundance;
		}
	}

	# Compute relative abundance for the contigs:
	foreach my $contig (keys %contigs) {
		$contigs{$contig} = $contigs{$contig} / $abundance_sum;
	}

	# Load the abundance file into memory:
	my $abundance_file;
	my $find_abundance = sub {
		$abundance_file = $File::Find::name if ($_ =~ /abundance$/);
	};
	find ($find_abundance, $directory);

	my @reference_abundances;
	open (my $abundance_in, '<', $abundance_file);
	push @reference_abundances, $_ while (<$abundance_in>);
	close ($abundance_in);

	$abundance_sum = sum (@reference_abundances);
	my @reference_relative_abundances = map { $_ / $abundance_sum } @reference_abundances;

	my $reference_file;
	my $find_reference = sub {
		$reference_file = $File::Find::name if ($_ =~ /\d+-\d+\.fna$/);
	};
	find ($find_reference, $directory);

	my %references;
	my $i = 0;
	my $seq_in = Bio::SeqIO->new(-file => $reference_file);
	while (my $seq = $seq_in->next_seq()) {
		my ($reference_id) = $seq->primary_id() =~ /(.*?)_?formatted/;
		$references{$reference_id} = $reference_relative_abundances[$i];
		$i++;
	}

	# Find the mappings file:
	my $mapping_file;
	my $find_mapping = sub {
		$mapping_file = $File::Find::name if ($_ =~ /assignment$/);
	};
	find ($find_mapping, $directory);
	my $residual_sum = 0;
	my $reference_max = -1;
	open (my $mapping_in, '<', $mapping_file);
	my $csv = Text::CSV->new();
	while (my $row = $csv->getline($mapping_in)) {
		my ($contig) = $row->[0] =~ /(\d+\(\d+bp\))/;
		my ($reference) = $row->[1] =~ /(.*?)\_/;
		next unless defined $contig and defined $reference;

		if ($method eq 'residual') {
			my $difference = ($references{$reference}) - ($contigs{$contig});
			$residual_sum += ($difference * $difference);
		} elsif ($method eq 'absolute') {
			my $difference = abs(($references{$reference}) - ($contigs{$contig}));
			$residual_sum += ($difference / $references{$reference});
		} elsif ($method eq 'max') {
			if ($references{$reference} > $reference_max) {
				$reference_max = $references{$reference};
				$residual_sum = $references{$reference} - $contigs{$contig};
			}
		}
	}
	close ($mapping_in);

	my ($residual_filename, $residual_directory, undef) = fileparse($abundance_file, '.abundance');
	open (my $residual_out, '>', "$residual_directory/$residual_filename.residual");
	say $residual_out $residual_sum;
	close $residual_out;
}
