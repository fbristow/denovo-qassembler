#!/usr/bin/perl

use Modern::Perl qw(2010);

use Getopt::Long qw(GetOptions);
use Text::ASCIITable;
use Hash::MultiValue;
use List::MoreUtils qw(uniq);
use Term::ANSIColor qw(:constants);

use File::Path qw(make_path);
use File::Basename qw(fileparse);

use Bio::AlignIO;
use Bio::SeqIO;

my $alignment_in = "ML2155.aln";
my $base_dir = "./";

my ($variants_min, $variants_max, $variants_step) = (2, 20, 2);
my ($length_min, $length_max, $length_step) = (1000, 10000, 1000);

my $report = 1;

say BOLD, BLUE, "Info: ", RESET, "Generating experiments.";

GetOptions(
	"alignment=s" => \$alignment_in,
	"variants-min=i" => \$variants_min,
	"variants-max=i" => \$variants_max,
	"variants-step=i" => \$variants_step,
	"length-min=i" => \$length_min,
	"length-max=i" => \$length_max,
	"length-step=i" => \$length_step,
	"report" => \$report,
	"base-dir=s" => \$base_dir
);

my $align_io = Bio::AlignIO->new(-file => $alignment_in);
my $alignment = $align_io->next_aln();

say BOLD, BLUE, "Info: ", RESET, "Alignment length: [".$alignment->length()."].";

if ($alignment->length() < $length_max) {
	say BOLD, YELLOW, "Warning: ", RESET, "Maximum slice length larger than alignment length, reducing length-max to ".$alignment->length().".";
	$length_max = $alignment->length();
}

my %experiments;

for (my $i = $length_min; $i <= $length_max; $i += $length_step) {
	for (my $j = $variants_min; $j <= $variants_max; $j += $variants_step) {
		say BOLD, BLUE, "Info: ", RESET, "Selecting [$j] unique variants with length [$i].";
		# select a random slice from the entire alignment
		my $slice_start = int(rand($alignment->length() - $i));
		my $slice = $alignment->slice($slice_start + 1, $slice_start + $i + 1);
		# select variants from the slice until you find a total of $j unique variants
		my @rows = 0..($slice->num_sequences() - 1);
		my @selected_variants;
		my @sequences_observed;
		while (@selected_variants < $j and @rows > 0) {
			my $variant_row = splice(@rows, int(rand(@rows)), 1);
			my $seq = $slice->get_seq_by_pos($variant_row + 1);
			my $gapless_sequence = uc($seq->seq());
			$gapless_sequence =~ s/[\.\-]//g;
			unless (grep { $_ eq $gapless_sequence } @sequences_observed) {
				push @sequences_observed, $gapless_sequence;
				push @selected_variants, $seq;
			}
		}

		unless (uniq(@sequences_observed) == @sequences_observed) {
			die BOLD, RED, "Error: ", RESET, "Selected duplicate sequences. Bailing.";
		}

		$experiments{$i}{$j} = {
			'slice-start' => $slice_start,
			'slice-end' => $slice_start + $i,
			'selected-variants' => \@selected_variants,
		};
	}
}

say BOLD, BLUE, "Info: ", RESET, "Writing experiments to disk.";
foreach my $length (keys %experiments) {
	foreach my $variants (keys %{$experiments{$length}}) {
		my $directory = "$base_dir/$length/$variants";
		make_path($directory);
		my ($filename, undef, undef) = fileparse($alignment_in, '.aln');
		my $seq_out = Bio::SeqIO->new(
			-file => ">$directory/$filename-$length-$variants.fna",
			-format => 'fasta',
		);

		foreach my $selected_variant (@{$experiments{$length}{$variants}->{'selected-variants'}}) {
			my $slice_range = "(".$experiments{$length}{$variants}->{'slice-start'} . " to " . $experiments{$length}{$variants}->{'slice-end'}.")";
			my $gapless_sequence = $selected_variant->seq();
			$gapless_sequence =~ s/[\.\-]//g;
			$selected_variant->seq($gapless_sequence);
			$selected_variant->desc($slice_range);
			$seq_out->write_seq($selected_variant);
		}
	}
}
