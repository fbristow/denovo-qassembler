#!/usr/bin/perl

use Modern::Perl qw(2010);

use File::Find qw(find);
use File::Basename qw(fileparse);
use Getopt::Long qw(GetOptions);
use Term::ANSIColor qw(:constants);
use Parallel::ForkManager;
use Bio::SeqIO;
use List::Util qw(sum);
use POSIX qw(ceil);

my $base_dir = "./experiments";
my $total_coverage = 1500;
my $processors = 1;
my $read_cmd = "/home/fbristow/Applications/ART/art_454";
my $experiment_file;
my @experiments;

say BOLD, BLUE, "Info: ", RESET, "Generating reads for experiments.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"experiment-file=s" => \$experiment_file,
	"processors=i" => \$processors,
	"read-cmd=s" => \$read_cmd,
	"total-coverage=i" => \$total_coverage,
);

if (defined $experiment_file) {
	push @experiments, $experiment_file;
} else {
	my $find_slices = sub {
		push @experiments, $File::Find::name if ($_ =~ /\.fna$/);
	};
	find ($find_slices, $base_dir);
}

my $pm = Parallel::ForkManager->new($processors);
foreach my $experiment (@experiments) {
	say BOLD, BLUE, "Info: ", RESET, "Generating reads for experiment [$experiment].";
	$pm->start() and next;
	my $seq_in = Bio::SeqIO->new(-file => $experiment);
	my ($filename, $directory, undef) = fileparse($experiment, '.fna');
	my ($num_variants) = $experiment =~ /(\d+)\.fna$/;

	# generate an abundance file in the appropriate directory
	my @abundances = map { $_ ** 2} (1..$num_variants);
	my @output_filenames;
	my $total_abundance = sum(@abundances);
	my $variant = 0;

	while (my $seq = $seq_in->next_seq()) {
		my $relative_abundance = $abundances[$variant++] / $total_abundance;
		my $relative_coverage = ceil($total_coverage * $relative_abundance);
		my $variant_filename = "$directory/$filename-$variant";

		# write the single variant out into its own file:
		my $seq_out = Bio::SeqIO->new(-file => ">$variant_filename.fna", -format => 'fasta');
		$seq_out->write_seq($seq);

		my $cmd = "$read_cmd -s $variant_filename.fna $variant_filename-reads-part $relative_coverage 2>&1 > /dev/null";
		push @output_filenames, "$variant_filename-reads-part.fq";
		qx($cmd);
	}

	# now cat all the files together:
	local $/;
	open (my $cat, '>', "$directory/$filename-reads.fq");
	for my $reads_file (@output_filenames) {
		open (my $reads, '<', $reads_file);
		while(<$reads>) {
			print $cat $_;
		}
		close $reads;
	}
	close $cat;

	$pm->finish();
}

$pm->wait_all_children();
