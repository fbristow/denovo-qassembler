#!/usr/bin/perl

use Modern::Perl qw(2010);

use File::Find qw(find);
use File::Basename qw(fileparse);
use Getopt::Long qw(GetOptions);
use Term::ANSIColor qw(:constants);
use Parallel::ForkManager;

my $base_dir = "./experiments";
my ($num_reads, $read_length, $std_dev) = (40000, 450, 40);
my $processors = 1;
my $read_cmd = "/home/fbristow/Projects/qassembler/readsim/readsim.pl";
my $experiment_file;
my $distribution = "squares";
my @experiments;

say BOLD, BLUE, "Info: ", RESET, "Generating reads for experiments.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"read-length=i" => \$read_length,
	"std-dev=i" => \$std_dev,
	"experiment-file=s" => \$experiment_file,
	"processors=i" => \$processors,
	"read-cmd=s" => \$read_cmd,
	"distribution=s" => \$distribution
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
	my ($num_variants) = $experiment =~ /(\d+)\.fna$/;
	my ($filename, $directory, undef) = fileparse($experiment, '.fna');

	# generate an abundance file in the appropriate directory
	my @abundances;

	if ($distribution eq 'squares') {
		@abundances = map { $_ ** 2 } (1..$num_variants);
	} elsif ($distribution eq 'cubes') {
		@abundances = map { $_ ** 3 } (1..$num_variants);
	} elsif ($distribution eq 'linear') {
		@abundances = (1..$num_variants);
	} elsif ($distribution eq 'nlogn') {
		@abundances = map { $_ * log($_) } (1..$num_variants);
	}
	open (my $abundance_fh, '>', "$directory/$filename.abundance");
	print $abundance_fh join "\n", @abundances;
	close ($abundance_fh);

	my $cmd = "$^X $read_cmd --reads $experiment --lengths $read_length --stdev $std_dev --numreads $num_reads ";
	$cmd .= "--no-errors --no-mates --output $directory/$filename-reads.fna --abundfile $directory/$filename.abundance ";
	$cmd .= "2>&1 > /dev/null";

	qx($cmd);

	$pm->finish();
}

$pm->wait_all_children();
