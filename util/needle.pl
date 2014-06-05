#!/usr/bin/perl

use Modern::Perl qw(2010);

use Bio::SeqIO;
use Bio::Factory::EMBOSS;

use Parallel::ForkManager;

use Getopt::Long qw(GetOptions);

use File::Find qw(find);
use File::Basename qw(fileparse);
use File::Temp;

use Term::ANSIColor qw(:constants);

my $base_dir = "./";
my $processors = 1;
my $count = 0;
my $experiment_file;
my @experiments;

say BOLD, BLUE, "Info: ", RESET, "Aligning contigs to references.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"processors=i" => \$processors,
	"experiment-file=s" => \$experiment_file,
);

if (defined $experiment_file) {
	push @experiments, $experiment_file;
} else {
	my $find_contigs = sub {
		push @experiments, $File::Find::name if ($_ =~ /^\d+\-filtered.fna$/);
	};
	find ($find_contigs, $base_dir);
}

my $pm = Parallel::ForkManager->new($processors);
my $emboss_factory = Bio::Factory::EMBOSS->new();
my $needle = $emboss_factory->program('needle');

foreach my $sequence_fna (@experiments) {
	$count++;
	say BOLD, BLUE, "Info: ", RESET, "Aligning experiment [$sequence_fna]."; 
	$pm->start() and next;
	my ($name, $path, $suffix) = fileparse($sequence_fna, '.fna');
	my $seq_in = Bio::SeqIO->new(-file=>$sequence_fna);
	while (my $seq = $seq_in->next_seq()) {
		my ($source) = $path =~ /.*\/(.*)-reads\.fna-contigs/;
		my ($id) = $seq->primary_id =~ /(\d+)\(/;

		my $forward_out = "$path/$name-$id-forward.needle";
		my $reverse_out = "$path/$name-$id-reverse.needle";

		# now run needle on both:
		my $output = $needle->run({
				-asequence => $seq,
				-bsequence => "$path/../$source.fna",
				-gapopen => '10.0',
				-gapextend => '0.5',
				-outfile => $forward_out
			}) unless -e $forward_out;
		say BOLD, BLUE, "Debug: ", RESET, "Output of needle: $output.";
		$output = $needle->run({
				-asequence => $seq->revcom,
				-bsequence => "$path/../$source.fna",
				-gapopen => '10.0',
				-gapextend => '0.5',
				-outfile => $reverse_out
			}) unless -e $reverse_out;
		say BOLD, BLUE, "Debug: ", RESET, "Output of needle: $output.";
	}

	$pm->finish();
}

$pm->wait_all_children();


