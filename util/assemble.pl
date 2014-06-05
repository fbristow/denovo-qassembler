#!/usr/bin/perl

use Modern::Perl qw(2010);
use File::Find qw(find);
use File::Basename qw(fileparse);
use Getopt::Long qw(GetOptions);
use Term::ANSIColor qw(:constants);
use Parallel::ForkManager;

my $base_dir = "./experiments";
my $processors = 1;
my $kmer_length = 327;
my $assembler_cmd = "denovo-qassembler";
my $assembler_cmd_flags = "";
my $experiment_file;
my @experiments;

say BOLD, BLUE, "Info: ", RESET, "Assembling reads.";

GetOptions (
	"base-dir=s" => \$base_dir,
	"experiment-file=s" => \$experiment_file,
	"processors=i" => \$processors,
	"kmer-length=i" => \$kmer_length,
	"assembler-cmd-flags=s" => \$assembler_cmd_flags,
	"assembler-cmd=s" => \$assembler_cmd,
);

if (defined $experiment_file) {
	push @experiments, $experiment_file;
} else {
	my $find_banks = sub {
		push @experiments, $File::Find::name if ($_ =~ /-reads\.fq$/);
	};
	find ($find_banks, $base_dir);
}

my $pm = Parallel::ForkManager->new($processors);
foreach my $experiment (@experiments) {
	say BOLD, BLUE, "Info: ", RESET, "Assembling experiment [$experiment].";
	$pm->start() and next;

	my ($filename, $directory, undef) = fileparse($experiment, '.fq');
	my $cmd = "$assembler_cmd_flags $assembler_cmd -i $experiment -k $kmer_length ";
	$cmd .= "-s --sequence-dir=$directory/$filename-contigs/ -g --graph-dir=$directory/$filename-graphs/ ";
	$cmd .= "--abundance-method markov-chain 2> $directory/denovo-qassembler.stderr > $directory/denovo-qassembler.stdout";
	qx($cmd);

	$pm->finish();
}

$pm->wait_all_children();
