#!/usr/bin/perl

use Modern::Perl '2012';

use Hash::MultiValue;

use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use File::Basename qw(fileparse basename);
use Term::ANSIColor qw(:constants);

use Bio::SeqIO;
use Parallel::ForkManager;

my $base_dir = "./";
my $min_length_pct = 0;
my $processors = 1;

say BOLD, BLUE, "Info: ", RESET, "Filtering reads.";

GetOptions (
	"base-dir=s" => \$base_dir,
	"processors=i" => \$processors,
	"min-length-pct=f" => \$min_length_pct,
);

my @sequence_files;
my $wanted = sub {
	push @sequence_files, $File::Find::name if $_ =~ /^\d+\.fna/;
};
find ($wanted, $base_dir);

my $sorted_files = Hash::MultiValue->new();

foreach my $file (@sequence_files) {
	my ($name, $path, $suffix) = fileparse ($file);
	$sorted_files->add($path => $name);
}

my $pm = Parallel::ForkManager->new($processors);

my $count = 0;
my $total = scalar(keys %$sorted_files);
foreach my $directory (keys %$sorted_files) {
	$count++;
	say BOLD, BLUE, "Info: ", RESET, "Working on directory: [$directory] ($count of $total)";
	$pm->start() and next;
	my ($min_length) = $directory =~ /(\d+)-\d+-reads-contigs/;
	$min_length = ($min_length // 0) * $min_length_pct;

	# 1) load all sequences in the directory into memory:
	my $sequences = Hash::MultiValue->new();
	foreach my $file ($sorted_files->get_all($directory)) {
		my $seq_in = Bio::SeqIO->new(-file => "$directory/$file");
		while (my $seq = $seq_in->next_seq()) {
			if ($seq->length() > $min_length) {
				$sequences->add($file => $seq);
			}
		}
	}

	my @values = $sequences->values();
	my @keys = $sequences->keys();
	my %delete;

	OUTER:	for (my $i = 0; $i < scalar(@values); $i++) {
		# don't bother checking a sequence if we've already chucked it out:
		next if $delete{$i};
		my $seq = $values[$i];
		for (my $j = 0; $j < scalar(@values); $j++) {
			my $cmp = $values[$j];
			# don't compare a sequence to itself:
			next if $i == $j;
			# don't bother checking a sequence if we've already chucked it out:
			next if $delete{$j};
			if ($seq->length() > $cmp->length()) {
				if (index ($seq->seq, $cmp->seq) != -1 or index($seq->revcom->seq, $cmp->seq) != -1) {
					$delete{$j}++;
				}
			} else {
				if (index($cmp->seq, $seq->seq) != -1 or index($cmp->revcom->seq, $seq->seq) != -1) {
					$delete{$i}++;
					next OUTER;
				}
			}
		}
	}

	say BOLD, BLUE, "Info: ", RESET, "Found ".scalar(keys %delete)." duplicates.";

	# remove duplicates
	foreach my $el (keys %delete) {
		delete $values[$el];
		delete $keys[$el];
	}

	# then print out each file again
	for (my $i = 0; $i < scalar(@values); $i++) {
		next unless defined $keys[$i];
		my $filename = basename ($keys[$i], '.fna');
		my $seq_out = Bio::SeqIO->new(-file => ">>$directory/$filename-filtered.fna", -format => 'fasta');
		$seq_out->write_seq($values[$i]);
	}
	
	$pm->finish();
}

$pm->wait_all_children();
