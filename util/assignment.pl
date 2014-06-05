#!/usr/bin/perl

use Modern::Perl '2010';
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use Text::CSV;
use List::Util qw(sum);
use Algorithm::Munkres qw(assign);
use Term::ANSIColor qw(:constants);

my $base_dir = "./";
my $count = 0;

say BOLD, BLUE, "Info: ", RESET, "Assigning contigs to references.";

GetOptions(
	"base-dir=s" => \$base_dir,
);

my @reports;
my $find_reports = sub {
	push @reports, $File::Find::name if ($_ =~ /report.csv$/);
};
find ($find_reports, $base_dir);

die BOLD, RED, "Error: ", RESET, "No reports could be found, run report.pl first." unless @reports;

foreach my $report (@reports) {
	$count++;
	say "Working on report $count/".scalar(@reports);
	my @score_matrix;

	my $csv = Text::CSV->new();
	open (my $csv_in, '<', $report);
	my @references = @{$csv->getline($csv_in)};
	shift @references; # discard the empty element in the header row

	my @contigs;

	# get the remaining rows in pairs (forward and reverse)
	while (my ($forward_ref, $reverse_ref) = ($csv->getline($csv_in), $csv->getline($csv_in))) {
		my @forward = @$forward_ref;
		my @reverse = @$reverse_ref;

		# get rid of the row headers
		my $forward_name = shift @forward;
		my $reverse_name = shift @reverse;

		my $forward_avg = sum(@forward) / scalar(@forward);
		my $reverse_avg = sum(@reverse) / scalar(@reverse);

		if ($forward_avg > $reverse_avg) {
			push @contigs, $forward_name;
			@forward = map { 100 - $_ } @forward;
			push @score_matrix, \@forward;
		} else {
			push @contigs, $reverse_name;
			@reverse = map { 100 - $_ } @reverse;
			push @score_matrix, \@reverse;
		}
	}
	my @assignment;
	assign (\@score_matrix, \@assignment);

	my $csv_reporter = Text::CSV->new();
	open (my $assignment_out, '>', "$report.assignment");
	for (my $i = 0; $i < scalar(@assignment); $i++) {
		my @assign_report;
		if ($assignment[$i] < scalar(@references) and defined $score_matrix[$i]) {
			my $score = 100 - ($score_matrix[$i][$assignment[$i]]);
			@assign_report = ($contigs[$i], $references[$assignment[$i]], $score);
		} elsif (defined $contigs[$i]) {
			@assign_report = ($contigs[$i], "", "");
		} else {
			@assign_report = ("", $references[$assignment[$i]], "");
		}
		$csv_reporter->combine(@assign_report);
		say $assignment_out $csv_reporter->string();
	}
	close($assignment_out);
	
}
