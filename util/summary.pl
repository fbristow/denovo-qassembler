#!/usr/bin/perl

use Modern::Perl '2010';
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use File::Basename qw(fileparse);
use Text::CSV;
use Text::ASCIITable;
use Term::ANSIColor qw(:constants);

my $base_dir = "./";
my $report_location = "summary";
my $output_format = "txt";
my $count = 0;
my $metric;

say BOLD, BLUE, "Info: ", RESET, "Summarizing reports.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"report-location=s" => \$report_location,
	"output-format=s" => \$output_format,
	"metric=s" => \$metric,
);

my @assignments;
my $find_assignments = sub {
	push @assignments, $File::Find::name if ($_ =~ /assignment$/);
};
find($find_assignments, $base_dir);

die BOLD, RED, "Error: ", RESET, "No assignments could be found, run assign.pl first." unless @assignments;

my %reports;
foreach my $assign_file (@assignments) {
	$count++;
	say "Working on assignment ($assign_file) $count/".scalar(@assignments);
	my ($name, $path, $suffix) = fileparse($assign_file, '.assignment');
	my ($length, $contigs) = $path =~ /(\d+)-(\d+)-reads.fna-contigs/;

	my $matches = 0;
	my $avg_match_score = 0;
	my $unmatched_refs = 0;
	my $unmatched_contigs = 0;
	my $perfect_contigs = 0;

	my $csv = Text::CSV->new();
	open (my $assignment_in, '<', $assign_file);
	while (my $assignment_ref = $csv->getline($assignment_in)) {
		my @assignment = @{$assignment_ref};
		if ($assignment[1] eq '') {
			$unmatched_contigs++;
		} elsif ($assignment[0] eq '') {
			$unmatched_refs++;
		} else {
			$matches++;
			$avg_match_score += $assignment[2];

			if ($assignment[2] == 100) {
				$perfect_contigs++;
			}
		}
	}

	$reports{$length}{$contigs} = {
		'matches' => $matches,
		'average_match_score' => $avg_match_score / $matches,
		'unmatched_refs' => $unmatched_refs,
		'unmatched_contigs' => $unmatched_contigs,
		'perfect_contigs' => $perfect_contigs,
	};
}

my @headers = ("", sort keys %reports);
my %row_ids;
my @rows;

foreach my $length (sort keys %reports) {
	foreach my $variants (sort keys %{$reports{$length}}) {
		$row_ids{$variants}++;
	}
}

foreach my $variants (sort {$a <=> $b} keys %row_ids) {
	my @row;
	push @row, $variants;
	foreach my $length (sort keys %reports) {
		my $report = $reports{$length}{$variants};
		my $summary;
		if ($metric) {
			$summary = $report->{$metric};
		} else {
			$summary = "Matches: ".$report->{'matches'}."\n";
			$summary .= "Average Score: ".$report->{'average_match_score'}."\n";
			$summary .= "Unmatched References: ".$report->{'unmatched_refs'}."\n";
			$summary .= "Unmatched Contigs: ".$report->{'unmatched_contigs'}."\n";
			$summary .= "100% matches: ".$report->{'perfect_contigs'}."\n";
		}
		push @row, $summary;
	}
	push @rows, \@row;
}

if ($output_format eq 'txt') {
	my $reporter = Text::ASCIITable->new({ headingText => "Summary Report"});
	$reporter->setCols(@headers);
	$reporter->addRowLine();
	foreach my $row (@rows) {
		$reporter->addRow(@{$row});
		$reporter->addRowLine();
	}
	open (my $summary_report, '>', "$report_location.$output_format");
	print $summary_report $reporter;
	close $summary_report;
} elsif ($output_format eq 'csv') {
	my $csv = Text::CSV->new( { binary => 1 } );
	unshift @rows, \@headers;
	$csv->eol("\r\n");
	open (my $fh, '>', "$report_location.$output_format");
	$csv->print($fh, $_) for @rows;
	close $fh;
} else {
	die "Output format not recognized.";
}
