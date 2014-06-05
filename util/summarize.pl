#!/usr/bin/perl

use Modern::Perl qw(2010);
use Text::CSV;
use Term::ANSIColor qw(:constants);
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use List::Util qw(sum);

my $base_directory = "./";
my $output_format = "csv";

say BOLD, BLUE, "Info: " , RESET, "Summarizing summaries.";

GetOptions(
	"base-dir=s" => \$base_directory,
	"output-format=s" => \$output_format,
);

my @summaries;
my $find_summaries = sub {
	push @summaries, $File::Find::name if ($_ =~ /csv$/);
};
find($find_summaries, $base_directory);

die BOLD, RED, "Error: ", RESET, "No summaries could be found, run summary.pl first." unless @summaries;

my %reports;
foreach my $summary (@summaries) {
	my @headers;
	open (my $csv_in, '<', $summary);
	my $csv = Text::CSV->new();
	my $header_line = $csv->getline($csv_in);
	foreach my $header (@{$header_line}) {
		push @headers, $header unless $header eq '';
	}
	while (my $row = $csv->getline($csv_in)) {
		my @values = @{$row};
		my $variants = shift @values;
		for (my $i = 0; $i < @values; $i++) {
			my $header = $headers[$i];
			push @{$reports{$header}{$variants}}, $values[$i];
		}
	}
}

foreach my $length (keys %reports) {
	foreach my $variants (keys %{$reports{$length}}) {
		my $sum = sum(@{$reports{$length}{$variants}});
		my $avg = $sum / scalar(@{$reports{$length}{$variants}});
		$reports{$length}{$variants} = $avg;
	}
}

my @headers = ("", sort {$a <=> $b} keys %reports);
my %label_count;

foreach my $header (@headers) {
	next if $header eq '';
	foreach my $row_label (keys %{$reports{$header}}) {
		$label_count{$row_label}++;
	}
}
my @row_labels = keys %label_count;

my @rows;
push @rows, \@headers;
foreach my $row_label (sort {$a <=> $b } @row_labels) {
	my @row = ($row_label);
	foreach my $column (sort {$a <=> $b} keys %reports) {
		push @row, $reports{$column}{$row_label};
	}
	push @rows, \@row;
}

my $csv = Text::CSV->new();
$csv->eol("\r\n");
open (my $csv_out, '>', "summary.csv");
$csv->print($csv_out, $_) for @rows;
close $csv_out;
