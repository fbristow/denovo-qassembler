#!/usr/bin/perl

use Modern::Perl qw(2010);
use Text::CSV;
use Term::ANSIColor qw(:constants);
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use File::Basename qw(fileparse);

my $base_directory = "./";

say BOLD, BLUE, "Info: ", RESET, "Summarizing residual sum of squares.";

GetOptions (
	"base-directory=s" => \$base_directory,
);

my @residual_files;
my $find_residuals = sub {
	push @residual_files, $File::Find::name if ($_ =~ /residual$/);
};
find ($find_residuals, $base_directory);

my %residuals;
foreach my $residual_file (@residual_files) {
	my ($filename, $directory, undef) = fileparse($residual_file, '.residual');
	my ($length, $variants) = $filename =~ /(\d+)-(\d+)$/;
	open (my $residual_in, '<', $residual_file);
	$residuals{$length}{$variants} += <$residual_in>;
	close ($residual_in);
}

my @headers = ("", sort {$a <=> $b} keys %residuals);
my %label_count;
foreach my $header (@headers) {
	next if $header eq '';
	foreach my $row_label (keys %{$residuals{$header}}) {
		$label_count{$row_label}++;
	}
}

my @row_labels = keys %label_count;

my @rows;
push @rows, \@headers;
foreach my $row_label (@row_labels) {
	my @row = ($row_label);
	foreach my $column (sort {$a <=> $b} keys %residuals) {
		push @row, $residuals{$column}{$row_label} / 50;
	}
	push @rows, \@row;
}

my $csv = Text::CSV->new();
$csv->eol("\r\n");
open (my $csv_out, '>', "residuals.csv");
$csv->print($csv_out, $_) for @rows;
close $csv_out;
