#!/usr/bin/perl

use Modern::Perl '2010';
use Bio::AlignIO;
use Parallel::ForkManager;
use Getopt::Long qw(GetOptions);
use File::Find qw(find);
use IO::Scalar;
use Text::ASCIITable;
use Text::CSV;
use Term::ANSIColor qw(:constants);

my $base_dir = "./";
my $processors = 1;
my $out_fmt = "txt";
my $count = 0;

say BOLD, BLUE, "Info: ", RESET, "Generating reports from needle output.";

GetOptions(
	"base-dir=s" => \$base_dir,
	"processors=i" => \$processors,
	"output-format=s" => \$out_fmt);

my @contig_dirs;

my $find_contigs = sub {
	push @contig_dirs, $File::Find::name if (-d $_ and $_ =~ /contigs$/);
};
find ($find_contigs, $base_dir);

my $pm = Parallel::ForkManager->new($processors);
foreach my $contig_dir (@contig_dirs) {
	$count++;
	$pm->start() and next;
	say BOLD, BLUE, "Info: ", RESET, "Working on file $count/".scalar(@contig_dirs);
	my @alignments;
	my $find_alignments = sub {
		push @alignments, $File::Find::name if ($_ =~ /needle$/);
	};
	find ($find_alignments, $contig_dir);
	my %scores;
	foreach my $needle_file (@alignments) {
		my ($direction) = $needle_file =~ /.*-(.*)\.needle$/;
		my $aln_in = Bio::AlignIO->new(
			-format => 'emboss',
			-file => $needle_file
		);
	
		while (my $aln = $aln_in->next_aln()) {
			my @sequences;
			for my $seq ($aln->each_seq()) {
				push @sequences, $seq->display_id();
			}
			$scores{$sequences[0]}{$sequences[1]}{$direction} = $aln->overall_percentage_identity();
		}
	}

	my $report;

	given ($out_fmt) {
		when (/csv/) {
			$report = csv_output(\%scores);
		}
		default {
			$report = ascii_output(\%scores);
		}
	}

	open (my $report_out, '>', "$contig_dir/report.$out_fmt");
	print $report_out $report;
	close $report_out;
	$pm->finish();
}

$pm->wait_all_children();

sub csv_output {
	my ($scores) = @_;
	my $reporter = Text::CSV->new();
	my $output;
	my $newline = "\r\n";

	my %headers;
	foreach my $seq (sort keys %$scores) {
		foreach my $cmp (sort keys %{$scores->{$seq}}) {
			$headers{$cmp}++;
		}
	}
	my @header = ("", (sort keys %headers));

	$reporter->combine(@header);
	$output .= $reporter->string() . $newline;
	foreach my $seq (sort keys %$scores) {
		my @forward;
		my @reverse;

		push @forward, "$seq (forward)";
		push @reverse, "$seq (reverse)";

		foreach my $cmp (sort keys %{$scores->{$seq}}) {
			push @forward, $scores->{$seq}{$cmp}{'forward'} // "0.0";
			push @reverse, $scores->{$seq}{$cmp}{'reverse'} // "0.0";
		}

		$reporter->combine(@forward);
		$output .= $reporter->string() . $newline;
		$reporter->combine(@reverse);
		$output .= $reporter->string() . $newline;
	}

	return $output;
}

sub ascii_output {
	my ($scores) = @_;
	my $reporter = Text::ASCIITable->new({ headingText => "Alignments (forward and reverse)" });
	my $output;
	my $out_string = IO::Scalar->new(\$output);
	my %headers;
	# construct the set of headers first
	foreach my $seq (sort keys %$scores) {
		foreach my $cmp (sort keys %{$scores->{$seq}}) {
			$headers{$cmp}++;
		}
	}

	my @header = ("", sort keys %headers);
	$reporter->setCols(@header);
	$reporter->addRowLine();

	foreach my $seq (sort keys %$scores) {
		my @forward;
		my @reverse;

		push @forward, "$seq (forward)";
		push @reverse, "$seq (reverse)";

		foreach my $cmp (sort keys %headers) {
			push @forward, $scores->{$seq}{$cmp}{'forward'} // "0.0";
			push @reverse, $scores->{$seq}{$cmp}{'reverse'} // "0.0";
		}
		$reporter->addRow(@forward);
		$reporter->addRow(@reverse);
		$reporter->addRowLine();
	}
	print $out_string $reporter;
	return $output;
}
