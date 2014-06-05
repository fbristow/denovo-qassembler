#!/bin/sh

PROCESSORS=4
BASE_DIR=experiments

perl ./slice.pl --base-dir $BASE_DIR --alignment data/ML2155.aln --variants-max 4 --length-max 2000
perl ./art454-reads.pl --base-dir $BASE_DIR --processors $PROCESSORS
perl ./assemble.pl --base-dir $BASE_DIR --processors $PROCESSORS --assembler-cmd qassembler.bin --kmer-length 57 
perl ./filter.pl --base-dir $BASE_DIR --processors $PROCESSORS --min-length-pct 0.75
perl ./needle.pl --base-dir $BASE_DIR --processors $PROCESSORS
perl ./report.pl --base-dir $BASE_DIR --processors $PROCESSORS --output-format csv
perl ./assignment.pl --base-dir $BASE_DIR
perl ./summary.pl --base-dir $BASE_DIR
