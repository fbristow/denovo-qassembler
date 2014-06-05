#!/bin/sh

perl ./slice.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --alignment data/ML2155.aln 
perl ./readsim-reads.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --processors 14 --read-cmd /home/fbristow/readsim/readsim.pl --distribution linear 
perl ./cluster-assemble.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --assembler-cmd /home/fbristow/qassembler/build/src/denovo-qassembler --assembler-wrapper-cmd /home/fbristow/qassembler/util/assemble.pl
./cluster-assemble.sh
perl ./filter.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --processors 14 --min-length-pct 0.75
perl cluster-needle.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --needle-wrapper-cmd /home/fbristow/qassembler/util/needle.pl
./cluster-needle.sh
perl ./report.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --processors 14 --output-format csv
perl ./assignment.pl --base-dir /home/fbristow/qassembler-experiments/experiments/
perl ./summary.pl --base-dir /home/fbristow/qassembler-experiments/experiments/ --report-location /home/fbristow/qassembler-experiments/experiments/report.summary
