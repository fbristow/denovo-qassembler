/**
 * File:   FastaStream.cc
 * Author: fbristow
 *
 * Created on June 20, 2012
 */

#ifndef FASTA_STREAM_CC
#define FASTA_STREAM_CC

#include "IO/FastaStream.hh"

FastaStream::FastaStream() {}

FastaStream::FastaStream(std::string filename) {
	this->filename = filename;
	handle = gzopen(filename.c_str(), "r");
	kseq_init(handle);
	k_seq = kseq_init(handle);
}

FastaStream::~FastaStream() {
	kseq_destroy(k_seq);
	gzclose(handle);
}


boost::shared_ptr<Sequence>
FastaStream::nextSeq() {
	boost::shared_ptr<Sequence> seq;
	int status;

	status = kseq_read(k_seq);

	if (status >= 0) {
		std::string sequence (k_seq->seq.s, k_seq->seq.l);
		std::string name (k_seq->name.s, k_seq->name.l);
		std::string comment (k_seq->comment.s, k_seq->comment.l);
		std::string qual (k_seq->qual.s, k_seq->qual.l);
		seq = boost::make_shared<Sequence>(sequence, name, comment, qual);
	} 

	return seq;
}

std::size_t
FastaStream::seqCount() {
	std::size_t count = 0;
	int status;

	while ((status = kseq_read(k_seq)) >= 0) {
		count++;
	}

	gzclose(handle);
	handle = gzopen(filename.c_str(), "r");
	kseq_init(handle);
	k_seq = kseq_init(handle);

	return count;
}

#endif // FASTA_STREAM_CC
