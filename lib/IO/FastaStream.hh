/**
 * File:   FastaStream.hh
 * Author: fbristow
 *
 * Created on June 20, 2012
 */

#ifndef FASTA_STREAM_HH
#define FASTA_STREAM_HH

#include <zlib.h>
#include "IO/kseq.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Sequence/Sequence.hh"
#include "QAssemblerConfig.h"

#define LOGGER_NAME "qassembler.FastaStream"
#include "Logging/Logging.hh"

KSEQ_INIT(gzFile, gzread)

class FastaStream {
public:
	/**
	 * Constructor, specifying fasta file to read.
	 * @param filename the fasta file to read.
	 */
	FastaStream(std::string filename);

	/**
	 * Destructor. Closes the file handle as required.
	 */
	~FastaStream();

	/**
	 * The next sequence in the fasta file. Returns NULL when
	 * the last records in the fasta file is read.
	 * @return the next fasta record, or NULL at eof.
	 */
	boost::shared_ptr<Sequence> nextSeq();

	/**
	 * How many FASTA records are in this file?
	 * @return the number of FASTA records in the file.
	 */
	std::size_t seqCount();
private:
	/**
	 * Default constructor, shouldn't be called.
	 */
	FastaStream();

	/** a file handle for the fasta file. */
	gzFile handle;

	/** a place for temporarily storing sequence records from kseq. */
	kseq_t *k_seq;

	/** the name of the file we're working with. */
	std::string filename;
};

#endif // FASTA_STREAM_HH
