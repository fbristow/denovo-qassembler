/* 
 * File:   KmerLengthException.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef KMER_LENGTH_EXCEPTION_HH
#define KMER_LENGTH_EXCEPTION_HH

#include <stdexcept>

class KmerLengthException : public std::runtime_error {
public:
	explicit KmerLengthException(const std::string& m) : std::runtime_error(m) {}
	virtual ~KmerLengthException() throw() {};
};

#endif // KMER_LENGTH_EXCEPTION_HH

