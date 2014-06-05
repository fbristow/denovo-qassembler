/* 
 * File:   FirstKmer.hh
 * Author: fbristow
 *
 * Created on January 9, 2012
 */
#ifndef FIRST_KMER_HH
#define FIRST_KMER_HH

#include "Kmer.hh"

class FirstKmer: public Kmer {
public:
	/**
	 * Copy constructor.
	 * @param copy the source to copy from.
	 */
	FirstKmer(FirstKmer *);
	/**
	 * Constructor specifying data members.
	 * @param hash the hash for this kmer.
	 * @param sequence the complete k-length sequence for this kmer.
	 * @param source the identifier for the read where this kmer came from.
	 * @param position the position in the read where this kmer was generated.
	 * @param strand the direction of the read when this kmer was generated.
	 */
	FirstKmer(std::size_t hash, std::string sequence, std::size_t source,
		       std::size_t position, Kmer::Strand strand);
	/** get the sequence for the first kmer */
	std::string getSequence();
	/** allow the caller to get the most significant base from this sequence */
	char getBase();
private:
	std::string sequence;
};

#endif // FIRST_KMER_HH
