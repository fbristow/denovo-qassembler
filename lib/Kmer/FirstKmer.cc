/* 
 * File:   FirstKmer.cc
 * Author: fbristow
 *
 * Created on January 6, 2012
 */
#ifndef FIRST_KMER_CC
#define FIRST_KMER_CC

#include "FirstKmer.hh"

FirstKmer::FirstKmer(FirstKmer *mer) : Kmer(mer) {
	this->sequence = mer->sequence;
}

FirstKmer::FirstKmer(std::size_t hash, std::string sequence, std::size_t source,
	       std::size_t position, Kmer::Strand strand) :
	Kmer(hash, 'x', source, position, strand) {
	this->sequence = sequence;
}

std::string 
FirstKmer::getSequence() {
	return this->sequence;
}

char
FirstKmer::getBase() {
	return this->sequence[this->sequence.length() - 1];
}

#endif // FIRST_KMER_CC
