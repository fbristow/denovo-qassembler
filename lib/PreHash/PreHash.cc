/* 
 * File:   PreHash.cc
 * Author: fbristow
 *
 * Created on January 30, 2012
 */
#ifndef PRE_HASH_CC
#define PRE_HASH_CC

#include "PreHash.hh"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

DECLARE_LOG(logger, "qassembler.PreHash");

PreHash::PreHash() {
	this->kmerLength = 31;
	this->hashes2reads.rehash(50000);
	this->reads2hashes.rehash(50000);
}

PreHash::PreHash(std::size_t kmerLength, std::size_t initialSize) {
	this->kmerLength = kmerLength;
	this->reads2hashes.rehash(initialSize);
	this->hashes2reads.rehash(initialSize);
}

void
PreHash::addRead(boost::shared_ptr<Sequence> read) {
	std::string sequence = read->getSequence();
	std::size_t id = read->getID();
	addRead(sequence, id, Kmer::FORWARD);
	sequence = read->getReverseComplement();
	addRead(sequence, id, Kmer::REVERSE);
}

void
PreHash::addRead(std::string sequence, std::size_t readId, Kmer::Strand direction) {
	std::string upperSequence = boost::to_upper_copy(sequence);
	std::string subseq = upperSequence.substr(0, kmerLength);
	std::size_t hash = qassembler::hash(subseq);
	addKmer(hash, readId, direction, kmerLength - 1);

	for (std::size_t i = kmerLength; i < upperSequence.length(); i++) {
		subseq.erase(0, 1);
		subseq.push_back(upperSequence[i]);
		hash = qassembler::hash(subseq);
		addKmer(hash, readId, direction, i);
	}
}

void
PreHash::addKmer(std::size_t hash, std::size_t readId, Kmer::Strand direction, std::size_t position) {
	TRACE(logger, "Adding hash: [" << std::hex << hash << std::dec << "] from read [" << readId << "]");
	this->hashes2reads[hash].insert(readId);
	this->reads2hashes[readId][direction].push_back(hash);
}

std::vector<std::size_t>
PreHash::getHashes(std::size_t read, Kmer::Strand direction) {
	return this->reads2hashes[read][direction];
}

boost::unordered_set<std::size_t>
PreHash::getReads(std::size_t hash) {
	return this->hashes2reads[hash];
}

std::size_t
PreHash::hashCount(std::size_t hash) {
	return this->hashes2reads[hash].size();
}

std::size_t
PreHash::kmerCount(std::string kmer) {
	std::size_t hash = qassembler::hash(kmer);
	return this->hashes2reads[hash].size();
}

boost::unordered_set<std::size_t>
PreHash::getAllHashes() {
	boost::unordered_set<std::size_t> hashes;
	boost::unordered_map<std::size_t, boost::unordered_set<std::size_t> >::iterator current;
	for (current = this->hashes2reads.begin(); current != this->hashes2reads.end(); current++) {
		hashes.insert(current->first);
	}

	return hashes;
}

#endif // PRE_HASH_CC
