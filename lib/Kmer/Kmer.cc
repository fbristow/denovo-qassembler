/* 
 * File:   Kmer.cc
 * Author: fbristow
 *
 * Created on December 15, 2011
 */
#ifndef KMER_CC
#define KMER_CC

#include "Kmer.hh"

Kmer::Kmer() {}

Kmer::Kmer(Kmer *mer) {
	this->hash = mer->hash;
	this->base = mer->base;
	this->sources = mer->sources;
}

Kmer::Kmer(std::size_t hash, char base, std::size_t source, std::size_t position, Strand strand) {
	this->hash = hash;
	this->base = base;
	this->sources[source] = std::make_pair(position, strand);
}

std::size_t
Kmer::getHash() {
	return this->hash;
}

void
Kmer::setHash(std::size_t hash) {
	this->hash = hash;
}

char
Kmer::getBase() {
	return this->base;
}

void
Kmer::setBase(char base) {
	this->base = base;
}

void
Kmer::addSource(std::size_t source, std::size_t position, Strand strand) {
	this->sources[source] = std::make_pair(position, strand);
}

boost::unordered_map<std::size_t /* readIdentifier */, Kmer::Source /* position */>
Kmer::getSources() {
	return this->sources;
}

void
Kmer::setSources(boost::unordered_map<std::size_t /* readIdentifier */, Kmer::Source /* position */> sources) {
	this->sources = sources;
}

std::size_t
Kmer::getCount() {
	return this->sources.size();
}

Kmer::Source
Kmer::getSource() {
	return this->sources.begin()->second;
}

std::string 
Kmer::getSequence() {
	return std::string(1, this->base);
}

void
Kmer::addTransition(char base) {
	this->transitions[base]++;
}

std::size_t
Kmer::getTransitionCount(char base) {
	return this->transitions[base];
}

#endif // KMER_CC
