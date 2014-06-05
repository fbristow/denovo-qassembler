/* 
 * File:   SequenceNode.cc
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef SEQUENCE_NODE_CC
#define SEQUENCE_NODE_CC

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include "SequenceNode.hh"

SequenceNode::SequenceNode() {}

SequenceNode::SequenceNode(SequenceNode *copy) {
	this->name = copy->name;
	this->id = copy->id;
	this->kmers = copy->kmers;
}

SequenceNode::SequenceNode(std::size_t id, std::string name) {
	this->name = name;
	this->id = id;
}

SequenceNode::~SequenceNode() {
	this->kmers.clear();
}

void
SequenceNode::merge(boost::shared_ptr<SequenceNode> source) {
	std::vector<boost::shared_ptr<Kmer> > nKmers = source->getKmers();
	this->kmers.insert(kmers.begin(),
			   nKmers.begin(),
			   nKmers.end());
	updateKmerLocations();
}

std::string
SequenceNode::fullSequence() {
	std::string sequence = "";
	
	if (kmers.size() > 0) {
		sequence = kmers[0]->getSequence();
	
		for (std::size_t k = 1; k < kmers.size(); k++) {
			sequence += kmers[k]->getBase();
		}
	}
	
	return sequence;
}

std::string
SequenceNode::sequence() {
	std::string sequence = "";
	BOOST_FOREACH (boost::shared_ptr<Kmer> k, kmers) {
		sequence += k->getBase();
	}
	return sequence;
}

std::string
SequenceNode::getName() {
	return this->name;
}

void
SequenceNode::setName(std::string name) {
	this->name = name;
}

std::size_t
SequenceNode::getId() {
	return this->id;
}

void
SequenceNode::setId(std::size_t id) {
	this->id = id;
}

std::vector<boost::shared_ptr<Kmer> >
SequenceNode::getKmers() {
	return this->kmers;
}

boost::shared_ptr<Kmer>
SequenceNode::getKmer(std::size_t position) {
	return this->kmers[position];
}

void
SequenceNode::setKmers(std::vector<boost::shared_ptr<Kmer> > kmers) {
	this->kmers = kmers;
	updateKmerLocations();
}

void
SequenceNode::addKmer(std::size_t hash, char base, std::size_t source, std::size_t offset, Kmer::Strand strand) {
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>(hash, base, source, offset, strand);
	this->kmers.push_back(k);
	this->kmerLocation[hash] = this->kmers.size() - 1;
}

void
SequenceNode::addKmer(boost::shared_ptr<Kmer> mer) {
	this->kmers.push_back(mer);
	this->kmerLocation[mer->getHash()] = this->kmers.size() - 1;
}

void
SequenceNode::addKmerAt(boost::shared_ptr<Kmer> mer, std::size_t position) {
	this->kmers.insert(this->kmers.begin() + position, mer);
	this->kmerLocation[mer->getHash()] = position;
	updateKmerLocations();
}

void
SequenceNode::addKmerSourceAt(std::size_t position, std::size_t source, std::size_t offset, Kmer::Strand strand) {
	this->kmers[position]->addSource(source, offset, strand);
}

int32_t
SequenceNode::findKmer(std::size_t hash) {
	int position = -1;
	
	//for (int i = 0; i < this->kmers.size(); i++) {
	//	if (this->kmers[i]->getHash() == hash) {
	//		position = i;
	//	}
	//}
	if (kmerLocation.count(hash) > 0) {
		position = kmerLocation[hash];
	}
	
	return position;
}

std::size_t
SequenceNode::kmerCount() {
	return this->kmers.size();
}

void
SequenceNode::updateKmerLocations() {
	for (std::size_t i = 0; i < this->kmers.size(); i++) {
		this->kmerLocation[this->kmers[i]->getHash()] = i;
	}
}

#endif // SEQUENCE_NODE
