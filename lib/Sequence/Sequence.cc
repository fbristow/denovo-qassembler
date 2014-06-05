/* 
 * File:   Sequence.cc
 * Author: fbristow
 *
 * Created on June 17, 2012
 */
#ifndef SEQUENCE_CC
#define SEQUENCE_CC

#include "Sequence.hh"
#include "Util/Util.hh"
#include "Exception/InvalidInputException.hh"
#include <boost/algorithm/string.hpp>

Sequence::Sequence() {}

Sequence::Sequence(Sequence *seq) {
	this->sequence = seq->sequence;
	this->name = seq->name;
	this->comment = seq->comment;
	this->qual = seq->qual;
	this->id = seq->id;
	revcom();
}

Sequence::Sequence(std::string sequence, std::string name, std::string comment, std::string qual) {
	this->sequence = boost::to_upper_copy(sequence);
	this->name = name;
	this->comment = comment;
	this->qual = qual;
	this->id = qassembler::hash(name);
	revcom();
}

void
Sequence::revcom() {
	for (int i = sequence.size() - 1; i >= 0; i--) {
		switch (sequence[i]) {
			case 'T': case 't':
				reverse += 'A';
				break;
			case 'G': case 'g':
				reverse += 'C';
				break;
			case 'C': case 'c':
				reverse += 'G';
				break;
			case 'A': case 'a':
				reverse += 'T';
				break;
			case 'Y': case 'y':
				reverse += 'R';
				break;
			case 'R': case 'r':
				reverse += 'Y';
				break;
			case 'S': case 's':
				reverse += 'S';
				break;
			case 'W': case 'w':
				reverse += 'W';
				break;
			case 'M': case 'm':
				reverse += 'K';
				break;
			case 'K': case 'k':
				reverse += 'M';
				break;
			case 'V': case 'v':
				reverse += 'B';
				break;
			case 'H': case 'h':
				reverse += 'D';
				break;
			case 'D': case 'd':
				reverse += 'H';
				break;
			case 'B': case 'b':
				reverse += 'V';
				break;
			case 'N': case 'n':
				reverse += 'N';
				break;
			default:
				throw InvalidInputException ("Non-DNA input detected.");
		}
	}
}

std::string
Sequence::getName() {
	return this->name;
}

void
Sequence::setName(std::string name) {
	this->name = name;
}

std::string
Sequence::getSequence() {
	return this->sequence;
}

void
Sequence::setSequence(std::string sequence) {
	this->sequence = sequence;
	revcom();
}

std::string
Sequence::getReverseComplement() {
	return this->reverse;
}

std::string
Sequence::getQual() {
	return this->qual;
}

void
Sequence::setQual(std::string qual) {
	this->qual = qual;
}

std::size_t
Sequence::getLength() {
	return this->sequence.size();
}

std::string
Sequence::getComment() {
	return this->comment;
}

void
Sequence::setComment(std::string comment) {
	this->comment = comment;
}

std::size_t
Sequence::getID() {
	return this->id;
}

void
Sequence::setID(std::size_t id) {
	this->id = id;
}

#endif // SEQUENCE_CC
