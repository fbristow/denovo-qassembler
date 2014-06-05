/* 
 * File:   WeightedEdge.cc
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef WEIGHTED_EDGE_CC
#define WEIGHTED_EDGE_CC

#include "WeightedEdge.hh"

WeightedEdge::WeightedEdge() : weight(1), lockedWeight(1) {};

WeightedEdge::WeightedEdge(std::size_t weight) : weight(weight), lockedWeight(weight) {};

std::size_t
WeightedEdge::getWeight() {
	return this->weight;
}

void
WeightedEdge::setWeight(std::size_t weight) {
	this->weight = weight;
}

WeightedEdge
WeightedEdge::operator++(int) {
	this->weight++;
	return *this;
}

void
WeightedEdge::increaseWeight(std::size_t amount) {
	this->weight += amount;
}

void
WeightedEdge::decreaseWeight(std::size_t amount) {
	if (amount < this->weight) {
		this->weight -= amount;
	} else {
		this->weight = 0;
	}
}

void
WeightedEdge::lockWeight() {
	this->lockedWeight = this->weight;
}

void
WeightedEdge::resetWeight() {
	this->weight = this->lockedWeight;
}

bool
WeightedEdge::removed() {
	return this->weight < 1;
}

#endif // WEIGHTED_EDGE_CC
