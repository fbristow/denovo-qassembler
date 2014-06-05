/* 
 * File:   Abundance.cc
 * Author: fbristow
 *
 * Created on February 29, 2012
 */
#ifndef ABUNDANCE_CC
#define ABUNDANCE_CC

#include <boost/foreach.hpp>

#include "Abundance.hh"

Abundance::Abundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths) {
	this->graph = graph;
	this->paths = paths;
}

boost::shared_ptr<HeftyGraph>
Abundance::getGraph() {
	return this->graph;
}

void
Abundance::setGraph(boost::shared_ptr<HeftyGraph> graph) {
	this->graph = graph;
}

boost::unordered_set<std::string>
Abundance::getPaths() {
	return this->paths;
}

void
Abundance::setPaths(boost::unordered_set<std::string> paths) {
	this->paths = paths;
}

#endif // ABUNDANCE_CC
