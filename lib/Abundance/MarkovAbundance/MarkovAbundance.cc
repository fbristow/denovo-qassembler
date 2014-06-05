/*
 * File:   MarkovAbundance.cc
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef MARKOV_ABUNDANCE_CC
#define MARKOV_ABUNDANCE_CC

#include <boost/foreach.hpp>

#include "MarkovAbundance.hh"

DECLARE_LOG(logger, "qassembler.MarkovAbundance");

MarkovAbundance::MarkovAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths) : Abundance(graph, paths) {
	this->beginStateTransitionSum = computeBeginStateTransitionSum();
}

std::size_t
MarkovAbundance::getBeginStateTransitionSum() {
	return this->beginStateTransitionSum;
}

std::size_t
MarkovAbundance::computeBeginStateTransitionSum() {
	std::size_t transitionSum = 0;

	BOOST_FOREACH (boost::shared_ptr<SkinnyGraph> g, this->graph->getGraphs()) {
		BOOST_FOREACH (SkinnyGraph::Vertex v, g->getVertexIterators()) {
			if (boost::in_degree (v, *g->graph()) == 0) {
				boost::shared_ptr<SequenceNode> node = g->node(v);
				boost::shared_ptr<Kmer> firstMer = node->getKmer(0);
				transitionSum += firstMer->getCount();
			}
		}
	}

	return transitionSum;
}

#endif // MARKOV_ABUNDANCE_CC
