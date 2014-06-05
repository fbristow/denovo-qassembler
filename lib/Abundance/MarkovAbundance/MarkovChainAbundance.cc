/*
 * File:   MarkovChainAbundance.cc
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef MARKOV_CHAIN_ABUNDANCE_CC
#define MARKOV_CHAIN_ABUNDANCE_CC

#include <boost/foreach.hpp>
#include <cmath>

#include "MarkovChainAbundance.hh"
#include "Exception/InvalidGraphStateException.hh"

DECLARE_LOG(logger, "qassembler.MarkovChainAbundance");

MarkovChainAbundance::MarkovChainAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths) : MarkovAbundance(graph, paths) {}

boost::unordered_map<std::string, double>
MarkovChainAbundance::computeAbundances() {
	TRACE(logger, "Beginning to generate abundances.");
	boost::unordered_map<std::string, double> records;
	uint16_t kmerLength = graph->getKmerLength();

	BOOST_FOREACH (std::string path, paths) {
		TRACE(logger, "Working on path [" << path << "]");
		std::string currentKmer = path.substr(0, kmerLength);
		std::size_t currentHash = qassembler::hash(currentKmer);
		std::string nextKmer = currentKmer;
		std::size_t nextHash = currentHash;
		boost::shared_ptr<SkinnyGraph> currentGraph; SkinnyGraph::Vertex currentVertex;
		boost::shared_ptr<SkinnyGraph> nextGraph; SkinnyGraph::Vertex nextVertex;
		boost::tie(currentGraph, currentVertex) = graph->getGraphAndVertexForHash(currentHash);
		// initial probability is the probability of transitioning from the 'start' state
		// to the current kmer. the probability of that happening is the number of instances
		// of the first kmer in the node where this hash came from compared to the sum of all
		// instances of first kmers.
		// TODO: what happens when the kmer comes from a vertex which has (or had) incoming edges?
		double probability = log(currentGraph->node(currentVertex)->getKmer(0)->getCount()) -
					log(getBeginStateTransitionSum());
		TRACE(logger, "Initial probability: [" << probability << "]");

		for (std::size_t i = kmerLength; i < path.size(); i++) {
			DEBUG(logger, "Current probability: [" << probability << "]");
			nextKmer.erase(0, 1);
			nextKmer.push_back(path[i]);
			nextHash = qassembler::hash(nextKmer);
			TRACE(logger, "Current kmer: [" << currentKmer << "]");
			TRACE(logger, "Next kmer:    [" << nextKmer << "]");
			TRACE(logger, "Current hash: [" << std::hex << currentHash << std::dec << "]");
			TRACE(logger, "Next hash: [" << std::hex << nextHash << std::dec << "]");
			TRACE(logger, "Getting next graph and vertex.");
			boost::tie(nextGraph, nextVertex) = graph->getGraphAndVertexForHash(nextHash);
			TRACE(logger, "Next graph     [" << nextGraph << "]");
			TRACE(logger, "Next vertex    [" << nextVertex << "]");
			TRACE(logger, "Current graph  [" << currentGraph << "]");
			TRACE(logger, "Current vertex [" << currentVertex << "]");

			if (nextGraph->getId() != currentGraph->getId()) {
				throw InvalidGraphStateException("A path can be generated from only one graph.");
			}

			if (nextVertex != currentVertex) {
				// if the two kmers share a vertex, then we don't need to change the probability ---
				// the probability of transitioning between two consecutive kmers in the same vertex is 100%.
				// we only need to update the probability when the two vertices are not the same.
				
				// We can update the probabilities by determing the sum of incoming edges to nextVertex
				// and dividing weight of the shared edge by the total sum, then multiplying by the current
				// probability.
				double outgoingWeight = sumOutgoingEdges(currentVertex, currentGraph);
				SkinnyGraph::Edge sharedEdge;
				boost::tie(sharedEdge, boost::tuples::ignore) = boost::edge(currentVertex, nextVertex, *currentGraph->graph());
				double sharedWeight = currentGraph->edge(sharedEdge)->getWeight();
				double transitionProbability = log(sharedWeight) - log(outgoingWeight);
				TRACE(logger, "Shared weight: [" << sharedWeight << "], total outgoing weight: [" << outgoingWeight << "]");
				probability += transitionProbability;
			}
			currentKmer = nextKmer;
			currentHash = nextHash;
			currentGraph = nextGraph;
			currentVertex = nextVertex;
		}

		DEBUG(logger, "Final probability for path [" << path << "] is [" << probability << "]");
		records[path] = probability;
	}

	return records;
}

double
MarkovChainAbundance::sumIncomingEdges(SkinnyGraph::Vertex v, boost::shared_ptr<SkinnyGraph> g) {
	double sum = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, boost::in_edges(v, *g->graph())) {
		sum += g->edge(e)->getWeight();
	}

	return sum;
}

double
MarkovChainAbundance::sumOutgoingEdges(SkinnyGraph::Vertex v, boost::shared_ptr<SkinnyGraph> g) {
	double sum = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, boost::out_edges(v, *g->graph())) {
		sum += g->edge(e)->getWeight();
	}

	return sum;
}

#endif // MARKOV_CHAIN_ABUNDANCE_CC
