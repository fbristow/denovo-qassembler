/*
 * File:   ForwardAlgorithmAbundance.cc
 * Author: fbristow
 *
 * Created on February 8, 2012
 */
#ifndef FORWARD_ALGORITHM_ABUNDANCE_CC
#define FORWARD_ALGORITHM_ABUNDANCE_CC

#include <boost/foreach.hpp>

#include "ForwardAlgorithmAbundance.hh"

DECLARE_LOG(logger, "qassembler.ForwardAlgorithmAbundance");

ForwardAlgorithmAbundance::ForwardAlgorithmAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths) : MarkovAbundance(graph, paths) {}

boost::unordered_map<std::string, double>
ForwardAlgorithmAbundance::computeAbundances() {
	boost::unordered_map<std::string, double> records;
	std::size_t kmerLength = graph->getKmerLength();

	BOOST_FOREACH (std::string s, paths) {
		TRACE(logger, "Starting to compute abundance for [" << s << "].");
		std::string firstKmer = s.substr(0, kmerLength);
		std::string currentKmer = firstKmer;
		std::string previousKmer = currentKmer;
		std::size_t hash = qassembler::hash (currentKmer);
		std::size_t firstHash = hash;
		std::size_t prevHash;
		boost::unordered_map<int, boost::unordered_map<std::size_t, double> > m;

		m[kmerLength - 1][hash] = 1.;

		for (std::size_t i = kmerLength; i < s.size(); i++) {
			previousKmer = currentKmer;
			currentKmer.erase(0, 1);
			currentKmer.push_back(s[i]);
			prevHash = hash;
			hash = qassembler::hash (currentKmer);
			TRACE(logger, "Currently computing probability for: [" << std::hex << hash << std::dec << "] (" << i << " of " << s.size() << ").");
			TRACE(logger, "Currently computing probability for: [" << currentKmer << "] with length: [" << currentKmer.size() << "]");
			// get the vertex and graph that this hash belongs to:
			boost::shared_ptr<SkinnyGraph> g; SkinnyGraph::Vertex v;
			boost::tie(g, v) = graph->getGraphAndVertexForHash(hash);
			TRACE(logger, "Hash [" << std::hex << hash << std::dec << "] is in [" << std::hex << v << std::dec << "]");
			TRACE(logger, "Vertex [" << g->node(v)->getName() << "] is in [" << std::hex << g << std::dec << "]");
			// once I know which vertex the hash belongs to, I need to check the location of
			// the hash within that vertex. if the hash is any place except the first kmer, then
			// that means that the hash itself only has one neighbour and the transition probability
			// from it's neighbour to it is 100% (i.e., the two nodes would have been on a linear
			// chain if we didn't compress the graph). If the hash is the very first kmer in the node,
			// then it DOES have multiple incoming edges. In that case, we actually have to figure out
			// the probability of transitioning from each of those nodes to this one.
			if (g->node(v)->findKmer(hash) == 0) {
				TRACE(logger, "Current hash has more than one neighbour, computing probability of arriving in this node.");
				// this k-mer has (by definition) more than one neighbour, need to explore the 
				// probability of arriving at this node from those neighbours.
				std::size_t sum = 0;
				BOOST_FOREACH (SkinnyGraph::Edge incoming, boost::in_edges(v, *g->graph())) {
					SkinnyGraph::Vertex neighbour = boost::source(incoming, *g->graph());
					double f, t, edgeSum = 0;
					BOOST_FOREACH (SkinnyGraph::Edge outgoing, boost::out_edges(neighbour, *g->graph())) {
						edgeSum += g->edge(outgoing)->getWeight();
					}
					t = g->edge(incoming)->getWeight() / edgeSum;
					// the "neighbour" of the current k-mer is the last k-mer in the node.
					std::size_t neighbouringHash = g->node(neighbour)->getKmer(g->node(neighbour)->kmerCount() - 1)->getHash();

					if ((i - 1 == kmerLength - 1) ^ (neighbouringHash == firstHash)) {
						f = 0.;
					} else if (m.count(i - 1) == 0) {
						f = 0.;
					} else if (m[i - 1].count(neighbouringHash) == 0) {
						f = 0.;
					} else {
						f = m[i - 1][neighbouringHash];
					}

					sum += f * t;
				}
				m[i][hash] = sum;
			} else {
				TRACE(logger, "Current has one neighbour, copying probability from that neighbour.");
				// this k-mer is in the middle of a node, so by definition this k-mer has one
				// neighbour (i.e., the k-mer immediately preceding it in the sequence node).
				// The total sum is just equal to the previous element in the matrix:
				m[i][hash] = m[i - 1][prevHash];
			}

			TRACE(logger, "Current probability: [" << m[i][hash] << "]");
		}

		TRACE(logger, "Final probability computed for [" << s << "]: " << m[s.size() - 1][hash] << ".");
		records[s] = m[s.size() - 1][hash];
	}

	return records;
}

#endif // FORWARD_ALGORITHM_ABUNDANCE_CC
