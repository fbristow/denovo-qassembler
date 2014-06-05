/*
 * File:   MarkovChainAbundance.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef MARKOV_CHAIN_ABUNDANCE_HH
#define MARKOV_CHAIN_ABUNDANCE_HH

#include "Abundance/MarkovAbundance/MarkovAbundance.hh"

#include "Logging/Logging.hh"

class MarkovChainAbundance : public MarkovAbundance {
public:
	/**
	 * Constructor.
	 * @param graph the graph to use to construct paths
	 * @param paths the paths to compute abundances for
	 */
	MarkovChainAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths);

	/**
	 * Compute abundances for the paths supplied.
	 * @return a set of sequences and their relative abundance in the graph
	 */
	boost::unordered_map<std::string, double> computeAbundances();
private:
	/**
	 * Compute the sum of the weights of the edges outgoing from the specified vertex.
	 * @param v the vertex to compute the sum for.
	 * @param g the graph in which v can be found.
	 * @return the sum of edge weights exiting v.
	 */
	double sumOutgoingEdges(SkinnyGraph::Vertex v, boost::shared_ptr<SkinnyGraph> g);
	/**
	 * Compute the sum of the weights of the edges incoming to the specified vertex.
	 * @param v the vertex to compute the sum for.
	 * @param g the graph in which v can be found.
	 * @return the sum of edge weights entering v.
	 */
	double sumIncomingEdges(SkinnyGraph::Vertex v, boost::shared_ptr<SkinnyGraph> g);
};

#endif // MARKOV_CHAIN_ABUNDANCE_HH
