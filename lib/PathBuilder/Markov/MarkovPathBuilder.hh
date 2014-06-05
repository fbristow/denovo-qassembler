/*
 * MarkovPathBuilder.hh
 *
 *  Created on: 2012-06-29
 *      Author: fbristow
 */

#ifndef MARKOVPATHBUILDER_HH_
#define MARKOVPATHBUILDER_HH_

#include "PathBuilder/PathBuilder.hh"
#include "Logging/Logging.hh"
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

class MarkovPathBuilder: public PathBuilder {
public:
	/**
	 * Constructor.
	 * @param graph the graph to use to construct paths
	 */
	MarkovPathBuilder(boost::shared_ptr<SkinnyGraph> graph);

	/**
	 * Construct paths from the supplied graph using a markov approach.
	 * @return a set of max-bandwidth paths from the graph.
	 */
	boost::unordered_set<PathBuilder::Path> buildPaths();
private:

	/**
	 * Sum the weights of the edges incoming to a certain vertex
	 * @param v the vertex in question
	 */
	double sumIncomingEdges(SkinnyGraph::Vertex v);

	/**
	 * Sum the weights of the edges leaving a certain vertex
	 * @param v the vertex in question
	 */
	double sumOutgoingEdges(SkinnyGraph::Vertex v);

	boost::uniform_real<> uniform_distribution;
	boost::mt19937 generator;
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rng;
};

typedef std::pair<SkinnyGraph::Edge, double> EdgeWeightPair;
bool compareEdgeWeightPairs(EdgeWeightPair a, EdgeWeightPair b);

#endif /* MARKOVPATHBUILDER_HH_ */
