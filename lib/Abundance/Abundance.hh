/* 
 * File:   Abundance.hh
 * Author: fbristow
 *
 * Created on February 29, 2012
 */
#ifndef ABUNDANCE_HH
#define ABUNDANCE_HH

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph/HeftyGraph.hh"
#include "Graph/Node/SequenceNode.hh"

class Abundance {
public:
	/** get the graph that this abundance computer uses to compute abundances */
	boost::shared_ptr<HeftyGraph> getGraph();
	/** set the graph that this abundance computer uses to compute abundances */
	void setGraph(boost::shared_ptr<HeftyGraph> graph);
	/** get the paths that this abundance computer is computing abundances for */
	boost::unordered_set<std::string> getPaths();
	/** set the paths that this abundance computer is computing abundances for */
	void setPaths(boost::unordered_set<std::string> paths);

	/** 
	 * Compute abundances for the paths provided.
	 */
	virtual boost::unordered_map<std::string, double> computeAbundances() = 0;
protected:
	/** Default constructor. */
	Abundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths);

	/** the graph that we'll search for paths in */
	boost::shared_ptr<HeftyGraph> graph;
	/** the paths that we're computing abundances for */
	boost::unordered_set<std::string> paths;
};

#endif // ABUNDANCE_HH
