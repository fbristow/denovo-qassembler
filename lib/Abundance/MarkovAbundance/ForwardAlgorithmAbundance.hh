/*
 * File:   ForwardAlgorithmAbundance.hh
 * Author: fbristow
 *
 * Created on February 29, 2012
 */
#ifndef FORWARD_ALGORITHM_ABUNDANCE_HH
#define FORWARD_ALGORITHM_ABUNDANCE_HH

#include "Abundance/MarkovAbundance/MarkovAbundance.hh"

#include "Logging/Logging.hh"

class ForwardAlgorithmAbundance : public MarkovAbundance {
public:
	/**
	 * Constructor.
	 * @param graph the graph to use to construct paths
	 * @param paths the paths to compute abundances for
	 */
	ForwardAlgorithmAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths);

	/**
	 * Compute abundances for the paths supplied.
	 * @return a set of sequences and their relative abundance in the graph
	 */
	boost::unordered_map<std::string, double> computeAbundances();
};

#endif // FORWARD_ALGORITHM_ABUNDANCE_HH
